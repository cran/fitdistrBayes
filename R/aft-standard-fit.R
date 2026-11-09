.aft_prepare_frame <- function(formula,data,na.action) {
  if (!inherits(formula,"formula") || length(formula)!=3L) stop("Supply a two-sided regression formula.")
  mf <- stats::model.frame(formula,data=data,
    na.action=if (na.action=="fail") stats::na.fail else stats::na.omit)
  tt <- attr(mf,"terms")
  if (attr(tt,"intercept")!=1L) stop("The AFT interface requires an intercept for its propriety checks.")
  if (length(attr(tt,"offset"))) stop("Offsets are not supported by this AFT interface.")
  response <- stats::model.response(mf)
  if (inherits(response,"Surv")) {
    if (!identical(attr(response,"type"),"right")) stop("Only ordinary right censoring is supported; left truncation, interval censoring and multi-state responses are not supported.")
    time <- as.numeric(response[,1L]); status <- as.integer(response[,2L])
  } else {
    if (!is.numeric(response) || is.complex(response) || !is.null(dim(response)))
      stop("Response must be numeric complete times or survival::Surv(time, status).")
    time <- as.numeric(response); status <- rep(1L,length(time))
  }
  if (!length(time) || any(!is.finite(time)) || any(time<0) || anyNA(status) ||
      any(!status %in% 0:1) || any(time[status==1L]<=0))
    stop("AFT times must be finite and nonnegative; failure times must be strictly positive, with status 0/1.")
  X <- stats::model.matrix(tt,mf)
  if (any(!is.finite(X))) stop("The AFT design contains nonfinite values.")
  list(frame=mf,terms=tt,X=X,time=time,status=status)
}

.aft_standard_physical <- function(v,p,model,prob) {
  beta <- v[seq_len(p)]
  if (model$shape_estimated) {
    k <- exp(v[p+1L]); model$shape <- k
    if (!is.finite(k) || k<=0) return(rep(NA_real_,p+1L))
    beta[1L] <- beta[1L]-.aft_error_quantile(prob,model)
    return(c(beta,k))
  }
  if (is.null(model$fixed_scale)) {
    sigma <- exp(v[p+1L])
    beta[1L] <- beta[1L]-sigma*.aft_error_quantile(prob,model)
    return(c(beta,sigma))
  }
  beta
}

.aft_standard_target <- function(v,y,status,X,model,prior,prob) {
  p <- ncol(X); th <- .aft_standard_physical(v,p,model,prob)
  ll <- .aft_standard_loglik(th,y,status,X,model)
  if (any(!is.finite(ll))) return(-Inf)
  lp <- if (model$shape_estimated) .aft_gamma_logprior(th[p+1L],prior,p)+v[p+1L] else
    if (is.null(model$fixed_scale)) (1-prior$power)*v[p+1L] else 0
  value <- sum(ll)+lp
  if (is.finite(value)) value else -Inf
}

.aft_standard_mode <- function(target,initial,Vfallback,control) {
  fn <- function(v) { a <- target(v); if (is.finite(a)) -a else 1e100 }
  candidates <- list(list(par=initial,value=fn(initial),convergence=NA_integer_))
  for (method in c("BFGS","Nelder-Mead")) {
    if (length(initial)==1L && method=="Nelder-Mead") next
    fit <- tryCatch(stats::optim(initial,fn,method=method,
      control=list(maxit=control$optimize_maxit,reltol=1e-10)),error=function(e) NULL)
    if (!is.null(fit) && is.finite(fit$value) && fit$value<1e90) candidates[[length(candidates)+1L]] <- fit
  }
  values <- vapply(candidates,function(z) z$value,0.)
  f <- candidates[[which.min(values)]]
  if (!is.finite(f$value) || f$value>=1e90) stop("No finite AFT posterior starting point was found.")
  H <- tryCatch(stats::optimHess(f$par,fn),error=function(e) NULL)
  good <- !is.null(H) && all(is.finite(H))
  if (good) {
    H <- (H+t(H))/2
    ev <- eigen(H,symmetric=TRUE,only.values=TRUE)$values
    good <- min(ev)>0 && max(ev)/min(ev)<1e10
  }
  # A non-smooth Laplace likelihood can have a singular numerical Hessian.
  # The positive fixed fallback changes geometry only, never the target.
  V <- if (good) chol2inv(chol(H)) else Vfallback
  if (is.null(tryCatch(chol(V),error=function(e)NULL)))
    V <- diag(diag(Vfallback),nrow=length(initial))
  list(center=f$par,V=V,H=if (good) H else NULL,convergence=f$convergence,
    objective=f$value,geometry=if (good) "posterior-mode Hessian" else "fixed design-based fallback",
    starts=lapply(candidates,function(z) z[c("par","value","convergence")]))
}

.aft_standard_slice <- function(target,mode,transform,iter,warmup,thin,chains,control,parameters) {
  d <- length(mode$center); nsave <- floor((iter-warmup)/thin)
  C <- t(chol(mode$V)); D <- diag(sqrt(diag(mode$V)),nrow=d)
  output <- starts <- vector("list",chains); evaluations <- movement <- numeric(chains)
  for (k in seq_len(chains)) {
    counter <- new.env(parent=emptyenv()); counter$n <- 0
    counted <- function(v) {counter$n <- counter$n+1;target(v)}
    initialized <- FALSE
    for (j in seq_len(100L)) {
      v <- mode$center+as.vector(C%*%stats::rnorm(d))*(.5+k/chains)
      if (is.finite(counted(v))) {initialized <- TRUE;break}
    }
    if (!initialized) stop("Could not generate dispersed finite AFT chain starts.")
    starts[[k]] <- v
    out <- matrix(NA_real_,nsave,length(parameters),dimnames=list(NULL,parameters))
    moved <- saved <- 0L
    for (i in seq_len(iter)) {
      previous <- v
      basis <- if (stats::runif(1)<.5) C else D
      for (j in sample.int(d)) {
        direction <- basis[,j]
        one <- function(s) counted(v+s*direction)
        step <- .cens_slice(0,one,control$slice_width,control$slice_steps,control$max_shrink)
        v <- v+step*direction
      }
      if (i>warmup) {
        moved <- moved+as.integer(any(v!=previous))
        if ((i-warmup)%%thin==0L) {
          saved <- saved+1L; th <- transform(v)
          if (any(!is.finite(th))) stop("AFT draw is outside the representable range; no clipping was applied.")
          out[saved,] <- th
        }
      }
    }
    output[[k]] <- out; evaluations[k] <- counter$n;movement[k] <- moved/(iter-warmup)
  }
  list(chains=output,starts=starts,movement=movement,evaluations=evaluations)
}

.aft_loggamma_draw <- function(n, shape) {
  # Gamma(shape)=Gamma(shape+1)*U^(1/shape). Keep this on the log scale
  # when shape<1: the gamma draw itself can underflow while log(draw) is finite.
  if (shape<1) log(stats::rgamma(n,shape=shape+1))+log(stats::runif(n))/shape else
    log(stats::rgamma(n,shape=shape))
}

.aft_design_factor <- function(X) {
  p <- ncol(X); qrX <- qr(X,tol=1e-10)
  if (qrX$rank<p) stop("The computational failure design is numerically rank deficient.")
  R <- qr.R(qrX)[seq_len(p),seq_len(p),drop=FALSE]
  factor <- matrix(0,p,p)
  factor[qrX$pivot,] <- backsolve(R,diag(p))
  factor
}

.aft_standard_exact <- function(y,status,Z,model,prior,scales,iter,warmup,thin,chains,parameters) {
  p <- ncol(Z); ns <- floor((iter-warmup)/thin)
  gaussian <- model$key=="lognormal" && all(status==1L)
  weibull1 <- model$key %in% c("exponential","rayleigh") && p==1L
  gamma1 <- model$key=="gamma" && !model$shape_estimated && p==1L && all(status==1L)
  if (!gaussian && !weibull1 && !gamma1) return(NULL)
  if (gaussian) {
    fit <- stats::lm.fit(Z,y,tol=1e-10)
    rss <- sum(fit$residuals^2)
    shape <- (length(y)+prior$power-p-1)/2
    C <- .aft_design_factor(Z)
  } else {
    scale <- if (gamma1) 1 else model$fixed_scale
    z <- y/scale; maximum <- max(z)
    lograte <- maximum+log(sum(exp(z-maximum)))
    shape <- if (gamma1) length(y)*model$shape else sum(status)
  }
  out <- lapply(seq_len(chains),function(k) {
    if (gaussian) {
      sigma <- sqrt((rss/2)/stats::rgamma(ns,shape=shape))
      beta <- matrix(stats::rnorm(ns*p),ns,p)%*%t(C)
      beta <- sweep(sweep(beta,1L,sigma,"*"),2L,fit$coefficients,"+")
      z <- cbind(sweep(beta,2L,scales,"/"),sigma)
    } else {
      loglambda <- .aft_loggamma_draw(ns,shape)-lograte
      z <- matrix(if (gamma1) log(model$shape)-loglambda else -scale*loglambda,ncol=1L)
    }
    if (any(!is.finite(z))) stop("An exact AFT draw is outside the representable range; no clipping was applied.")
    colnames(z) <- parameters;z
  })
  list(chains=out,starts=list(),movement=rep(NA_real_,chains),evaluations=rep(0,chains),
       algorithm=if (gaussian) "exact normal/inverse-gamma simulation" else "exact transformed-gamma simulation")
}

.aft_fit_standard <- function(formula,data,model,prior,iter,warmup,thin,chains,seed,
                              na.action,control,criteria,call) {
  na.action <- match.arg(na.action,c("fail","omit"));control <- .aft_control(control)
  iter <- .fdb_scalar_count(iter,"iter",20L);warmup <- .fdb_scalar_count(warmup,"warmup",0L)
  thin <- .fdb_scalar_count(thin,"thin");chains <- .fdb_scalar_count(chains,"chains",2L)
  if (warmup>=iter || floor((iter-warmup)/thin)<10L) stop("Keep at least 10 draws per chain after warmup/thinning.")
  if (!is.null(seed) && (!is.numeric(seed) || length(seed)!=1L || !is.finite(seed) ||
      seed<0 || seed>.Machine$integer.max || seed!=floor(seed))) stop("Invalid seed.")
  frame <- .aft_prepare_frame(formula,data,na.action)
  X <- frame$X;p <- ncol(X);x <- frame$time;status <- frame$status;y <- log(x)
  pr <- .aft_standard_prior(prior,p,model)
  check <- .aft_standard_propriety(y,status,X,model,pr)
  if (!check$ok) stop("AFT posterior propriety is not certified: ",check$reason,call.=FALSE)
  if (!identical(criteria,FALSE)) .fdb_criteria_preflight(criteria,control)
  scales <- apply(abs(X),2L,max);scales[scales==0] <- 1
  Z <- sweep(X,2L,scales,"/")
  parameters <- c(paste0("beta:",colnames(X)),if (model$shape_estimated) "shape" else
    if (is.null(model$fixed_scale)) "sigma")
  prob <- max(.05,min(.5,mean(status)/2))
  sampled <- .fdb_with_seed(seed,.aft_standard_exact(y,status,Z,model,pr,scales,
    iter,warmup,thin,chains,parameters))
  independent <- !is.null(sampled)
  initialization <- list(method="analytic posterior",chain_starts=list())
  if (!independent) {
    event <- status==1L
    lm <- stats::lm.fit(Z[event,,drop=FALSE],y[event],tol=1e-10)
    beta <- as.numeric(lm$coefficients)
    spread <- .fdb_stable_sd(lm$residuals)
    if (!is.finite(spread) || spread<=0) spread <- 1
    if (model$shape_estimated) {
      k <- 1/spread^2
      if (!is.finite(k) || k<=0) stop("Gamma shape initialization exceeds the numerical range; inspect time/design scaling.")
      beta[1L] <- beta[1L]-(digamma(k)-log(k))
      initmodel <- model;initmodel$shape <- k
      beta[1L] <- beta[1L]+.aft_error_quantile(prob,initmodel)
      initial <- c(beta,log(k))
    } else if (is.null(model$fixed_scale)) {
      unit_iqr <- .aft_error_quantile(.75,model)-.aft_error_quantile(.25,model)
      scale <- stats::IQR(lm$residuals)/unit_iqr
      if (!is.finite(scale) || scale<=0) scale <- spread
      beta[1L] <- beta[1L]-scale*.aft_error_quantile(.5,model)
      beta[1L] <- beta[1L]+scale*.aft_error_quantile(prob,model)
      initial <- c(beta,log(scale))
    } else initial <- beta
    Vb <- tcrossprod(.aft_design_factor(Z[event,,drop=FALSE])) * max(spread,.1)^2
    d <- length(initial); V <- diag(1,d);V[seq_len(p),seq_len(p)] <- Vb
    target <- function(v) .aft_standard_target(v,y,status,Z,model,pr,prob)
    mode <- .aft_standard_mode(target,initial,V,control)
    transform <- function(v) {
      th <- .aft_standard_physical(v,p,model,prob)
      th[seq_len(p)] <- th[seq_len(p)]/scales
      th
    }
    sampled <- .fdb_with_seed(seed,.aft_standard_slice(target,mode,transform,iter,warmup,
      thin,chains,control,parameters))
    sampled$algorithm <- "quantile-centered mixed-basis slice sampling"
    initialization <- c(mode,list(method="failure-only least squares and scale/shape start",
      chain_starts=sampled$starts,design_scales=scales))
  }
  posterior <- .fdb_summarize(sampled$chains)
  moments <- .aft_standard_moments(parameters,p,sum(status),model,pr)
  posterior$mean_exists <- moments$mean_exists;posterior$variance_exists <- moments$variance_exists
  posterior$moment_note <- moments$note
  posterior$mean[!moments$mean_exists %in% TRUE] <- NA_real_
  posterior$sd[!moments$variance_exists %in% TRUE] <- NA_real_
  posterior$mcse_mean[!(moments$mean_exists %in% TRUE & moments$variance_exists %in% TRUE)] <- NA_real_
  precision <- .aft_quantile_precision(sampled$chains)
  good <- is.finite(posterior$rhat) & posterior$rhat<=control$rhat_threshold &
    is.finite(posterior$ess_bulk) & posterior$ess_bulk>=control$ess_threshold &
    is.finite(posterior$ess_tail) & posterior$ess_tail>=control$ess_threshold
  message <- if (all(good)) "Rank-normalized R-hat and bulk/tail ESS targets met." else
    "AFT sampling diagnostic targets were not met. Inspect trace plots and increase iterations; a returned fit is not proof of convergence."
  mf <- frame$frame
  ans <- list(call=call,formula=formula,terms=frame$terms,design=X,
    contrasts=attr(X,"contrasts"),xlevels=lapply(mf[vapply(mf,function(z) is.factor(z)||is.character(z),TRUE)],
      function(z) if (is.factor(z)) levels(z) else levels(factor(z))),
    model=list(name=model$name,parameters=parameters,n=length(x),p=p,aft_family=model),
    prior=c(pr,list(posterior_propriety=check$reason)),propriety=check,
    data=x,status=status,omitted=as.integer(attr(mf,"na.action")),na.action=attr(mf,"na.action"),
    chains=sampled$chains,draws=.fdb_long_draws(sampled$chains),summary=posterior,
    estimates=stats::setNames(posterior$median,parameters),moment_status=moments,
    quantile_precision=precision,initialization=initialization,
    engine=list(algorithm=sampled$algorithm,iter=iter,warmup=warmup,thin=thin,chains=chains,
      center_probability=prob,updates_per_iteration=if (independent) 0L else length(parameters)),
    diagnostics=list(converged=all(good),exact_or_independent=independent,messages=message,
      movement_fraction=sampled$movement,target_evaluations=sampled$evaluations,
      max_quantile_mcse_ratio=if (all(is.finite(precision$max_mcse_interval_ratio)))
        max(precision$max_mcse_interval_ratio) else NA_real_,
      max_rhat=max(posterior$rhat),min_ess_bulk=min(posterior$ess_bulk),min_ess_tail=min(posterior$ess_tail)),
    control=control,censoring=list(events=sum(status),censored=sum(status==0L),fraction=mean(status==0L)))
  if (control$store_callables) ans$.loglik <- .aft_standard_loglik_factory(X,y,status,model)
  class(ans) <- "fitAFTBayes"
  if (!identical(criteria,FALSE)) ans$criteria <- .fdb_criteria_after_fit(ans,criteria)
  if (!all(good) && control$warn_convergence) warning(message,call.=FALSE)
  ans
}
