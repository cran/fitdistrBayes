# Fixed-error AFT families and gamma regression with unknown/known shape.
# Priors are measures in the reported parameters, not restrictions of the
# generalized-gamma prior. See inst/doc/aft-models-theory.txt.
.aft_distribution <- function(distr, df = NULL, shape = NULL) {
  if (!is.character(distr) || length(distr) != 1L || is.na(distr))
    stop("Supply one AFT distribution name.")
  key <- gsub("[ _-]", "", tolower(distr))
  aliases <- c(gg="gg", gengamma="gg", generalizedgamma="gg", weibull="weibull",
    lognormal="lognormal", loggaussian="lognormal", loglogistic="loglogistic",
    loglaplace="loglaplace", logt="logt", logstudentt="logt", t="logt",
    exponential="exponential", rayleigh="rayleigh", gamma="gamma")
  if (!key %in% names(aliases)) stop("Unsupported AFT distribution. See fitAFTBayes_models().")
  key <- unname(aliases[key])
  positive <- function(z, name) {
    if (!is.numeric(z) || is.complex(z) || length(z) != 1L || !is.finite(z) || z <= 0)
      stop("'", name, "' must be one finite positive number.")
  }
  if (key == "logt") {
    if (is.null(df)) stop("Supply fixed 'df' for log-Student-t AFT; degrees of freedom are not estimated.")
    positive(df, "df")
  } else if (!is.null(df)) stop("'df' is used only with log-Student-t AFT.")
  if (!is.null(shape)) {
    if (key != "gamma") stop("'shape' fixes the shape only for gamma AFT.")
    positive(shape, "shape")
  }
  list(key=key, name=switch(key, gg="generalized gamma AFT", logt="log-Student-t AFT",
    lognormal="lognormal AFT", loglogistic="log-logistic AFT", loglaplace="log-Laplace AFT",
    paste(key, "AFT")), df=df, shape=shape,
    fixed_scale=switch(key, exponential=1, rayleigh=.5, gamma=if (!is.null(shape)) 1 else NULL, NULL),
    shape_estimated=key == "gamma" && is.null(shape))
}

fitAFTBayes_models <- function() {
  data.frame(distribution=c("gg","weibull","lognormal","loglogistic","loglaplace",
    "logt","exponential","rayleigh","gamma"),
    priors=c("reference; simple",rep("reference; jeffreys",8)),
    estimated=c("beta, sigma, Q",rep("beta, sigma",4),"beta, sigma","beta","beta",
                "beta, shape (or beta when shape is supplied)"),
    fixed_arguments=c(rep("none",5),"df is required","sigma = 1","sigma = 0.5",
                      "optional shape; NULL estimates it"),
    stringsAsFactors=FALSE)
}

.aft_standard_prior <- function(prior, p, model) {
  if (!is.character(prior) || length(prior) != 1L || is.na(prior)) stop("Specify one AFT prior name.")
  key <- gsub("[ _-]", "", tolower(prior))
  key <- switch(key, reference="reference", independencejeffreys="reference",
                jeffreys="jeffreys", jointjeffreys="jeffreys", "")
  if (!nzchar(key)) stop("This AFT model supports 'reference' and 'jeffreys'; MDI and 'simple' are not registered here.")
  list(key=key, label=key, power=if (key == "reference") 1 else p+1,
       fixed_scale=!is.null(model$fixed_scale))
}

.aft_log1mexp <- function(x) {
  ans <- numeric(length(x)); far <- x < -log(2)
  ans[far] <- log1p(-exp(x[far])); ans[!far] <- log(-expm1(x[!far])); ans
}

.aft_gamma_logprior <- function(k, prior, p) {
  if (!is.finite(k) || k <= 0) return(-Inf)
  if (k > 1e4) {
    inv <- 1/k
    out <- -log(k)-log(2)/2 + log1p(inv/3-inv^3/15+inv^5/21-inv^7/15)/2
  } else if (k < 1e-4) {
    out <- -log(k) + log1p(-k+k*k*trigamma(1+k))/2
  } else out <- log(trigamma(k)-1/k)/2
  out + if (prior$key == "jeffreys") p*log(k)/2 else 0
}

.aft_gamma_logdensity <- function(z, k) {
  norm <- if (k >= 1000) .5*(log(k)-log(2*pi))-1/(12*k)+1/(360*k^3) else
    k*log(k)-lgamma(k)-k
  h <- expm1(z)-z
  small <- is.finite(z) & abs(z)<.001
  v <- z[small]
  h[small] <- v*v*(.5+v*(1/6+v*(1/24+v*(1/120+v/720))))
  penalty <- k*h
  large <- is.finite(z) & z>30
  penalty[large] <- exp(log(k)+z[large])-k*(1+z[large])
  out <- norm-penalty
  out[!is.finite(z)] <- -Inf
  out
}

.aft_error_logdensity <- function(z, model) {
  ans <- switch(model$key,
    weibull=z-exp(z), exponential=z-exp(z), rayleigh=z-exp(z),
    lognormal=stats::dnorm(z,log=TRUE),
    loglogistic=stats::dlogis(z,log=TRUE),
    loglaplace=-log(2)-abs(z),
    logt=stats::dt(z,df=model$df,log=TRUE),
    gamma=.aft_gamma_logdensity(z,model$shape))
  ans[!is.finite(z)] <- -Inf
  ans
}

.aft_error_logsurvival <- function(z, model) {
  switch(model$key,
    weibull=-exp(z), exponential=-exp(z), rayleigh=-exp(z),
    lognormal=stats::pnorm(z,lower.tail=FALSE,log.p=TRUE),
    loglogistic=stats::plogis(z,lower.tail=FALSE,log.p=TRUE),
    loglaplace={ans <- -log(2)-z; lo <- z<0; ans[lo] <- log1p(-exp(z[lo])/2); ans},
    logt=stats::pt(z,df=model$df,lower.tail=FALSE,log.p=TRUE),
    gamma=.aft_logsurvival(z*sqrt(model$shape),1/sqrt(model$shape)))
}

.aft_error_quantile <- function(prob, model) {
  switch(model$key,
    weibull=log(-log1p(-prob)), exponential=log(-log1p(-prob)), rayleigh=log(-log1p(-prob)),
    lognormal=stats::qnorm(prob), loglogistic=stats::qlogis(prob),
    loglaplace=if (prob<.5) log(2*prob) else -log(2*(1-prob)),
    logt=stats::qt(prob,df=model$df),
    gamma={
      k <- model$shape
      if (k > 1e8) {
        z <- stats::qnorm(prob)
        log1p(z/sqrt(k)+(z*z-1)/(3*k))
      } else {
        lg <- (log(prob)+lgamma(1+k))/k
        if (lg < -30) lg-log(k) else log(stats::qgamma(prob,k))-log(k)
      }
    })
}

.aft_standard_loglik <- function(theta, y, status, X, model) {
  p <- ncol(X)
  if (any(!is.finite(theta))) return(rep(-Inf,length(y)))
  if (model$shape_estimated) {
    k <- theta[p+1L]
    if (k<=0) return(rep(-Inf,length(y)))
    model$shape <- k
    sigma <- 1
  } else sigma <- if (is.null(model$fixed_scale)) theta[p+1L] else model$fixed_scale
  if (!is.finite(sigma) || sigma<=0) return(rep(-Inf,length(y)))
  eta <- as.vector(X %*% theta[seq_len(p)])
  if (any(!is.finite(eta))) return(rep(-Inf,length(y)))
  z <- (y-eta)/sigma
  event <- status==1L
  out <- numeric(length(y))
  out[event] <- .aft_error_logdensity(z[event],model)-log(sigma)-y[event]
  out[!event] <- .aft_error_logsurvival(z[!event],model)
  out
}

# Max-absolute scaling precedes sums of squares. This is safe without long double.
.aft_unit_columns <- function(X) {
  scales <- apply(abs(X),2L,max)
  scales[scales==0] <- 1
  Z <- sweep(X,2L,scales,"/")
  norms <- sqrt(colSums(Z*Z)); norms[norms==0] <- 1
  sweep(Z,2L,norms,"/")
}

# A deliberately sufficient Student-t certificate, including repeated data.
# Reserve an invertible p-row design, then extract disjoint noninterpolable
# (p+1)-row blocks. Their residual lower bounds control sigma -> 0.
.aft_student_certificate <- function(y, X, df, power, tol=1e-10) {
  d <- length(y); p <- ncol(X)
  Z <- .aft_unit_columns(X)
  pivot <- qr(t(Z),tol=tol)$pivot[seq_len(p)]
  rest <- setdiff(seq_len(d),pivot)
  yc <- y-mean(y)
  A <- .aft_unit_columns(cbind(Z,yc))
  blocks <- list()
  while (length(rest)>=p+1L) {
    qrA <- qr(t(A[rest,,drop=FALSE]),tol=tol)
    if (qrA$rank<p+1L) break
    chosen <- rest[qrA$pivot[seq_len(p+1L)]]
    if (qr(A[chosen,,drop=FALSE],tol=tol)$rank<p+1L) break
    blocks[[length(blocks)+1L]] <- chosen
    rest <- setdiff(rest,chosen)
  }
  gain <- length(blocks)*(df+1)
  needed <- d+power-p-1
  list(ok=is.finite(gain) && gain>needed,
    reason=if (is.finite(gain) && gain>needed)
      "Disjoint noninterpolable failure blocks certify the small-scale Student-t integral." else
      "The conservative Student-t small-scale certificate failed; use another validated model or obtain a sharper data-specific proof. This is not a proof of impropriety.",
    blocks=blocks, basis=pivot, tail_gain=gain, required_gain=needed,
    sufficient_only=TRUE)
}

.aft_standard_propriety <- function(y,status,X,model,prior) {
  event <- which(status==1L); p <- ncol(X); d <- length(event)
  if (d<p) return(list(ok=FALSE,reason="There are fewer failures than regression coefficients."))
  Z <- .aft_unit_columns(X[event,,drop=FALSE])
  if (qr(Z,tol=1e-10)$rank<p)
    return(list(ok=FALSE,reason="The design matrix of the observed failures must have full column rank."))
  if (!is.null(model$fixed_scale))
    return(list(ok=TRUE,reason="Fixed-scale error density has exponential tails and the failure design has full column rank.",sufficient_only=TRUE))
  out <- .aft_propriety(y,status,X)
  if (!out$ok) return(out)
  if (model$key=="logt") return(.aft_student_certificate(y[event],X[event,,drop=FALSE],model$df,prior$power))
  if (model$shape_estimated)
    out$reason <- paste(out$reason,"The unrestricted gamma shape tails are integrable under the selected prior.")
  out$sufficient_only <- TRUE
  out
}

.aft_standard_moments <- function(parameters,p,d,model,prior) {
  out <- data.frame(parameter=parameters,mean_exists=NA,variance_exists=NA,
    note="Moment existence is not certified.",stringsAsFactors=FALSE)
  if (!is.null(model$fixed_scale)) {
    out$mean_exists <- out$variance_exists <- TRUE
    out$note <- "All polynomial coefficient moments are finite under the fixed-scale failure-design certificate."
  } else if (model$shape_estimated) {
    bound <- if (prior$key=="reference") d-p else d-p/2
    out$mean_exists[seq_len(p)] <- 1<bound
    out$variance_exists[seq_len(p)] <- 2<bound
    out$mean_exists[p+1L] <- out$variance_exists[p+1L] <- TRUE
    out$note[seq_len(p)] <- paste0("Absolute coefficient moments exist exactly for order < ",bound," under this certified gamma design.")
    out$note[p+1L] <- "All positive polynomial shape moments are finite under the noninterpolation condition."
  } else {
    bound <- d+prior$power-p-1
    out$mean_exists <- 1<bound; out$variance_exists <- 2<bound
    out$note <- paste0("Polynomial moments of order < ",bound," are certified; at or above this bound the corresponding absolute moments diverge.")
    if (model$key=="logt") {
      if (model$df<=1) out$mean_exists[seq_len(p)] <- if (1<bound) NA else FALSE
      if (model$df<=2) out$variance_exists[seq_len(p)] <- if (2<bound) NA else FALSE
      out$note[seq_len(p)] <- paste0("Sufficient coefficient moment condition: order < min(df, ",bound,
        "). Above the df bound existence is not certified; order >= ",bound," diverges.")
    }
  }
  out
}

.aft_standard_loglik_factory <- function(X,y,status,model) {
  force(X);force(y);force(status);force(model)
  function(theta) .aft_standard_loglik(theta,y,status,X,model)
}
