library(fitdistrBayes)
if (requireNamespace("loo", quietly = TRUE)) local({
  E <- asNamespace("fitdistrBayes")
  close <- function(a, b, tol = 1e-10) stopifnot(isTRUE(all.equal(a, b,
    tolerance = tol, check.attributes = FALSE)))
  set.seed(471)
  N <- 600L; C <- 4L; S <- N*C; ids <- rep(seq_len(C), each=N)
  ar <- replicate(C, as.numeric(stats::filter(rnorm(N), .95, method="recursive")))
  u <- as.numeric(plogis(ar)); ordinary <- rnorm(S)
  for (scale in c(1e-30, 1e-120, 1e-200)) {
    x <- -scale*(1+u)
    ll <- cbind(constant=rep(-.3,S), tiny=x, ordinary=ordinary)
    stopifnot(length(unique(x))>1L, all(exp(x-max(x))==1))
    eff <- E$.fdb_loo_efficiency(ll,ids,FALSE,1L)
    stopifnot(identical(eff$flat,1:2), identical(eff$constant,1L),
              eff$r_eff[1]==1, eff$r_eff[2]<.5)
    expected <- loo::relative_eff(matrix((max(x)-x)/(max(x)-min(x)),ncol=1),chain_id=ids,cores=1)
    close(eff$r_eff[2],expected)
    close(eff$r_eff[3],loo::relative_eff(matrix(exp(ordinary-max(ordinary)),ncol=1),chain_id=ids,cores=1))
    pw <- E$.fdb_loo_flat_pointwise(ll[,1:2,drop=FALSE],eff$r_eff[1:2])
    close(pw[2,"elpd_loo"]/scale,mean(x/scale))
    close(pw[2,"mcse_elpd_loo"]/scale,sd(x/scale)/sqrt(S*expected))
    stopifnot(pw[2,"elpd_loo"]<0,pw[2,"mcse_elpd_loo"]>0,
              all(pw[,"influence_pareto_k"]==-Inf),pw[1,"mcse_elpd_loo"]==0)
    if(scale>1e-150) close(pw[2,"p_loo"]/scale^2,var(x/scale)*(S-1)/S)
  }

  # Synthetic callbacks below test numerical contracts, not new statistical
  # models or posterior-propriety claims for a constant likelihood.
  f <- fitdistrBayes(1:20,"exponential","reference",iter=500,warmup=0,chains=2,seed=18)
  callback <- f$.loglik
  fake <- f
  fake$.loglik <- function(theta) {
    v <- callback(theta);v[1] <- -.3;v[2] <- -1e-120*(1+theta[[1]]);v
  }
  state <- .Random.seed; before <- f$draws
  result <- suppressWarnings(criteria(fake,c("waic","looic")))
  stopifnot(all(result$estimates$available),identical(state,.Random.seed),identical(before,f$draws),
    identical(result$details$loo_flat_weight_observations,1:2),
    identical(result$details$loo_constant_loglik_observations,1L),
    identical(result$details$loo_numerically_flat_observations,2L),
    !length(result$details$loo_zero_information_observations),
    result$details$looic$pointwise[2,"elpd_loo"]<0,
    result$details$looic$pointwise[2,"mcse_elpd_loo"]>0,
    !any(1:2 %in% result$details$loo_problematic_observations))
  native <- suppressWarnings(loo::loo(log_lik(fake)[,-(1:2)],r_eff=1,cores=1))
  close(result$details$looic$pointwise[-(1:2),],native$pointwise)
  invisible(capture.output(print(result$details$looic)))
  invisible(loo::loo_compare(result$details$looic,result$details$looic))

  fake$.loglik <- function(theta) rep(-.3,20)
  constant <- suppressWarnings(LOOIC(fake))
  stopifnot(constant$estimates$available,loo::is.psis_loo(constant$details$looic),
    identical(dim(constant$details$looic),c(1000L,20L)))
  close(constant$estimates$estimate,12)
  close(constant$estimates$se,0)
  close(constant$estimates$p_eff,0)
  close(loo::mcse_loo(constant$details$looic),0)
  stopifnot(all(loo::pareto_k_values(constant$details$looic)==-Inf))
  invisible(capture.output(print(constant$details$looic)))
  invisible(loo::loo_compare(constant$details$looic,constant$details$looic))

  # Ordinary nondegenerate columns, including genuine high-k problems,
  # must keep the native PSIS values and warning threshold.
  fake$.loglik <- function(theta) {v <- callback(theta);v[1] <- -1e5*theta[[1]]^2;v}
  bad <- suppressWarnings(LOOIC(fake))
  native <- suppressWarnings(loo::loo(log_lik(fake),r_eff=1,cores=1))
  close(bad$details$looic$pointwise,native$pointwise)
  stopifnot(!bad$estimates$reliable,1L %in% bad$details$loo_problematic_observations,
            !length(bad$details$loo_flat_weight_observations))

  # Real positive early censoring, with no synthetic likelihood callback.
  set.seed(180)
  d <- data.frame(x=rnorm(60));d$time <- exp(.4+.2*d$x+.25*rnorm(60));d$event <- 1L
  d <- rbind(d,data.frame(x=0,time=.01,event=0L))
  quiet <- list(warn_convergence=FALSE)
  a <- fitAFTBayes(survival::Surv(time,event)~x,d,distr="lognormal",
    iter=160,warmup=80,chains=2,seed=933,control=quiet)
  b <- suppressWarnings(fitAFTBayes(survival::Surv(time,event)~x,d,distr="lognormal",
    iter=160,warmup=80,chains=2,seed=933,control=quiet,criteria=c("waic","looic")))
  stopifnot(identical(a$draws,b$draws),identical(a$summary,b$summary),
    all(b$criteria$estimates$available),is.null(b$criteria$diagnostics$computation_error),
    61L %in% b$criteria$details$loo_numerically_flat_observations)
  direct <- suppressWarnings(criteria(a,c("waic","looic")))
  stopifnot(identical(direct,b$criteria))
  for(family in fitAFTBayes_models()$distribution) {
    cat("Checking early-censoring LOO:",family,"\n")
    dd <- d;dd$time[nrow(dd)] <- if(family=="lognormal") .01 else 1e-30
    fit <- fitAFTBayes(survival::Surv(time,event)~x,dd,distr=family,
      df=if(family=="logt")8 else NULL,iter=100,warmup=50,chains=2,seed=913,control=quiet)
    state <- .Random.seed;before <- fit$chains
    scores <- suppressWarnings(criteria(fit,c("waic","looic")))
    stopifnot(all(scores$estimates$available),identical(state,.Random.seed),identical(before,fit$chains))
  }
  # Same shared post-processing for complete and censored univariate fits.
  x <- c(seq(.5,3,length.out=40),1e-30);status <- c(rep(1L,40),0L)
  cens <- fitcensBayes(x,status,"weibull","reference",iter=120,warmup=60,chains=2,seed=83)
  stopifnot(all(suppressWarnings(criteria(cens,c("waic","looic")))$estimates$available))
  cat("PASS: stable LOO efficiency, nonzero flat contributions, all-flat API, genuine Pareto warnings and unchanged chains.\n")
})
