library(fitdistrBayes)
local({
  E <- asNamespace("fitdistrBayes")
  must_error <- function(expr,pattern) {
    msg <- tryCatch({ force(expr); "NO ERROR" },error=conditionMessage)
    stopifnot(grepl(pattern,msg,ignore.case=TRUE))
  }
  set.seed(100); state <- .Random.seed
  must_error(fitdistrBayes(1:20,"exponential","reference",criterias=TRUE),"Unused arguments")
  must_error(fitdistrBayes(array(1:8,c(2,2,2)),"exponential","reference"),"numeric vector")
  must_error(fitcensBayes(1:3+0i,c(1,1,1),"exponential","reference"),"numeric vector")
  stopifnot(identical(state,.Random.seed))
  # Stable inverse of the weighted Lindley mean over 400 orders of magnitude.
  for (mu in c(1e-200,1e-160,1,1e160,1e200)) for (phi in c(.2,2,1e100)) {
    lambda <- E$.fdb_weighted_lindley_lambda_from_mean(mu,phi)
    stopifnot(is.finite(lambda),lambda>0)
    log_mean <- log(phi) + E$.fdb_logsumexp(c(log(lambda),log(phi),0)) -
      log(lambda) - E$.fdb_logsumexp(c(log(lambda),log(phi)))
    stopifnot(abs(log_mean-log(mu))<1e-10)
  }
  tr <- E$.cens_transform("weighted lindley",c(lambda=1e200,phi=2))
  u <- tr$to(c(lambda=1e200,phi=2))
  stopifnot(all(is.finite(u)),max(abs(log(tr$from(u)$theta/c(lambda=1e200,phi=2))))<1e-10)
  tr <- E$.cens_transform("beta",c(shape1=1,shape2=1))
  for (u in list(c(40,40),c(-40,40),c(750,750),c(-750,750))) {
    z <- tr$from(u)
    if (abs(u[1])==40) stopifnot(all(is.finite(z$theta)),all(z$theta>0),is.finite(z$jacobian))
  }
  stopifnot(abs(E$.cens_transform("gamma",c(shape=1e200,rate=1e-200))$to(c(shape=1e200,rate=1e-200))[1]-400*log(10))<1e-10)

  # Scale-invariant mean ESS and representable standard deviations/SEs.
  set.seed(200)
  mat <- replicate(4,as.numeric(arima.sim(list(ar=.95),n=1000)))
  base <- E$.fdb_ess_matrix(mat)
  stopifnot(base<1000)
  for (factor in c(1e-200,1e-10,1e200)) {
    stopifnot(abs(E$.fdb_ess_matrix(mat*factor)/base-1)<1e-10,
      abs(E$.fdb_stable_sd(as.numeric(mat)*factor)/(sd(as.numeric(mat))*factor)-1)<1e-10)
  }
  stopifnot(abs(E$.fdb_ess_matrix(mat+1e6)/base-1)<1e-5,
    abs(E$.fdb_ic_se(c(1e200,2e200,3e200))/(sqrt(3)*1e200)-1)<1e-12,
    abs(E$.fdb_ic_se(c(1e-200,2e-200,3e-200))/(sqrt(3)*1e-200)-1)<1e-12)
  frozen <- E$.fdb_summarize(list(matrix(1,40,1,dimnames=list(NULL,"p")),
                                 matrix(1,40,1,dimnames=list(NULL,"p"))),FALSE)
  stopifnot(is.infinite(frozen$rhat),frozen$ess_bulk==0,frozen$ess_mean==0,
    frozen$ess_tail==0,is.na(frozen$mcse_mean))

  a <- fitdistrBayes(1:20,"exponential","reference",iter=500,warmup=250,chains=2,seed=9)
  big <- a; big$.loglik <- function(theta) rep(-1e307,20)
  result <- suppressWarnings(WAIC(big))
  stopifnot(!result$estimates$available,!result$estimates$reliable,is.na(result$estimates$estimate),
    nzchar(result$estimates$reason))
  broken <- a; broken$.loglik <- function(theta) rep(NaN,20)
  result <- suppressWarnings(E$.fdb_criteria_after_fit(broken,c("waic","dic")))
  stopifnot(inherits(result,"fitdistrBayes_criteria"),!any(result$estimates$available),
            length(result$diagnostics$computation_error)==1L,identical(broken$chains,a$chains))

  spec <- fitdistrBayes_model(
    density=function(x,theta,log=FALSE) {
      v <- if(abs(theta)<.5) rep(-Inf,length(x)) else dexp(x,1,log=TRUE)
      if(log) v else exp(v)
    },
    prior=function(theta,log=FALSE) {
      v <- if(abs(theta)>=.5 && abs(theta)<=2) -log(3) else -Inf
      if(log) v else exp(v)
    },start=c(theta=1),lower=-2,upper=2,name="disconnected",engine="custom",independent=TRUE,
    sampler=function(chains,n_save,...) lapply(seq_len(chains),function(i)
      matrix(sample(c(-1,1),n_save,TRUE)*runif(n_save,.5,2),ncol=1,dimnames=list(NULL,"theta"))),
    propriety=TRUE,moments=data.frame(parameter="theta",mean_exists=TRUE,variance_exists=TRUE))
  fit <- suppressWarnings(fitdistrBayes(c(.5,1,2),spec,iter=200,warmup=100,chains=2,seed=3,
                                       criteria=c("waic","dic")))
  stopifnot(inherits(fit,"fitdistrBayes"),identical(fit$criteria$estimates$available,c(TRUE,FALSE)))

  # A sure censoring event has an exactly zero contribution, not a bad Pareto tail.
  if(requireNamespace("loo",quietly=TRUE)) {
    b <- fitcensBayes(c(1:20,0),c(rep(1,20),0),"exponential","reference",
                      iter=500,warmup=250,chains=2,seed=9)
    messages <- character()
    cb <- withCallingHandlers(LOOIC(b),warning=function(w) {
      messages <<- c(messages,conditionMessage(w)); invokeRestart("muffleWarning")
    })
    ca <- LOOIC(a)
    stopifnot(!length(messages),cb$estimates$reliable,
      abs(cb$estimates$estimate-ca$estimates$estimate)<1e-10,
      identical(cb$details$loo_zero_information_observations,21L),
      all(log_lik(b)[,21]==0),
      cb$details$looic$pointwise[21,"looic"]==0,
      cb$details$looic$pointwise[21,"mcse_elpd_loo"]==0,
      identical(dim(cb$details$looic),c(500L,21L)))
    invisible(capture.output(print(cb$details$looic)))
    invisible(loo::loo_compare(cb$details$looic,cb$details$looic))
    # Count models at zero are not sure censoring events and must not be removed.
    discrete <- fitcensBayes(c(1,2,3,0),c(1,1,1,0),"Poisson","jeffreys",
      iter=100,warmup=50,chains=2,seed=8,control=list(warn_convergence=FALSE))
    stopifnot(!length(E$.fdb_zero_information(discrete)))
  }
  cat("PASS: numerical extremes, unit-invariant ESS, frozen chains, optional failure isolation, and zero-information LOO.\n")
})
