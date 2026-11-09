library(fitdistrBayes)
E <- asNamespace("fitdistrBayes")
assert <- function(x,msg) if (!isTRUE(x)) stop(msg,call.=FALSE)
fails <- function(expr,pattern) {
  message <- tryCatch({force(expr);"NO ERROR"},error=conditionMessage)
  assert(grepl(pattern,message,ignore.case=TRUE),message)
}
quiet <- list(warn_convergence=FALSE)
assert(nrow(fitAFTBayes_models())==9L,"AFT catalogue")

# Independent ordinary-range event and survival formulae.
time <- exp(seq(-5,5,length.out=21)); eta <- .3;sigma <- .8
X <- matrix(1,length(time),1)
for (name in c("weibull","lognormal","loglogistic","loglaplace","logt")) {
  model <- E$.aft_distribution(name,df=if(name=="logt")4 else NULL)
  f <- E$.aft_standard_loglik(c(eta,sigma),log(time),rep(1,length(time)),X,model)
  s <- E$.aft_standard_loglik(c(eta,sigma),log(time),rep(0,length(time)),X,model)
  z <- (log(time)-eta)/sigma
  expected <- switch(name,
    weibull=dweibull(time,1/sigma,exp(eta),log=TRUE),
    lognormal=dlnorm(time,eta,sigma,log=TRUE),
    loglogistic=dlogis(z,log=TRUE)-log(sigma)-log(time),
    loglaplace=-log(2*sigma)-abs(z)-log(time),
    logt=dt(z,4,log=TRUE)-log(sigma)-log(time))
  survival <- switch(name,weibull=pweibull(time,1/sigma,exp(eta),FALSE,TRUE),
    lognormal=plnorm(time,eta,sigma,FALSE,TRUE),
    loglogistic=plogis(z,lower.tail=FALSE,log.p=TRUE),
    loglaplace=log(ifelse(z>=0,exp(-z)/2,1-exp(z)/2)),
    logt=pt(z,4,lower.tail=FALSE,log.p=TRUE))
  assert(max(abs(f-expected))<1e-9,paste(name,"event density"))
  assert(max(abs(s-survival))<1e-9,paste(name,"observed survival"))
}
for (k in c(.05,.5,1,3,100,1e5)) {
  model <- E$.aft_distribution("gamma")
  th <- c(eta,k)
  f <- E$.aft_standard_loglik(th,log(time),rep(1,length(time)),X,model)
  s <- E$.aft_standard_loglik(th,log(time),rep(0,length(time)),X,model)
  expected <- dgamma(time,k,rate=k*exp(-eta),log=TRUE)
  ref <- pgamma(time,k,rate=k*exp(-eta),lower.tail=FALSE,log.p=TRUE)
  assert(max(abs(f-expected)/pmax(1,abs(expected)))<1e-10,"Gamma event kernel")
  assert(max(abs(s-ref)/pmax(1,abs(ref)))<1e-8,"Gamma survival kernel")
  for(pr in c("reference","jeffreys")) {
    prior <- E$.aft_standard_prior(pr,3L,model)
    refprior <- log(trigamma(k)-1/k)/2+if(pr=="jeffreys")1.5*log(k) else 0
    assert(abs(E$.aft_gamma_logprior(k,prior,3L)-refprior)<1e-9,"Gamma Fisher factor")
  }
}

# Regressions for old GG errors near Q=0, original units, and underflow.
for(q in c(-1e-200,-1e-12,0,1e-12,1e-200)) {
  z <- c(-40,-12,0,12,40)
  assert(max(abs(E$.aft_logdensity(z,0,q)-dnorm(z,log=TRUE)))<1e-6,"GG density lognormal limit")
  assert(max(abs(E$.aft_logsurvival(z,q)-pnorm(z,lower.tail=FALSE,log.p=TRUE)))<1e-6,"GG tail lognormal limit")
}
# Tail probability independently integrated from the density after rescaling
# the integration variable; no pgamma CDF is used as its own oracle.
for(q in c(-1e-9,-1e-8,1e-8,1e-9)) {
  w <- 100; d0 <- E$.aft_logdensity(w,0,q)
  integral <- integrate(function(v)exp(E$.aft_logdensity(w+v/w,0,q)-d0),
                        0,Inf,rel.tol=1e-9)$value
  expected <- d0-log(w)+log(integral)
  assert(abs(E$.aft_logsurvival(w,q)-expected)<1e-7,"GG small-Q remote tail quadrature")
}
Z <- cbind(1,c(-2,-1,0,1,2)*1e200)
assert(E$.aft_propriety(c(0,.2,.1,1,.8),rep(1,5),Z)$ok,"Overflow-safe design rank")

# Student-t point observations: the old full-rank/no-total-fit check alone
# accepted this improper reference posterior (three ties, df=1).
td <- data.frame(t=exp(c(0,0,0,1)))
fails(fitAFTBayes(t~1,td,distr="logt",df=1),"not certified")
fails(fitAFTBayes(t~1,td,distr="logt"),"fixed 'df'")
fails(fitAFTBayes(t~1,td,distr="gamma",df=4),"only with")
fails(fitAFTBayes(t~1,td,distr="lognormal",shape=2),"only for gamma")
fails(fitAFTBayes(t~1,td,distr="weibull",prior="simple"),"not registered")

set.seed(408)
d <- data.frame(z=rnorm(55),g=factor(rep(c("a","b"),length.out=55)))
d$t <- exp(.4+.3*d$z-.2*(d$g=="b")+rnorm(55,sd=.7))
d$event <- as.integer(d$t<3);d$time <- pmin(d$t,3)
for (family in c("weibull","lognormal","loglogistic","loglaplace","logt","exponential","rayleigh","gamma")) {
  for(pr in c("reference","jeffreys")) {
    f <- fitAFTBayes(survival::Surv(time,event)~z+g,d,distr=family,prior=pr,
      df=if(family=="logt")4 else NULL,iter=60,warmup=30,chains=2,seed=80,control=quiet)
    assert(all(is.finite(log_lik(f,draws=4))),paste(family,pr,"finite likelihood"))
    pp <- predict(f,d[1:2,],times=c(0,.5,2,Inf),draws=4,summary=FALSE)
    assert(identical(dim(pp),c(4L,2L,4L)) && all(pp[,,1]==1) && all(pp[,,4]==0),"Survival prediction contract")
    assert(all(apply(pp,c(1,2),function(z) all(diff(z)<=0))),"Monotone survival")
    assert(length(coef(f,type="time_ratio"))==2L,"Time-ratio names")
    assert(all(E$.aft_loo_guard(f)$certified),paste(family,"LOO certificate"))
    assert(suppressWarnings(DIC(f))$estimates$available,paste(family,"certified DIC"))
    assert(suppressWarnings(WAIC(f))$estimates$available,paste(family,"WAIC"))
  }
}

# Exact complete-lognormal posterior with arbitrary full-rank design.
for(pr in c("reference","jeffreys")) {
  f <- fitAFTBayes(t~z+g,d,distr="lognormal",prior=pr,iter=4000,warmup=0,chains=2,seed=28,control=quiet)
  l <- lm(log(t)~z+g,d);p <- length(coef(l));n <- nrow(d);a <- if(pr=="reference")1 else p+1
  shape <- (n+a-p-1)/2;rate <- sum(residuals(l)^2)/2
  expected_s2 <- rate/(shape-1)
  assert(f$diagnostics$exact_or_independent,"Exact lognormal engine")
  b <- as.matrix(f$draws[,f$model$parameters[seq_len(p)]])
  se <- sqrt(diag(vcov(l))*expected_s2/(sum(residuals(l)^2)/(n-p))/nrow(b))
  assert(all(abs(colMeans(b)-coef(l))<6*se),"Normal posterior beta means")
  expected <- sqrt(rate/qgamma(c(.975,.5,.025),shape))
  observed <- quantile(f$draws$sigma,c(.025,.5,.975),names=FALSE)
  assert(max(abs(observed-expected))<.04,"Inverse-gamma posterior scale quantiles")
}

# Exact exponential/Rayleigh censored posterior and equivalence of flat priors.
for(family in c("exponential","rayleigh")) {
  f <- fitAFTBayes(survival::Surv(time,event)~1,d,distr=family,iter=3000,warmup=0,chains=2,seed=15,control=quiet)
  g <- fitAFTBayes(survival::Surv(time,event)~1,d,distr=family,prior="jeffreys",iter=3000,warmup=0,chains=2,seed=15,control=quiet)
  scale <- if(family=="exponential")1 else .5
  expected <- -scale*(digamma(sum(d$event))-log(sum(d$time^(1/scale))))
  se <- scale*sqrt(trigamma(sum(d$event))/nrow(f$draws))
  assert(abs(mean(f$draws[["beta:(Intercept)"]])-expected)<6*se,"Exact transformed-gamma beta mean")
  assert(identical(f$draws,g$draws),"Identical fixed-scale prior measures")
}

# Gamma moment boundaries, fixed-shape route, and LOO refusal.
small <- data.frame(t=c(1,2))
f <- fitAFTBayes(t~1,small,distr="gamma",iter=40,warmup=20,chains=2,seed=19,control=quiet)
assert(!f$moment_status$mean_exists[1] && f$moment_status$mean_exists[2],"Gamma reference small-sample moments")
assert(!suppressWarnings(DIC(f))$estimates$available,"Gamma DIC withheld for nonexistent coefficient mean")
assert(!all(E$.aft_loo_guard(f)$certified),"Gamma deletion cannot discard noninterpolation")
fixed <- fitAFTBayes(t~1,small,distr="gamma",shape=2,iter=40,warmup=20,chains=2,seed=19,control=quiet)
assert(length(fixed$model$parameters)==1 && fixed$diagnostics$exact_or_independent,"Known-shape gamma engine")
assert(all(E$.aft_loo_guard(fixed)$certified),"Known-shape gamma deletion is proper")
fails(confint(fixed,type="time_ratio",level=2),"level")
tiny <- fitAFTBayes(t~1,small,distr="gamma",shape=1e-200,
  iter=100,warmup=0,chains=2,seed=12,control=quiet)
assert(all(is.finite(as.matrix(tiny$draws))),"Known-shape log-gamma draws avoid natural-scale underflow")
assert(all(is.finite(log_lik(tiny,draws=3))),"Extreme known-shape gamma likelihood")
large <- data.frame(t=c(1e250,2e250,3e250))
ray <- fitAFTBayes(t~1,large,distr="rayleigh",iter=100,warmup=0,chains=2,seed=12,control=quiet)
assert(all(is.finite(log_lik(ray,draws=3))),"Exact Rayleigh rate computed on log scale")

# Strongly rescaled covariates: original coefficients, unchanged predictors.
du <- d;du$z <- du$z*1e200
f <- fitAFTBayes(t~z+g,d,distr="lognormal",iter=100,warmup=0,chains=2,seed=3,control=quiet)
g <- fitAFTBayes(t~z+g,du,distr="lognormal",iter=100,warmup=0,chains=2,seed=3,control=quiet)
assert(max(abs(predict(f,type="linear",summary=FALSE)-predict(g,type="linear",summary=FALSE)))<1e-10,"Design-unit invariance")
huge <- f;huge$estimates[2] <- -1000
fails(coef(huge,type="time_ratio"),"representable")
# Use QR, not crossprod followed by Cholesky, near a design singularity.
set.seed(391)
near <- data.frame(x=rnorm(60),e=rnorm(60))
near$z <- near$x+1e-8*near$e
near$t <- exp(.2+.4*near$x+rnorm(60))
f <- fitAFTBayes(t~x+z,near,distr="lognormal",iter=100,warmup=0,chains=2,seed=3,control=quiet)
assert(all(is.finite(log_lik(f,draws=3))),"Consistent QR tolerance and stable covariance factor")
near$z <- near$x
fails(fitAFTBayes(t~x+z,near,distr="lognormal"),"full column rank")
cat("AFT families, objective measures, exact oracles, safeguards and integration passed.\n")
