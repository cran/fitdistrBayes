library(fitdistrBayes)
library(survival)
E <- asNamespace("fitdistrBayes")
assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
fails <- function(expr, pattern) {
  msg <- tryCatch({force(expr); "NO ERROR"}, error = conditionMessage)
  assert(grepl(pattern, msg, ignore.case = TRUE), msg)
}

# Independent limiting/submodel checks, not comparisons of wrappers to themselves.
x <- exp(seq(-8, 8, length.out = 41)); eta <- .4; sigma <- 1.3
X <- matrix(1, length(x), 1)
normal <- E$.aft_loglik(c(eta, sigma, 0), log(x), rep(1L,length(x)), X)
assert(max(abs(normal-dlnorm(x, eta, sigma, log=TRUE))) < 1e-12, "Lognormal density limit")
assert(max(abs(E$.aft_loglik(c(eta,sigma,1),log(x),rep(1L,length(x)),X)-dweibull(x,1/sigma,exp(eta),log=TRUE))) < 1e-10, "Weibull submodel density")
assert(max(abs(E$.aft_loglik(c(eta,sigma,1),log(x),rep(0L,length(x)),X)-pweibull(x,1/sigma,exp(eta),lower.tail=FALSE,log.p=TRUE))) < 1e-10, "Weibull survival")
for(q in c(-3,-.7,.2,2)) {
  w <- seq(-3,3,length.out=31); a <- q^-2; g <- a*exp(q*w)
  logf <- dgamma(g,shape=a,log=TRUE)+log(abs(q))+log(g)
  assert(max(abs(E$.aft_logdensity(w,0,q)-logf)) < 1e-10, "Gamma transformation density")
  assert(max(abs(E$.aft_logsurvival(w,q)-pgamma(g,a,lower.tail=q<0,log.p=TRUE))) < 1e-10, "Signed Q survival")
}
for(q in c(-1e-8,1e-8)) {
  w<-seq(-3,3,length.out=31)
  assert(max(abs(E$.aft_logdensity(w,0,q)-dnorm(w,log=TRUE))) < 1e-6, "Small Q density")
  assert(max(abs(E$.aft_logsurvival(w,q)-pnorm(w,lower.tail=FALSE,log.p=TRUE))) < 1e-6, "Small Q tail")
}
assert(abs(E$.aft_logomega(0)+log(24)/2)<1e-14,"Reference at zero")
assert(all(is.finite(vapply(c(-1e8,-10,-1,-.1,0,.1,1,10,1e8), E$.aft_logomega, 0.))), "Reference tail stability")

set.seed(98)
d<-data.frame(z=rnorm(60), group=rep(c("a","b"),30))
d$time<-exp(.5+.4*d$z-.3*(d$group=="b")+rnorm(60,sd=.65))
d$status<-as.integer(d$time<=3);d$time<-pmin(d$time,3)
quiet<-list(warn_convergence=FALSE)
fit<-function(...) fitAFTBayes(Surv(time,status)~z+group,data=d,iter=140,warmup=60,chains=2,seed=21,control=quiet,...)
f<-fit()
assert(identical(f$draws,fit()$draws),"Seeded reproducibility")
assert(identical(names(coef(f)),c("(Intercept)","z","groupb")),"Regression names")
assert(all.equal(unname(coef(f,type="time_ratio")),exp(unname(coef(f)[-1])))==TRUE,"Time ratios")
assert(nobs(f)==60L && identical(model.matrix(f),f$design),"Model metadata")
assert(!f$moment_status$mean_exists[f$moment_status$parameter=="Q"] && is.na(f$summary$mean[f$summary$parameter=="Q"]),"Divergent Q moment suppressed")
assert(all(is.finite(log_lik(f,draws=5))),"Observed likelihood")
assert(identical(dim(log_lik(f,draws=5)),c(5L,60L)),"Likelihood dimensions")
pred<-predict(f,newdata=data.frame(z=c(0,1),group=c("a","b")),times=c(0,1,2,Inf),summary=FALSE,draws=10,seed=3)
assert(identical(dim(pred),c(10L,2L,4L)),"Prediction dimensions")
assert(all(pred[,,1]==1)&all(pred[,,4]==0)&all(apply(pred,c(1,2),function(z)all(diff(z)<=0))),"Survival support and monotonicity")
one<-predict(f,newdata=data.frame(z=0,group="a"),times=1,draws=10,seed=3)
assert(nrow(one)==1 && one$estimate>0 && one$estimate<1,"Single new character level is retained")
fails(predict(f,newdata=data.frame(z=0,group="c"),times=1),"new level")
fails(predict(f,times=-1),"nonnegative")
fails(predict(f,newdata=data.frame(z=NA_real_,group="a"),times=1),"missing")
fails(fit(prior="jeffreys"),"improper")
fails(fit(prior_args=list(Q_sd=2)),"Unused")
fails(fit(prior="proper"),"only the objective")
fails(fit(prior="normal"),"only the objective")
fails(fitAFTBayes(Surv(time,status)~z-1,d),"intercept")
fails(fitAFTBayes(Surv(time,status)~z+I(2*z),d),"full column rank")
fails(fitAFTBayes(Surv(rep(0,60),time,status)~z,d),"left truncation")
fails(fitAFTBayes(Surv(time,status)~z+offset(z),d),"Offsets")
dc<-d;dc$status<-0L
fails(fitAFTBayes(Surv(time,status)~z,dc),"failures")
de<-data.frame(time=exp(1:6),z=1:6)
fails(fitAFTBayes(time~z,de),"interpolated")
dm<-d;dm$z[1]<-NA
fails(fitAFTBayes(Surv(time,status)~z,dm),"missing")
fm<-fitAFTBayes(Surv(time,status)~z,dm,na.action="omit",iter=40,warmup=20,chains=2,control=quiet,seed=1)
assert(nobs(fm)==59 && identical(fm$omitted,1L),"Aligned omission")

# No-censoring equivalence, zero-time censoring and intercept-only models.
du<-d;du$status<-1L
fu<-fitAFTBayes(time~z,du,iter=60,warmup=30,chains=2,seed=7,control=quiet)
fs<-fitAFTBayes(Surv(time,status)~z,du,iter=60,warmup=30,chains=2,seed=7,control=quiet)
assert(identical(fu$draws,fs$draws),"Complete/right interface equivalence")
dz<-d;dz$time[1]<-0;dz$status[1]<-0
fz<-fitAFTBayes(Surv(time,status)~1,dz,iter=40,warmup=20,chains=2,seed=3,control=quiet)
assert(all(log_lik(fz,draws=3)[,1]==0),"S(0)=1 likelihood")
assert(length(coef(fz,type="time_ratio"))==0 && nrow(confint(fz,type="time_ratio"))==0,"No slopes in intercept-only model")

# Both objective routes retain their fixed targets and moment guards.
g<-fit(prior="simple")
assert(!identical(g$draws,f$draws),"Distinct objective targets")
assert(is.na(g$summary$mean[g$summary$parameter=="Q"]),"Simple Q mean is suppressed")
assert(!suppressWarnings(DIC(f))$estimates$available,"DIC reference moment guard")
assert(suppressWarnings(WAIC(f))$estimates$available,"AFT WAIC integration")
assert(!suppressWarnings(DIC(g))$estimates$available,"AFT Simple DIC guard")
assert(all(E$.aft_loo_guard(f)$certified),"LOO failures checked")
changed<-f;changed$status[1]<-1L-changed$status[1]
fails(suppressWarnings(compare_models(f,changed,criterion="waic")),"identical retained")
assert(!identical(E$.fdb_criteria_identity(f),E$.fdb_criteria_identity(changed)),"Censoring belongs to criterion identity")

# Regression transformations must preserve training prediction metadata.
fp<-fitAFTBayes(Surv(time,status)~poly(z,2),d,iter=40,warmup=20,chains=2,seed=4,control=quiet)
assert(max(abs(predict(fp,type="linear",summary=FALSE)-predict(fp,newdata=d,type="linear",summary=FALSE)))<1e-12,"poly prediction basis")
assert(max(abs(predict(fp,type="linear",summary=FALSE)[,1,drop=FALSE]-predict(fp,newdata=d[1,,drop=FALSE],type="linear",summary=FALSE)))<1e-12,"poly basis in a single new row")
assert(grepl(paste0("Failures: ",sum(f$status)),paste(capture.output(print(summary(f))),collapse=" "),fixed=TRUE),"Summary event count")
cat("AFT kernel, formula, sampling, prior, moment, prediction and criterion contracts passed.\n")
