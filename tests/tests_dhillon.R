library(fitdistrBayes)
E <- asNamespace("fitdistrBayes")
assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
fails <- function(expr, pattern) {
  msg <- tryCatch({ force(expr); "NO ERROR" }, error = conditionMessage)
  assert(grepl(pattern, msg, ignore.case = TRUE), msg)
}
quiet <- list(warn_convergence = FALSE, init_jitter = 0.35, slice_steps = 100L)

# Distribution kernels independently checked against base R logistic functions.
for (b in c(0.3, 1, 4, 25)) for (t in c(1e-8, 0.5, 1, 5, 1e8)) {
  th <- c(beta = b, theta = t)
  lx <- seq(-20, 20, length.out = 41)
  ref <- dlogis(lx, location = -log(t) / b, scale = 1 / b, log = TRUE) - lx
  assert(max(abs(E$.fdb_dhillon_log_density(exp(lx), th) - ref)) < 1e-10, "Density mismatch")
  surv <- plogis(lx, location = -log(t) / b, scale = 1 / b, lower.tail = FALSE, log.p = TRUE)
  assert(max(abs(E$.cens_survival("dhillon", exp(lx), th, list()) - surv)) < 1e-10, "Survival mismatch")
  assert(identical(E$.cens_survival("dhillon", 0, th, list()), 0), "S(0) must be exactly one")
  pr <- c(0.001, 0.2, 0.5, 0.9, 0.999)
  inverse <- E$.fdb_dhillon_truncated(rep(0, length(pr)), th, log1p(-pr))
  assert(max(abs(E$.fdb_dhillon_survival(inverse, th) - log1p(-pr))) < 1e-10, "Inverse tail mismatch")
  tr <- E$.cens_transform("dhillon", th)
  back <- tr$from(tr$to(th))
  assert(isTRUE(all.equal(unname(back$theta), unname(th), tolerance = 1e-10)), "Transform roundtrip")
  assert(abs(back$jacobian - sum(log(th))) < 1e-10, "Transform Jacobian")
}

# Exact Fisher determinant and the two ordered reference measures.
for (th in list(c(beta = 2, theta = 3), c(beta = 0.5, theta = 0.01))) {
  b <- th[1]; t <- th[2]
  I <- matrix(c((pi^2+3+3*log(t)^2)/(9*b^2), -log(t)/(3*b*t),
                -log(t)/(3*b*t), 1/(3*t^2)), 2, 2)
  assert(abs(det(I) * (b*t)^2 - (pi^2+3)/27) < 1e-12, "Fisher determinant")
  assert(abs(E$.fdb_dhillon_log_prior(th, "jeffreys") + sum(log(th))) < 1e-12, "Jeffreys measure")
  expected <- -sum(log(th)) - .5*log(pi^2+3+3*log(t)^2)
  assert(abs(E$.fdb_dhillon_log_prior(th, "reference-theta") - expected) < 1e-12, "Reverse reference measure")
}

# Exhaustive small samples: independent knot minimization and leave-one-out.
grid <- expand.grid(a = c(0.5, 1, 2), b = c(0.5, 1, 2), c = c(0.5, 1, 2))
for (i in seq_len(nrow(grid))) for (mask in 0:7) {
  x <- as.numeric(grid[i, ]); d <- as.integer(intToBits(mask))[1:3]
  proper <- E$.fdb_dhillon_proper(x, d)
  lx <- log(x); knots <- -unique(lx)
  A0 <- vapply(knots, function(z) sum(abs(z+lx[d==1])) + sum(pmax(z+lx[d==0], 0)), 0.0)
  assert(identical(proper, sum(d) > 0 && min(A0) > 0), "Propriety boundary")
  looref <- vapply(1:3, function(k) E$.fdb_dhillon_proper(x[-k], d[-k]), TRUE)
  assert(identical(E$.fdb_dhillon_loo(x, d==1), looref), "LOO condition")
  if (proper) for (q in c(1, 2)) {
    Aq <- vapply(knots, function(z) sum(abs(z+lx[d==1])) + sum(pmax(z+lx[d==0],0)) - q*z, 0.0)
    expected <- q < length(x) && min(Aq) > 1e-12
    got <- E$.fdb_dhillon_moment(x, d, q)
    if (!is.na(got)) assert(identical(got, expected), "Theta moment criterion")
  }
}
assert(!E$.fdb_dhillon_moment(c(.5,1), c(1,1), 1), "Boundary moment must not be finite")
assert(E$.fdb_dhillon_moment(c(1,2), c(1,0), 1), "Censored theta mean")
assert(!E$.fdb_dhillon_moment(c(1,2), c(1,0), 2), "Censored theta second moment")
assert(identical(E$.fdb_dhillon_moments(c(1,2),c(1,0)),
                 E$.fdb_dhillon_moments(c(1,2,0),c(1,0,0))), "Zero censoring changes moments")

for (pr in c("jeffreys", "reference-beta", "reference-theta")) {
  fails(fitdistrBayes(rep(2, 4), "dhillon", pr), "improper")
  fails(fitcensBayes(c(1,2), c(0,0), "dhillon", pr), "improper")
  fails(fitcensBayes(c(2,1), c(1,0), "dhillon", pr), "improper")
  fails(fitdistrBayes(c(0,1,2), "dhillon", pr), "positive")
}
fails(fitdistrBayes(c(1,2,3), "dhillon", "mdi"), "improper")
fails(fitcensBayes(c(1,2),c(1,0),"dhillon","mdi"), "improper")
fails(fitdistrBayes(c(1,2,3), "dhillon", "reference", fixed = list(beta=2)), "fixed")
fails(fitdistrBayes(c(1,2,3), "dhillon", "reference", start = c(theta=-1)), "positive")
fails(fitdistrBayes(c(1,2,3), "dhillon", "reference", start = c(foo=2)), "named")

x <- exp((qlogis(seq(.03,.97,length.out=30))-log(3))/2)
fit <- function(pr, model = "dhillon") fitdistrBayes(x,model,pr,
  iter=100,warmup=40,chains=2,seed=900,control=quiet)
J <- fit("jeffreys"); R <- fit("reference-beta")
assert(identical(J$draws,R$draws), "Equivalent priors must have identical seeded chains")
assert(identical(J$draws,fit("reference", "log-logistic")$draws), "Dhillon aliases")
T <- fit("reference-theta")
assert(!identical(J$draws,T$draws), "The reverse reference prior must not be replaced by Jeffreys")
unc <- fitcensBayes(x, rep(1L,length(x)), "dhillon", "jeffreys",
  iter=100,warmup=40,chains=2,seed=900,control=quiet)
assert(identical(J$draws,unc$draws), "No-censoring seeded equivalence")
assert(identical(J$predictive_moments$mean_exists,FALSE), "Predictive warning metadata")
assert(all(is.finite(log_lik(J,draws=10))), "Complete pointwise likelihood")

y <- pmin(x, 1); d <- as.integer(x<=1)
for (method in c("direct","augmentation")) for (pr in c("jeffreys","reference-theta")) {
  f <- fitcensBayes(y,d,"dhillon",pr,iter=100,warmup=40,chains=2,
    seed=901,method=method,control=quiet,start=c(beta=2,theta=3))
  assert(all(as.matrix(f$draws[,c("beta","theta")]) > 0), "Parameter support")
  ll <- f$.loglik(c(beta=2,theta=3))
  ref <- ifelse(d==1,dlogis(log(y),-log(3)/2,.5,log=TRUE)-log(y),
                plogis(log(y),-log(3)/2,.5,lower.tail=FALSE,log.p=TRUE))
  assert(max(abs(ll-ref)) < 1e-10, "Observed rather than imputed likelihood")
  im <- predict(f,type="impute",draws=5,seed=10)
  assert(all(sweep(im,2,y[d==0],"-")>0), "Imputation respects censoring")
  s <- predict(f,type="survival",times=c(0,.5,1,2),draws=5,seed=11)
  assert(all(s[,1]==1) && all(apply(s,1,diff)<=0), "Survival monotonicity")
}

minimal <- fitcensBayes(c(1,2),c(1,0),"dhillon","jeffreys",
  iter=100,warmup=40,chains=2,seed=902,control=quiet)
assert(identical(minimal$moment_status$variance_exists,c(TRUE,FALSE)), "Known divergent moment masking")
assert(is.na(minimal$summary$sd[2]), "Divergent SD must not be reported")
if (requireNamespace("loo", quietly = TRUE)) {
  ic <- suppressWarnings(LOOIC(minimal))
  assert(!ic$estimates$available[1], "LOO must not use improper training posteriors")
}
infinite <- fitdistrBayes(c(.5,1),"dhillon","reference-theta",
  iter=100,warmup=40,chains=2,seed=903,control=quiet)
assert(is.na(infinite$summary$mean[2]), "Infinite theta mean under reverse reference")
assert(!suppressWarnings(DIC(infinite))$estimates$available[1], "DIC requires finite means")
# A latent lifetime may overflow even when its logarithm and parameters do not.
th <- c(beta=.001, theta=1)
latent <- E$.fdb_dhillon_truncated_log(1, th, -2)
assert(latent > log(.Machine$double.xmax), "Extreme latent-time fixture")
assert(is.finite(E$.fdb_dhillon_logtime_density(latent,th)), "Log-time conditional remains finite")
fails(E$.fdb_dhillon_truncated(1, th, -2), "floating-point")
assert(is.finite(E$.fdb_dhillon_target(c(log(2),0),-500,c(-501,-499),c(TRUE,TRUE),"jeffreys")),
  "Posterior target must remain finite when log(theta)=1000")
transform <- E$.cens_transform("dhillon",c(beta=1,theta=exp(500)))
fails(transform$from(c(log(2),0)), "floating-point")
# Independent target check, including the centered-coordinate Jacobian.
for(pr in c("jeffreys","reference-theta")) {
  theta <- c(beta=3,theta=.7); lx <- c(-1,.2,1); event <- c(TRUE,FALSE,TRUE)
  center <- .3; u <- c(log(theta[1]),log(theta[2])+theta[1]*center)
  actual <- E$.fdb_dhillon_target(u,center,lx,event,pr)
  expected <- sum(E$.fdb_dhillon_log_density(exp(lx[event]),theta)) +
    sum(E$.fdb_dhillon_survival(exp(lx[!event]),theta)) +
    E$.fdb_dhillon_log_prior(theta,pr) + sum(log(theta)) + sum(lx[event])
  assert(abs(actual-expected)<1e-10,"Log-coordinate posterior density")
}
settings <- expand.grid(beta=c(.5,2,8),theta=c(.1,1,10),n=c(5L,30L,100L),fraction=c(0,.5,.8))
for (i in c(58L,61L)) {
  s <- settings[i,]; set.seed(7000+i)
  y <- exp((rlogis(s$n)-log(s$theta))/s$beta)
  limit <- as.numeric(quantile(y,1-s$fraction))
  f <- fitcensBayes(pmin(y,limit),as.integer(y<=limit),"dhillon","reference-theta",
    method="augmentation",iter=160,warmup=80,chains=2,
    seed=if(i==58L)9178 else 9190,control=list(warn_convergence=FALSE))
  assert(all(is.finite(as.matrix(f$draws[,c("beta","theta")]))), "Heavy-censoring log-latent regression")
}
cat("PASS: Dhillon kernels, priors, transformations, exact boundaries, moments, aliases, sampling, censoring, prediction and criteria.\n")
