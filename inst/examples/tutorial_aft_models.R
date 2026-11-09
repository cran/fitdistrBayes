# Run one model block at a time. All output goes to the console.
# No datasets or results are written to disk.
library(fitdistrBayes)
library(survival)
packageVersion("fitdistrBayes")
fitAFTBayes_models()

set.seed(20261108)
n <- 120
d <- data.frame(age10 = rnorm(n), group = factor(rep(c("A", "B"), n/2)))
eta <- .5 - .3*d$age10 + .4*(d$group == "B")
censor <- rexp(n, rate = .1)

# 1. Weibull: lifetime shape=1/sigma, scale=exp(eta).
t <- rweibull(n, shape = 1.5, scale = exp(eta))
d$time <- pmin(t, censor)
d$event <- as.integer(t <= censor)
reference <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                        distr="weibull", prior="reference", seed=101)
jeffreys <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                       distr="weibull", prior="jeffreys", seed=102)
rbind(reference=coef(reference,type="all"), jeffreys=coef(jeffreys,type="all"))
summary(reference)
reference$diagnostics
reference$quantile_precision
coef(reference, type="time_ratio")
confint(reference, type="time_ratio")
predict(reference, data.frame(age10=0, group=factor(c("A","B"),levels=levels(d$group))),
        times=c(1,3,5))
# plot(reference)

# 2. Lognormal. A numeric response fits COMPLETE data by exact simulation.
t <- rlnorm(n, meanlog=eta, sdlog=.7)
d$lifetime <- t
complete <- fitAFTBayes(lifetime ~ age10 + group, d, distr="lognormal", seed=103)
summary(complete)
complete$engine
# The Surv response below instead fits the censored observations.
d$time <- pmin(t,censor)
d$event <- as.integer(t<=censor)
reference <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                        distr="lognormal", prior="reference", seed=104)
jeffreys <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                       distr="lognormal", prior="jeffreys", seed=105)
rbind(reference=coef(reference,type="all"), jeffreys=coef(jeffreys,type="all"))
reference$diagnostics
jeffreys$diagnostics

# 3. Log-logistic: logistic error in log time.
t <- exp(eta + .6*rlogis(n))
d$time <- pmin(t,censor)
d$event <- as.integer(t<=censor)
reference <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                        distr="loglogistic", prior="reference", seed=106)
jeffreys <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                       distr="loglogistic", prior="jeffreys", seed=107)
rbind(reference=coef(reference,type="all"), jeffreys=coef(jeffreys,type="all"))
reference$diagnostics
jeffreys$diagnostics

# 4. Log-Laplace: signed exponential error in log time.
error <- sample(c(-1,1),n,replace=TRUE)*rexp(n)
t <- exp(eta + .6*error)
d$time <- pmin(t,censor)
d$event <- as.integer(t<=censor)
reference <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                        distr="loglaplace", prior="reference", seed=108)
jeffreys <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                       distr="loglaplace", prior="jeffreys", seed=109)
rbind(reference=coef(reference,type="all"), jeffreys=coef(jeffreys,type="all"))
reference$diagnostics
jeffreys$diagnostics

# 5. Log-Student-t. df=4 is fixed, NOT estimated.
t <- exp(eta + .6*rt(n,df=4))
d$time <- pmin(t,censor)
d$event <- as.integer(t<=censor)
reference <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                        distr="logt", df=4, prior="reference", seed=110)
jeffreys <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                       distr="logt", df=4, prior="jeffreys", seed=111)
rbind(reference=coef(reference,type="all"), jeffreys=coef(jeffreys,type="all"))
reference$propriety
reference$diagnostics
jeffreys$diagnostics

# 6. Exponential. Reference and Jeffreys coincide (flat beta).
t <- rexp(n,rate=exp(-eta))
d$time <- pmin(t,censor)
d$event <- as.integer(t<=censor)
reference <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                        distr="exponential", prior="reference", seed=112)
jeffreys <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                       distr="exponential", prior="jeffreys", seed=112)
rbind(reference=coef(reference), jeffreys=coef(jeffreys))
reference$diagnostics

# 7. Rayleigh: Weibull shape=2, AFT sigma=1/2.
t <- rweibull(n,shape=2,scale=exp(eta))
d$time <- pmin(t,censor)
d$event <- as.integer(t<=censor)
reference <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                        distr="rayleigh", prior="reference", seed=113)
jeffreys <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                       distr="rayleigh", prior="jeffreys", seed=113)
rbind(reference=coef(reference), jeffreys=coef(jeffreys))
reference$diagnostics

# 8. Gamma: conditional mean=exp(eta); shape is estimated by default.
t <- rgamma(n,shape=2,rate=2*exp(-eta))
d$time <- pmin(t,censor)
d$event <- as.integer(t<=censor)
reference <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                        distr="gamma", prior="reference", seed=114)
jeffreys <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
                       distr="gamma", prior="jeffreys", seed=115)
rbind(reference=coef(reference,type="all"), jeffreys=coef(jeffreys,type="all"))
summary(reference)
reference$moment_status
jeffreys$diagnostics
# Only use shape=2 if the shape is KNOWN; otherwise omit the argument.
# known <- fitAFTBayes(Surv(time,event) ~ age10 + group, d,
#                     distr="gamma", shape=2, seed=116)

# Model comparison on the SAME observed gamma-generated dataset.
weibull <- fitAFTBayes(Surv(time,event) ~ age10 + group,d,distr="weibull",seed=117)
lognormal <- fitAFTBayes(Surv(time,event) ~ age10 + group,d,distr="lognormal",seed=118)
WAIC(reference)
DIC(reference) # Available only when all necessary parameter means exist.
compare_models(gamma=reference,weibull=weibull,lognormal=lognormal,criterion="waic")
if (requireNamespace("loo",quietly=TRUE)) LOOIC(reference)
# Inspect convergence, p_waic and Pareto-k warnings before interpreting scores.

# 9. Prentice generalized gamma: retains Reference/Simple priors.
Q <- .7
w <- (log(rgamma(n,shape=Q^-2))-log(Q^-2))/Q
t <- exp(eta+.7*w)
d$time <- pmin(t,censor)
d$event <- as.integer(t<=censor)
reference <- fitAFTBayes(Surv(time,event) ~ age10 + group,d,
                        distr="gg",prior="reference",seed=119)
simple <- fitAFTBayes(Surv(time,event) ~ age10 + group,d,
                     distr="gg",prior="simple",seed=120)
rbind(reference=coef(reference,type="all"),simple=coef(simple,type="all"))
reference$diagnostics
simple$diagnostics
reference$quantile_precision
# For GG, E|Q| diverges: Q means and DIC are intentionally unavailable.
