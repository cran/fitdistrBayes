# Generalized gamma AFT regression, synthetic data only.
library(fitdistrBayes)
library(survival)
set.seed(20261107)
n <- 120
d <- data.frame(age10 = rnorm(n), group = factor(rep(c("A", "B"), n/2)))
Q <- .7
w <- (log(rgamma(n, shape = Q^-2))-log(Q^-2))/Q
t <- exp(1-.3*d$age10+.5*(d$group=="B")+.8*w)
censor <- rexp(n, rate = .15)
d$time <- pmin(t, censor)
d$event <- as.integer(t <= censor)

reference <- fitAFTBayes(Surv(time, event) ~ age10 + group, data = d,
  prior = "reference", iter = 10000, warmup = 2000, seed = 40)
print(summary(reference))
print(coef(reference, type = "time_ratio"))
print(confint(reference, type = "time_ratio"))
profiles <- data.frame(age10 = 0, group = factor(c("A", "B"), levels = levels(d$group)))
print(predict(reference, profiles, times = c(1, 3, 5)))

# Fixed objective alternative: pi(beta,sigma,Q) proportional to
# 1 / (sigma * (1 + Q^2)). No prior hyperparameters are required.
simple <- fitAFTBayes(Surv(time, event) ~ age10 + group, d,
  prior = "simple", iter = 10000, warmup = 2000, seed = 41)
print(rbind(reference = coef(reference, type = "time_ratio"),
            simple = coef(simple, type = "time_ratio")))
print(reference$diagnostics)
print(simple$diagnostics)
print(reference$quantile_precision)
# Endpoint MCSE / interval width quantifies Monte Carlo accuracy separately
# from mixing diagnostics. Increase iterations for more accurate endpoints.

# Optional post-processing: inspect all warnings and sampling diagnostics.
print(WAIC(reference))
if (requireNamespace("loo", quietly = TRUE)) print(LOOIC(reference))
# DIC(reference) is intentionally unavailable: E|Q| diverges.

# Complete-data formula interface:
d$lifetime <- t
complete <- fitAFTBayes(lifetime ~ age10 + group, d,
  iter = 10000, warmup = 2000, seed = 43)
print(coef(complete))
print(complete$diagnostics)
# plot(reference, pars = c("beta:groupB", "Q"))
