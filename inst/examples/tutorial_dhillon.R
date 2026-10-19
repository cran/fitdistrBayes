# Dhillon in fitdistrBayes 0.6.0. Run one block at a time in the console.
# No files are saved, and no packages are installed by this tutorial.
library(fitdistrBayes)
packageVersion("fitdistrBayes")

# 1. Generate complete lifetimes: beta=3, theta=2.
# theta is positive, NOT restricted to (0,1).
set.seed(20261019)
x <- exp((rlogis(100) - log(2)) / 3)
summary(x)

# 2. Group the three registered prior choices for the same data.
fit_J <- fitdistrBayes(x, "dhillon", "jeffreys", seed = 101)
fit_Rbeta <- fitdistrBayes(x, "dhillon", "reference-beta", seed = 101)
fit_Rtheta <- fitdistrBayes(x, "dhillon", "reference-theta", seed = 102)
print(list(Jeffreys = fit_J, Reference_beta = fit_Rbeta,
           Reference_theta = fit_Rtheta))

# Jeffreys and reference-beta are the same prior: identical seeded chains.
identical(fit_J$draws, fit_Rbeta$draws)
fit_Rtheta$diagnostics
fit_Rtheta$initialization
fit_Rtheta$moment_status
confint(fit_Rtheta)
head(fit_Rtheta$draws)

# 3. Independent administrative right censoring at time 1.
observed <- pmin(x, 1)
status <- as.integer(x <= 1)  # 1 = exact event, 0 = right censored
table(status)

cens_J <- fitcensBayes(observed, status, "dhillon", "jeffreys", seed = 201)
cens_Rbeta <- fitcensBayes(observed, status, "dhillon", "reference-beta", seed = 201)
cens_Rtheta <- fitcensBayes(observed, status, "dhillon", "reference-theta", seed = 202)
print(list(Jeffreys = cens_J, Reference_beta = cens_Rbeta,
           Reference_theta = cens_Rtheta))

# 4. The same posterior through data augmentation.
cens_DA <- fitcensBayes(observed, status, "dhillon", "reference-theta",
  method = "augmentation", iter = 8000, warmup = 2000, seed = 203)
print(list(Observed_likelihood = cens_Rtheta, Augmentation = cens_DA))
cens_DA$diagnostics
# Different chains should agree within Monte Carlo uncertainty, not exactly.

# 5. Prediction: quantiles and survival, NOT lifetime mean or SD.
fit_J$predictive_moments
future <- predict(fit_J, draws = 2000, size = 1, seed = 301)
quantile(future, c(.025, .5, .975))
survival <- predict(cens_Rtheta, type = "survival", times = c(0, .5, 1, 2),
  draws = 2000, seed = 302)
apply(survival, 2, quantile, probs = c(.025, .5, .975))
imputed <- predict(cens_Rtheta, type = "impute", draws = 5, seed = 303)
imputed[, seq_len(min(5, ncol(imputed))), drop = FALSE]

# 6. Optional criteria: not computed during chain generation by default.
WAIC(cens_Rtheta)
DIC(cens_Rtheta)  # May be unavailable when a parameter mean is not finite.
if (requireNamespace("loo", quietly = TRUE)) LOOIC(cens_Rtheta)
# Comparison requires the SAME observed data/status and measurement unit.
cens_weibull <- fitcensBayes(observed, status, "weibull", "reference", seed = 204)
compare_models(cens_Rtheta, cens_weibull, criterion = "waic")

# 7. Mathematical safeguards (examples of intended rejection).
try(fitdistrBayes(x, "dhillon", "mdi"))             # Improper posterior
try(fitcensBayes(c(1,2), c(0,0), "dhillon", "jeffreys")) # No events
try(fitdistrBayes(rep(1,10), "dhillon", "reference"))   # Constant sample

# Parameter theta moments depend on the time unit; medians always exist
# when the posterior is proper. Inspect diagnostics and lengthen chains as
# needed. Mathematical validity does not imply numerical convergence.
if (interactive()) plot(cens_Rtheta, type = "trace")
help("Dhillon")
