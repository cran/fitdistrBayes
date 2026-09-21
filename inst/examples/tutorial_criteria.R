# fitdistrBayes 0.5.0 -- optional Bayesian model comparison
# Run section by section. Results are printed in the console; no files are saved.
library(fitdistrBayes)
packageVersion("fitdistrBayes")

# 1. The usual call is unchanged: no information criteria are computed.
set.seed(101)
x <- rgamma(100, shape = 2.5, rate = 1.3)
gamma_ref <- fitdistrBayes(x, "gamma", "reference-shape", seed = 102)
gamma_ref
is.null(gamma_ref$criteria)  # TRUE

# 2. Compute WAIC and DIC afterwards, using the existing chains.
# This does not fit the model again and requires no additional package.
criteria(gamma_ref, methods = c("waic", "dic"))

# 3. Ask for criteria in the fitting call. They are calculated AFTER sampling.
gamma_j <- fitdistrBayes(x, "gamma", "jeffreys", seed = 103,
                        criteria = c("waic", "dic"))
gamma_j$criteria

# 4. Fit another distribution to EXACTLY the same observations.
weibull_ref <- fitdistrBayes(x, "weibull", "reference", seed = 104)
compare_models(Gamma_reference = gamma_ref, Gamma_Jeffreys = gamma_j,
               Weibull_reference = weibull_ref, criterion = "waic")

# DIC is intentionally unavailable for Weibull: the posterior mean of scale
# does not exist under these priors. A finite Monte Carlo average cannot fix it.
DIC(weibull_ref)

# 5. PSIS-LOO uses the optional package 'loo'. Install it once if necessary:
# install.packages("loo")
if (requireNamespace("loo", quietly = TRUE)) {
  loo_gamma <- LOOIC(gamma_ref)
  loo_gamma
  loo_gamma$details$looic$diagnostics
  compare_models(Gamma = gamma_ref, Weibull = weibull_ref, criterion = "looic")

  # TRUE requests all three criteria in one post-processing pass.
  gamma_all <- fitdistrBayes(x, "gamma", "jeffreys", seed = 103,
                            criteria = TRUE)
  gamma_all$criteria
  identical(gamma_all$chains, gamma_j$chains)  # TRUE
}

# 6. Independent right censoring: x first, status second.
set.seed(105)
censoring <- rexp(length(x), rate = 0.3)
observed <- pmin(x, censoring)
status <- as.integer(x <= censoring)
table(status)

cens_gamma <- fitcensBayes(observed, status, "gamma", "reference-shape",
                          seed = 106, criteria = c("waic", "dic"))
cens_gamma$criteria
cens_weibull <- fitcensBayes(observed, status, "weibull", "reference", seed = 107)

# The likelihood includes log survival for censored observations, not the
# density of their imputed lifetimes. Only compare these two censored fits.
compare_models(Gamma = cens_gamma, Weibull = cens_weibull, criterion = "waic")
if (requireNamespace("loo", quietly = TRUE)) {
  comparison <- compare_models(Gamma = cens_gamma, Weibull = cens_weibull,
                               criterion = "looic")
  comparison
  comparison$table
}

# Lower IC is better. Also inspect R-hat/ESS, Pareto-k and warnings: a numerical
# ordering alone is not a reliable conclusion. Do not compare the original
# complete-data fit against the censored-data fit: they use different data.
# help("criteria") explains assumptions, paired SEs, and result components.
