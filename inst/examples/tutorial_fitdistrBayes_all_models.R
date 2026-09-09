# =============================================================================
# fitdistrBayes: simple tutorial covering all models and priors
# =============================================================================
#
# This file:
#   1. generates a sample;
#   2. fits by maximum likelihood with MASS::fitdistr, when comparable;
#   3. fits the same model under every available objective prior;
#   4. prints everything to the console and saves no analysis results.
#
# In RStudio, run one section at a time with Ctrl+Enter, or run the entire file:
# source("tutorial_fitdistrBayes_all_models.R")
#
# If the package is not installed, run this ONCE with the tarball path:
# install.packages("fitdistrBayes_0.2.3.tar.gz", repos = NULL, type = "source")

if (!requireNamespace("fitdistrBayes", quietly = TRUE)) {
  stop("Install fitdistrBayes first.")
}
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("The recommended MASS package is not available in this R installation.")
}

library(fitdistrBayes)

# TRUE: quick demonstration. FALSE: settings intended for examining the chains
# more carefully. Even in the full mode, always inspect the diagnostics.
quick_mode <- TRUE
if (quick_mode) {
  tutorial_iter <- 900L
  tutorial_warmup <- 300L
  tutorial_chains <- 2L
} else {
  tutorial_iter <- 4000L
  tutorial_warmup <- 2000L
  tutorial_chains <- 4L
}

# Avoid repeated warnings in quick mode. The final table explicitly identifies
# fits that require more iterations.
tutorial_control <- list(warn_convergence = FALSE)

summary_columns <- c(
  "parameter", "mean", "median", "q2.5", "q97.5",
  "rhat", "ess_bulk", "ess_tail"
)

# An NA in mean does not indicate failure: the package does not report a mean
# when that posterior moment has not been certified to be finite. The median
# and quantiles remain available.

compact_summary <- function(fit) {
  fit$summary[, summary_columns, drop = FALSE]
}

show_bayes <- function(fits) {
  print(lapply(fits, compact_summary))
  invisible(fits)
}

cat("\nCONFIGURATION\n")
cat("iter =", tutorial_iter, "| warmup =", tutorial_warmup,
    "| chains =", tutorial_chains, "\n")

# =============================================================================
# 1. BETA
# =============================================================================
cat("\n\n==================== BETA ====================\n")
set.seed(1001)
x_beta <- rbeta(80, shape1 = 2.5, shape2 = 5)

cat("\nMaximum likelihood:\n")
beta_mean <- mean(x_beta)
beta_precision <- beta_mean * (1 - beta_mean) / var(x_beta) - 1
print(suppressWarnings(MASS::fitdistr(
  x_beta, "beta",
  start = list(
    shape1 = beta_mean * beta_precision,
    shape2 = (1 - beta_mean) * beta_precision
  )
)))

beta_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_beta, "beta", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1002,
    control = tutorial_control
  ),
  Reference = fitdistrBayes(
    x_beta, "beta", "reference", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1003,
    control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(beta_bayes)

# =============================================================================
# 2. CAUCHY
# =============================================================================
cat("\n\n==================== CAUCHY ==================\n")
set.seed(1101)
x_cauchy <- rcauchy(100, location = 1, scale = 2)

cat("\nMaximum likelihood:\n")
print(suppressWarnings(MASS::fitdistr(x_cauchy, "cauchy")))

cauchy_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_cauchy, "cauchy", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1102,
    control = tutorial_control
  ),
  Reference = fitdistrBayes(
    x_cauchy, "cauchy", "reference", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1103,
    control = tutorial_control
  ),
  MDI = fitdistrBayes(
    x_cauchy, "cauchy", "mdi", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1104,
    control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(cauchy_bayes)

# =============================================================================
# 3. CHI-SQUARED
# =============================================================================
cat("\n\n=================== CHI-SQUARED ==============\n")
set.seed(1201)
x_chisq <- rchisq(80, df = 6)

cat("\nMaximum likelihood:\n")
print(suppressWarnings(MASS::fitdistr(
  x_chisq, "chi-squared", start = list(df = mean(x_chisq))
)))

chisq_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_chisq, "chi-squared", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1202,
    control = tutorial_control
  ),
  Reference = fitdistrBayes(
    x_chisq, "chi-squared", "reference", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1203,
    control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(chisq_bayes)

# =============================================================================
# 4. EXPONENTIAL
# =============================================================================
cat("\n\n================= EXPONENTIAL =================\n")
set.seed(1301)
x_exponential <- rexp(80, rate = 1.5)

cat("\nMaximum likelihood:\n")
print(suppressWarnings(MASS::fitdistr(x_exponential, "exponential")))

exponential_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_exponential, "exponential", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1302
  ),
  Reference = fitdistrBayes(
    x_exponential, "exponential", "reference", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1303
  ),
  MDI = fitdistrBayes(
    x_exponential, "exponential", "mdi", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1304
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(exponential_bayes)

# =============================================================================
# 5. GAMMA
# =============================================================================
cat("\n\n===================== GAMMA ===================\n")
set.seed(1401)
x_gamma <- rgamma(80, shape = 2.5, rate = 1.3)

cat("\nMaximum likelihood:\n")
print(suppressWarnings(MASS::fitdistr(x_gamma, "gamma")))

gamma_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_gamma, "gamma", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1402,
    control = tutorial_control
  ),
  First_rule = fitdistrBayes(
    x_gamma, "gamma", "first-rule", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1403,
    control = tutorial_control
  ),
  Reference_shape = fitdistrBayes(
    x_gamma, "gamma", "reference-shape", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1404,
    control = tutorial_control
  ),
  Reference_rate = fitdistrBayes(
    x_gamma, "gamma", "reference-rate", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1405,
    control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(gamma_bayes)

# =============================================================================
# 6. GEOMETRIC: number of failures before the first success
# =============================================================================
cat("\n\n=================== GEOMETRIC =================\n")
set.seed(1501)
x_geometric <- rgeom(80, prob = 0.35)

cat("\nMaximum likelihood:\n")
print(suppressWarnings(MASS::fitdistr(x_geometric, "geometric")))

geometric_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_geometric, "geometric", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1502
  ),
  Reference = fitdistrBayes(
    x_geometric, "geometric", "reference", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1503
  ),
  MDI = fitdistrBayes(
    x_geometric, "geometric", "mdi", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1504,
    control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(geometric_bayes)

# =============================================================================
# 7. LOGNORMAL
# =============================================================================
cat("\n\n=================== LOGNORMAL =================\n")
set.seed(1601)
x_lognormal <- rlnorm(80, meanlog = 1, sdlog = 0.6)

cat("\nMaximum likelihood:\n")
print(suppressWarnings(MASS::fitdistr(x_lognormal, "lognormal")))

lognormal_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_lognormal, "lognormal", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1602
  ),
  Reference = fitdistrBayes(
    x_lognormal, "lognormal", "reference", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1603
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(lognormal_bayes)

# =============================================================================
# 8. LOGISTIC
# =============================================================================
cat("\n\n==================== LOGISTIC ==================\n")
set.seed(1701)
x_logistic <- rlogis(100, location = 1, scale = 2)

cat("\nMaximum likelihood:\n")
print(suppressWarnings(MASS::fitdistr(x_logistic, "logistic")))

logistic_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_logistic, "logistic", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1702,
    control = tutorial_control
  ),
  Reference = fitdistrBayes(
    x_logistic, "logistic", "reference", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1703,
    control = tutorial_control
  ),
  MDI = fitdistrBayes(
    x_logistic, "logistic", "mdi", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1704,
    control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(logistic_bayes)

# =============================================================================
# 9. NEGATIVE BINOMIAL
# =============================================================================
cat("\n\n================ NEGATIVE BINOMIAL =============\n")
set.seed(1801)
known_size <- 5
x_negbin <- rnbinom(100, size = known_size, mu = 8)

cat("\nMaximum likelihood (MASS estimates size and mu):\n")
print(suppressWarnings(MASS::fitdistr(x_negbin, "negative binomial")))
cat("In fitdistrBayes, size is known and only mu is estimated.\n")

negbin_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_negbin, "negative binomial", "jeffreys",
    fixed = list(size = known_size), iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1802
  ),
  Reference = fitdistrBayes(
    x_negbin, "negative binomial", "reference",
    fixed = list(size = known_size), iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1803
  ),
  MDI = fitdistrBayes(
    x_negbin, "negative binomial", "mdi",
    fixed = list(size = known_size), iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1804,
    control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(negbin_bayes)

# =============================================================================
# 10. NORMAL
# =============================================================================
cat("\n\n===================== NORMAL ==================\n")
set.seed(1901)
x_normal <- rnorm(100, mean = 3, sd = 2)

cat("\nMaximum likelihood:\n")
print(suppressWarnings(MASS::fitdistr(x_normal, "normal")))

normal_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_normal, "normal", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1902
  ),
  Reference = fitdistrBayes(
    x_normal, "normal", "reference", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1903
  ),
  MDI = fitdistrBayes(
    x_normal, "normal", "mdi", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 1904
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(normal_bayes)

# =============================================================================
# 11. POISSON
# =============================================================================
cat("\n\n===================== POISSON =================\n")
set.seed(2001)
x_poisson <- rpois(100, lambda = 4)

cat("\nMaximum likelihood:\n")
print(suppressWarnings(MASS::fitdistr(x_poisson, "Poisson")))

poisson_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_poisson, "Poisson", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2002
  ),
  Reference = fitdistrBayes(
    x_poisson, "Poisson", "reference", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2003
  ),
  MDI = fitdistrBayes(
    x_poisson, "Poisson", "mdi", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2004,
    control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(poisson_bayes)

# =============================================================================
# 12. STUDENT-t
# =============================================================================
cat("\n\n===================== STUDENT-t ===============\n")
set.seed(2101)
true_df <- 7
x_t <- 1 + 2 * rt(100, df = true_df)

cat("\nMaximum likelihood (location, scale, and df estimated):\n")
print(suppressWarnings(MASS::fitdistr(
  x_t, "t", start = list(m = median(x_t), s = sd(x_t), df = 5)
)))

cat("\nIn the first three Bayesian routes, df is known and fixed at 7.\n")
t_bayes <- list(
  Jeffreys_fixed_df = fitdistrBayes(
    x_t, "t", "jeffreys", fixed = list(df = true_df),
    iter = tutorial_iter, warmup = tutorial_warmup,
    chains = tutorial_chains, seed = 2102, control = tutorial_control
  ),
  Reference_fixed_df = fitdistrBayes(
    x_t, "t", "reference", fixed = list(df = true_df),
    iter = tutorial_iter, warmup = tutorial_warmup,
    chains = tutorial_chains, seed = 2103, control = tutorial_control
  ),
  MDI_fixed_df = fitdistrBayes(
    x_t, "t", "mdi", fixed = list(df = true_df),
    iter = tutorial_iter, warmup = tutorial_warmup,
    chains = tutorial_chains, seed = 2104, control = tutorial_control
  ),
  Independence_Jeffreys_estimated_df = fitdistrBayes(
    x_t, "t", "independence-jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2105,
    control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(t_bayes)

# =============================================================================
# 13. WEIBULL
# =============================================================================
cat("\n\n===================== WEIBULL =================\n")
set.seed(2201)
x_weibull <- rweibull(100, shape = 1.6, scale = 2)

cat("\nMaximum likelihood:\n")
print(suppressWarnings(MASS::fitdistr(x_weibull, "weibull")))

weibull_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_weibull, "weibull", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2202,
    control = tutorial_control
  ),
  Reference = fitdistrBayes(
    x_weibull, "weibull", "reference", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2203,
    control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(weibull_bayes)

# =============================================================================
# 14. GUMBEL
# =============================================================================
cat("\n\n====================== GUMBEL =================\n")
set.seed(2301)
x_gumbel <- 1 - 2 * log(-log(runif(100)))
cat("MASS::fitdistr has no built-in named route for Gumbel.\n")

gumbel_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_gumbel, "gumbel", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2302,
    control = tutorial_control
  ),
  Reference = fitdistrBayes(
    x_gumbel, "gumbel", "reference", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2303,
    control = tutorial_control
  ),
  MDI = fitdistrBayes(
    x_gumbel, "gumbel", "mdi", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2304,
    control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(gumbel_bayes)

# =============================================================================
# 15. FRECHET: F(x) = exp(-scale * x^(-shape))
# =============================================================================
cat("\n\n===================== FRECHET =================\n")
set.seed(2401)
x_frechet <- (4 / rexp(100))^(1 / 2.5)
cat("MASS::fitdistr has no built-in named route for Frechet.\n")

frechet_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_frechet, "frechet", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2402,
    control = tutorial_control
  ),
  Reference = fitdistrBayes(
    x_frechet, "frechet", "reference", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2403,
    control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(frechet_bayes)

# =============================================================================
# 16. LOMAX
# =============================================================================
cat("\n\n====================== LOMAX ==================\n")
set.seed(2501)
shape_lomax <- 3
scale_lomax <- 2
x_lomax <- scale_lomax * expm1(-log(runif(100)) / shape_lomax)
cat("MASS::fitdistr has no built-in named route for Lomax.\n")

lomax_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_lomax, "lomax", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2502,
    control = tutorial_control
  )
)
cat("\nObjective Bayes: the known reference prior yields an improper posterior.\n")
show_bayes(lomax_bayes)

# =============================================================================
# 17. NAKAGAMI-m: spread = E(X^2)
# =============================================================================
cat("\n\n=================== NAKAGAMI-m ================\n")
set.seed(2601)
shape_nakagami <- 2.5
spread_nakagami <- 4
x_nakagami <- sqrt(rgamma(
  100, shape = shape_nakagami,
  rate = shape_nakagami / spread_nakagami
))
cat("MASS::fitdistr has no built-in named route for Nakagami-m.\n")

nakagami_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_nakagami, "nakagami-m", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2602,
    control = tutorial_control
  ),
  Reference = fitdistrBayes(
    x_nakagami, "nakagami-m", "reference", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2603,
    control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(nakagami_bayes)

# =============================================================================
# 18. EXPONENTIAL-LOGARITHMIC
# =============================================================================
cat("\n\n============= EXPONENTIAL-LOGARITHMIC ==========\n")
set.seed(2701)
theta_el <- 0.4
rate_el <- 1.3
u_el <- runif(100)
x_el <- -(
  log(-expm1((1 - u_el) * log(theta_el))) - log1p(-theta_el)
) / rate_el
cat("MASS::fitdistr has no built-in named route for this model.\n")

el_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_el, "exponential-logarithmic", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2702,
    control = tutorial_control
  ),
  MDI = fitdistrBayes(
    x_el, "exponential-logarithmic", "mdi", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2703,
    control = tutorial_control
  ),
  Reference_theta = fitdistrBayes(
    x_el, "exponential-logarithmic", "reference-theta",
    iter = tutorial_iter, warmup = tutorial_warmup,
    chains = tutorial_chains, seed = 2704, control = tutorial_control
  ),
  Reference_rate = fitdistrBayes(
    x_el, "exponential-logarithmic", "reference-rate",
    iter = tutorial_iter, warmup = tutorial_warmup,
    chains = tutorial_chains, seed = 2705, control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(el_bayes)

# =============================================================================
# 19. RICIAN
# =============================================================================
cat("\n\n====================== RICIAN =================\n")
set.seed(2801)
noncentrality_rician <- 5
scale_rician <- 2
x_rician <- sqrt(
  rnorm(100, noncentrality_rician, scale_rician)^2 +
    rnorm(100, 0, scale_rician)^2
)
cat("MASS::fitdistr has no built-in named route for Rician.\n")

rician_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_rician, "rician", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2802,
    control = tutorial_control
  )
)
cat("\nObjective Bayes:\n")
show_bayes(rician_bayes)

# =============================================================================
# 20. WEIGHTED LINDLEY
# =============================================================================
cat("\n\n================ WEIGHTED LINDLEY ============\n")
set.seed(2901)
lambda_wl <- 2.5
phi_wl <- 0.8
component_wl <- runif(100) < lambda_wl / (lambda_wl + phi_wl)
x_wl <- rgamma(
  100, shape = phi_wl + as.numeric(!component_wl), rate = lambda_wl
)

dweighted_lindley <- function(x, lambda, phi, log = FALSE) {
  log_density <- (phi + 1) * log(lambda) - log(lambda + phi) -
    lgamma(phi) + (phi - 1) * log(x) + log1p(x) - lambda * x
  if (log) log_density else exp(log_density)
}

cat("\nMaximum likelihood:\n")
print(suppressWarnings(MASS::fitdistr(
  x_wl, dweighted_lindley, start = list(lambda = 2, phi = 1)
)))

weighted_lindley_bayes <- list(
  Jeffreys = fitdistrBayes(
    x_wl, "weighted lindley", "jeffreys", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2902,
    control = tutorial_control
  ),
  One_group_reference = fitdistrBayes(
    x_wl, "weighted lindley", "reference", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2903,
    control = tutorial_control
  ),
  Jeffreys_first_rule = fitdistrBayes(
    x_wl, "weighted lindley", "first-rule", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2904,
    control = tutorial_control
  ),
  Independence_Jeffreys = fitdistrBayes(
    x_wl, "weighted lindley", "independence-jeffreys",
    iter = tutorial_iter, warmup = tutorial_warmup,
    chains = tutorial_chains, seed = 2905, control = tutorial_control
  ),
  Reference_lambda = fitdistrBayes(
    x_wl, "weighted lindley", "reference-lambda", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2906,
    control = tutorial_control
  ),
  Reference_phi = fitdistrBayes(
    x_wl, "weighted lindley", "reference-phi", iter = tutorial_iter,
    warmup = tutorial_warmup, chains = tutorial_chains, seed = 2907,
    control = tutorial_control
  )
)
cat("\nObjective Bayes, grouped by prior:\n")
show_bayes(weighted_lindley_bayes)
cat("\nThe weighted Lindley MDI posterior is improper and is not fitted.\n")

# =============================================================================
# FINAL TABLE: all fits and diagnostics
# =============================================================================
all_fits <- c(
  Beta = beta_bayes,
  Cauchy = cauchy_bayes,
  Chi_squared = chisq_bayes,
  Exponential = exponential_bayes,
  Gamma = gamma_bayes,
  Geometric = geometric_bayes,
  Lognormal = lognormal_bayes,
  Logistic = logistic_bayes,
  Negative_binomial = negbin_bayes,
  Normal = normal_bayes,
  Poisson = poisson_bayes,
  Student_t = t_bayes,
  Weibull = weibull_bayes,
  Gumbel = gumbel_bayes,
  Frechet = frechet_bayes,
  Lomax = lomax_bayes,
  Nakagami_m = nakagami_bayes,
  Exponential_logarithmic = el_bayes,
  Rician = rician_bayes,
  Weighted_lindley = weighted_lindley_bayes
)

diagnostics <- do.call(rbind, lapply(
  names(all_fits),
  function(name) {
    fit <- all_fits[[name]]
    data.frame(
      fit = name,
      model = fit$model$name,
      prior = fit$prior$input,
      algorithm = fit$engine$algorithm,
      converged = fit$diagnostics$converged,
      max_rhat = round(fit$diagnostics$max_rhat, 3),
      min_ess_bulk = round(fit$diagnostics$min_ess_bulk, 1),
      min_ess_tail = round(fit$diagnostics$min_ess_tail, 1),
      row.names = NULL
    )
  }
))

cat("\n\n================== FINAL DIAGNOSTICS =================\n")
print(diagnostics, row.names = FALSE)

if (any(!diagnostics$converged)) {
  cat(
    "\nSome quick fits did not satisfy the R-hat/ESS thresholds. ",
    "This does not indicate a function error.\n",
    "Set quick_mode to FALSE and rerun before interpreting these fits.\n",
    sep = ""
  )
}

# =============================================================================
# HOW TO EXAMINE AN INDIVIDUAL FIT
# =============================================================================
cat("\n\n==================== A FIT IN DETAIL ==================\n")

# Complete summary:
print(summary(gamma_bayes$Jeffreys))

# Credible intervals:
print(confint(gamma_bayes$Jeffreys))

# Starting values and automatic classical method used:
print(gamma_bayes$Jeffreys$initialization)
print(weibull_bayes$Jeffreys$initialization)
print(weighted_lindley_bayes$Reference_lambda$initialization)

# First posterior draws:
print(head(as.data.frame(gamma_bayes$Jeffreys)))

# Five predictive samples, each containing 10 observations:
print(predict(gamma_bayes$Jeffreys, draws = 5, size = 10, seed = 3001))

# Pointwise log-likelihood for three draws. Only the first six observations
# are displayed to keep the console output compact.
loglik_gamma <- log_lik(gamma_bayes$Jeffreys, draws = 3)
print(loglik_gamma[, 1:6, drop = FALSE])

# Plots are opened only in an interactive session.
if (interactive()) {
  plot(gamma_bayes$Jeffreys, type = "trace")
  plot(gamma_bayes$Jeffreys, type = "density")
}

cat("\nTutorial completed.\n")
