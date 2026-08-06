library("fitdistrBayes")

assert_mc <- function(estimate, target, mcse, label, multiplier = 8) {
  tolerance <- multiplier * mcse + 2e-4 * max(1, abs(target))
  if (!is.finite(estimate + target + mcse) ||
      abs(estimate - target) > tolerance) {
    stop(sprintf(
      "%s failed: estimate=%g, quadrature=%g, MCSE=%g, tolerance=%g.",
      label, estimate, target, mcse, tolerance
    ), call. = FALSE)
  }
}

log_sum_exp <- function(z) {
  m <- max(z)
  m + log(sum(exp(z - m)))
}

quadrature_mean_exp <- function(log_density, search = c(-12, 12)) {
  objective <- function(eta) log_density(eta)
  mode <- optimize(objective, interval = search, maximum = TRUE)
  height <- mode$objective
  scaled <- function(eta) {
    vapply(eta, function(z) {
      value <- log_density(z)
      if (is.finite(value)) exp(value - height) else 0
    }, numeric(1))
  }
  denominator <- integrate(scaled, -Inf, Inf, rel.tol = 2e-9,
                           subdivisions = 2000L)$value
  scaled_moment <- function(eta) {
    vapply(eta, function(z) {
      value <- log_density(z)
      log_value <- z + value - height
      if (is.finite(log_value)) exp(log_value) else 0
    }, numeric(1))
  }
  numerator <- integrate(scaled_moment,
                         -Inf, Inf, rel.tol = 2e-9,
                         subdivisions = 2000L)$value
  numerator / denominator
}

settings <- list(
  iter = 8000L, warmup = 2000L, chains = 4L,
  control = list(rhat_threshold = 1.02, ess_threshold = 800,
                 warn_convergence = FALSE)
)

# Chi-squared: direct quadrature of the posterior for eta = log(df).
x_chisq <- c(0.8, 1.4, 2.1, 3.0, 4.2, 5.5, 6.1, 7.8)
log_chisq <- function(eta) {
  if (!is.finite(eta) || abs(eta) > 600) return(-Inf)
  df <- exp(eta)
  sum(dchisq(x_chisq, df = df, log = TRUE)) +
    0.5 * log(suppressWarnings(trigamma(df / 2))) + eta
}
target_chisq <- quadrature_mean_exp(log_chisq)
fit_chisq <- fitdistrBayes(
  x_chisq, "chi-squared", "Jeffreys", seed = 801,
  iter = settings$iter, warmup = settings$warmup,
  chains = settings$chains, control = settings$control
)
row <- fit_chisq$summary[fit_chisq$summary$parameter == "df", ]
assert_mc(row$mean, target_chisq, row$mcse_mean,
          "Chi-squared posterior mean")

# Gamma/reference-shape: quadrature of the analytically marginalized shape.
x_gamma <- c(0.35, 0.62, 0.91, 1.20, 1.57, 2.05, 2.80, 3.60)
n_gamma <- length(x_gamma)
S_gamma <- sum(x_gamma)
Plog_gamma <- sum(log(x_gamma))
gamma_joint_term <- function(a) {
  if (a > 1e5) {
    return(1 / (2 * a) + 1 / (6 * a^2) - 1 / (30 * a^4))
  }
  a * suppressWarnings(trigamma(a)) - 1
}
log_gamma_shape <- function(eta) {
  if (!is.finite(eta) || abs(eta) > 600) return(-Inf)
  a <- exp(eta)
  0.5 * (log(gamma_joint_term(a)) - log(a)) +
    lgamma(n_gamma * a) - n_gamma * lgamma(a) +
    (a - 1) * Plog_gamma - n_gamma * a * log(S_gamma) + eta
}
target_gamma <- quadrature_mean_exp(log_gamma_shape)
fit_gamma <- fitdistrBayes(
  x_gamma, "gamma", "reference-shape", seed = 802,
  iter = settings$iter, warmup = settings$warmup,
  chains = settings$chains, control = settings$control
)
row <- fit_gamma$summary[fit_gamma$summary$parameter == "shape", ]
assert_mc(row$mean, target_gamma, row$mcse_mean,
          "Gamma marginal shape mean")

# Frechet: direct quadrature of the marginal induced by
# L(shape, scale)/(shape*scale).  The natural-scale marginal contains
# shape^(n-1); this check prevents loss of that factor during transformation.
x_frechet <- c(0.72, 0.94, 1.18, 1.46, 1.91, 2.55, 3.42, 5.10)
n_frechet <- length(x_frechet)
lx_frechet <- log(x_frechet)
log_frechet_shape <- function(eta) {
  if (!is.finite(eta) || abs(eta) > 600) return(-Inf)
  shape <- exp(eta)
  log_A <- log_sum_exp(-shape * lx_frechet)
  n_frechet * eta - shape * sum(lx_frechet) - n_frechet * log_A
}
target_frechet <- quadrature_mean_exp(log_frechet_shape)
fit_frechet <- fitdistrBayes(
  x_frechet, "frechet", "jeffreys", seed = 804,
  iter = settings$iter, warmup = settings$warmup,
  chains = settings$chains, control = settings$control
)
row <- fit_frechet$summary[fit_frechet$summary$parameter == "shape", ]
assert_mc(row$mean, target_frechet, row$mcse_mean,
          "Frechet marginal shape mean")

# Weibull/reference: quadrature of the marginalized shape posterior.
x_weibull <- c(0.42, 0.68, 0.95, 1.30, 1.75, 2.20, 2.95, 3.80)
n_weibull <- length(x_weibull)
lx <- log(x_weibull)
log_weibull_shape <- function(eta) {
  if (!is.finite(eta) || abs(eta) > 700) return(-Inf)
  k <- exp(eta)
  log_A <- log_sum_exp(k * lx)
  (n_weibull - 2) * log(k) + (k - 1) * sum(lx) -
    n_weibull * log_A + eta
}
target_weibull <- quadrature_mean_exp(log_weibull_shape)
fit_weibull <- fitdistrBayes(
  x_weibull, "weibull", "reference", seed = 803,
  iter = settings$iter, warmup = settings$warmup,
  chains = settings$chains, control = settings$control
)
row <- fit_weibull$summary[fit_weibull$summary$parameter == "shape", ]
assert_mc(row$mean, target_weibull, row$mcse_mean,
          "Weibull marginal shape mean")

cat(sprintf(
  paste0("Quadrature validations passed. Posterior means: ",
         "chi-squared df %.6f; Gamma shape %.6f; Frechet shape %.6f; ",
         "Weibull shape %.6f.\n"),
  target_chisq, target_gamma, target_frechet, target_weibull
))
