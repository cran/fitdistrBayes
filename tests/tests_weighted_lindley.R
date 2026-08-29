library("fitdistrBayes")

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

expect_error <- function(code, pattern) {
  message <- tryCatch({
    force(code)
    NA_character_
  }, error = function(e) conditionMessage(e))
  if (is.na(message) || !grepl(pattern, message, ignore.case = TRUE)) {
    stop(sprintf("Expected an error matching '%s'; got '%s'.",
                 pattern, message), call. = FALSE)
  }
}

dwl <- function(x, lambda, phi, log = FALSE) {
  log_density <- (phi + 1) * log(lambda) - log(lambda + phi) -
    lgamma(phi) + (phi - 1) * log(x) + log1p(x) - lambda * x
  if (log) log_density else exp(log_density)
}

rwl <- function(n, lambda, phi) {
  first <- runif(n) < lambda / (lambda + phi)
  rgamma(n, shape = phi + as.numeric(!first), rate = lambda)
}

# The density and two-Gamma mixture have the claimed normalization and moments.
for (theta in list(c(lambda = 0.4, phi = 0.35),
                   c(lambda = 1, phi = 2),
                   c(lambda = 8, phi = 12))) {
  integral <- integrate(
    dwl, 0, Inf, lambda = theta["lambda"], phi = theta["phi"],
    rel.tol = 1e-10
  )$value
  assert(abs(integral - 1) < 2e-8,
         "Weighted Lindley density did not integrate to one.")
}

set.seed(280826)
lambda <- 1.7
phi <- 0.65
large_sample <- rwl(200000, lambda, phi)
truth_mean <- phi * (lambda + phi + 1) / (lambda * (lambda + phi))
truth_variance <- (phi + 1) / lambda^2 - 1 / (lambda + phi)^2
assert(abs(mean(large_sample) - truth_mean) < 0.015,
       "Weighted Lindley RNG mean is inconsistent with the mixture moment.")
assert(abs(var(large_sample) - truth_variance) < 0.025,
       "Weighted Lindley RNG variance is inconsistent with the mixture moment.")

# Stable prior expressions agree with the natural-scale formulas.
log_prior_stable <- getFromNamespace(
  ".fdb_weighted_lindley_log_prior", "fitdistrBayes"
)
lambda_from_mean <- getFromNamespace(
  ".fdb_weighted_lindley_lambda_from_mean", "fitdistrBayes"
)
grid <- expand.grid(lambda = c(0.08, 0.7, 4, 30),
                    phi = c(0.09, 0.8, 3, 25))
for (i in seq_len(nrow(grid))) {
  lambda <- grid$lambda[i]
  phi <- grid$phi[i]
  ratio_log <- log(lambda / phi)
  phi_log <- log(phi)
  tri <- trigamma(phi)
  A <- (lambda + phi)^2 + 2 * lambda + phi
  c_phi <- phi * tri - 1
  s_phi <- sqrt(phi + 1) + sqrt(phi)
  w_zero <- sqrt(phi + 1) / s_phi
  w_infinity <- sqrt(phi) / s_phi
  q_zero <- trigamma(phi + 1) - 1 / (phi + 1)
  q_infinity <- tri - 1 / phi
  natural <- c(
    jeffreys = 0.5 * log(A * c_phi - 1) - log(lambda) -
      log(lambda + phi),
    reference = 0.5 * log(A * c_phi - 1) - log(lambda) -
      log(lambda + phi),
    `first-rule` = -log(lambda) - log(phi),
    `independence-jeffreys` = 0.5 * log(
      ((phi + 1) * (lambda + phi)^2 - lambda^2) *
        (tri * (lambda + phi)^2 - 1)
    ) - log(lambda) - 2 * log(lambda + phi),
    `reference-lambda` = -log(lambda) +
      0.5 * log(tri - 1 / (lambda + phi)^2),
    `reference-phi` = 0.5 * log(phi * A) - log(lambda) -
      log(lambda + phi) - log(s_phi) +
      0.5 * w_zero * log(q_zero) +
      0.5 * w_infinity * log(q_infinity)
  )
  for (prior in names(natural)) {
    stable <- log_prior_stable(lambda, phi, prior)
    assert(abs(stable - natural[[prior]]) < 2e-10,
           sprintf("Stable %s prior is inconsistent.", prior))
  }
  mean_parameter <- phi * (lambda + phi + 1) /
    (lambda * (lambda + phi))
  assert(abs(lambda_from_mean(mean_parameter, phi) - lambda) <
           2e-12 * (1 + lambda),
         "The Fisher-orthogonal mean transformation is not invertible.")
}

# All six proper objective-prior routes produce finite fitted objects,
# log-likelihoods, predictive draws, and certified posterior moments.
set.seed(280827)
x <- rwl(45, lambda = 2.5, phi = 0.7)
quick <- list(rhat_threshold = 5, ess_threshold = 1,
              warn_convergence = FALSE)
priors <- c("jeffreys", "reference", "first-rule",
            "independence-jeffreys", "reference-lambda", "reference-phi")
for (position in seq_along(priors)) {
  fit <- fitdistrBayes(
    x, "weighted Lindley", priors[position],
    iter = 500, warmup = 250, chains = 2, seed = 2800 + position,
    control = quick
  )
  assert(inherits(fit, "fitdistrBayes"),
         "Weighted Lindley route did not return a fitted object.")
  assert(identical(fit$model$parameters, c("lambda", "phi")),
         "Weighted Lindley output parameter names changed.")
  assert(all(is.finite(as.matrix(fit$draws[, c("lambda", "phi")]))),
         "Weighted Lindley posterior draws are not finite.")
  assert(all(fit$draws$lambda > 0 & fit$draws$phi > 0),
         "Weighted Lindley posterior draws left the positive support.")
  assert(all(is.finite(log_lik(fit, draws = 20, seed = 2))),
         "Weighted Lindley pointwise log likelihood is invalid.")
  predictive <- predict(fit, draws = 20, size = 3, seed = 3)
  assert(all(is.finite(predictive) & predictive > 0),
         "Weighted Lindley predictive draws are invalid.")
  assert(all(fit$moment_status$mean_exists) &&
           all(fit$moment_status$variance_exists),
         "Weighted Lindley moment certification is incomplete.")
  assert(grepl("weighted Lindley", fit$initialization$method,
               ignore.case = TRUE),
         "Weighted Lindley automatic initialization was not recorded.")
  assert(identical(
    fit$engine$algorithm,
    "adaptive Metropolis in Fisher-orthogonal mean coordinates"
  ), "Weighted Lindley did not use the orthogonal MCMC parameterization.")
}

# Aliases, partial user starts, and all pre-sampling guards remain executable.
fit_alias <- fitdistrBayes(
  x, "wl", "reference", start = c(lambda = 2),
  iter = 240, warmup = 120, chains = 2, seed = 2810,
  control = quick
)
assert(identical(fit_alias$model$name, "weighted lindley"),
       "The WL alias was not normalized.")
assert(identical(fit_alias$initialization$source,
                 "automatic with user overrides"),
       "A partial weighted Lindley start was not recorded.")

expect_error(
  fitdistrBayes(x, "weighted lindley", "mdi",
                iter = 100, warmup = 50, chains = 1),
  "improper"
)
expect_error(
  fitdistrBayes(c(1, 1), "weighted lindley", "reference",
                iter = 100, warmup = 50, chains = 1),
  "nonconstant"
)
expect_error(
  fitdistrBayes(c(0, 1), "weighted lindley", "reference",
                iter = 100, warmup = 50, chains = 1),
  "positive"
)
expect_error(
  fitdistrBayes(1, "weighted lindley", "reference",
                iter = 100, warmup = 50, chains = 1),
  "n >= 2"
)

cat("Weighted Lindley analytic, posterior, predictive, and guard tests passed.\n")
