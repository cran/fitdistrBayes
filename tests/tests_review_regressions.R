library("fitdistrBayes")

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  assert(!is.null(message) && grepl(pattern, message, ignore.case = TRUE),
         sprintf("Expected error matching '%s'; received: %s",
                 pattern, if (is.null(message)) "<no error>" else message))
}

ctrl <- list(rhat_threshold = 10, ess_threshold = 1,
             warn_convergence = FALSE)

# Constant positive Nakagami-m samples have improper marginal shape tails.
for (prior in c("jeffreys", "reference")) {
  expect_error(
    fitdistrBayes(rep(1, 4), "nakagami-m", prior,
                  iter = 100, warmup = 50, chains = 1, seed = 1,
                  control = ctrl),
    "nonconstant"
  )
}
x_near_constant <- c(1, 1 + .Machine$double.eps, 1, 1)
fit_near_constant <- fitdistrBayes(
  x_near_constant, "nakagami-m", "reference",
  iter = 200, warmup = 100, chains = 2, seed = 2,
  control = ctrl
)
assert(inherits(fit_near_constant, "fitdistrBayes"),
       "A representably nonconstant Nakagami-m sample was rejected.")

# A custom function with ... but no explicit log formal is ordinary-density
# mode by default. It must never be silently interpreted as log-density mode.
d_ignore_log <- function(x, location, ...) stats::dnorm(x, location, 1)
p_ignore_log <- function(location, ...) stats::dnorm(location, 0, 10)
fit_ordinary <- suppressWarnings(fitdistrBayes(
  c(-1, 0, 1), d_ignore_log, p_ignore_log,
  start = c(location = 0), iter = 120, warmup = 60,
  chains = 2, seed = 3, control = ctrl
))
assert(!fit_ordinary$model$custom_contract$density_is_log &&
       !fit_ordinary$model$custom_contract$prior_is_log,
       "Custom ordinary-density mode was not recorded correctly.")
expect_error(
  suppressWarnings(fitdistrBayes(
    c(-1, 0, 1), d_ignore_log, p_ignore_log,
    start = c(location = 0), iter = 100, warmup = 50,
    chains = 1, seed = 4,
    control = c(ctrl, list(density_is_log = TRUE))
  )),
  "log-density contract"
)
expect_error(
  suppressWarnings(fitdistrBayes(
    c(-1, 0, 1), d_ignore_log, p_ignore_log,
    start = c(location = 0), iter = 100, warmup = 50,
    chains = 1, seed = 5,
    control = c(ctrl, list(prior_is_log = TRUE))
  )),
  "log-density contract"
)

d_explicit_log <- function(x, location, log = FALSE) {
  stats::dnorm(x, location, 1, log = log)
}
p_explicit_log <- function(location, log = FALSE) {
  stats::dnorm(location, 0, 10, log = log)
}
fit_logged <- suppressWarnings(fitdistrBayes(
  c(-1, 0, 1), d_explicit_log, p_explicit_log,
  start = c(location = 0), iter = 120, warmup = 60,
  chains = 2, seed = 6, control = ctrl
))
assert(fit_logged$model$custom_contract$density_is_log &&
       fit_logged$model$custom_contract$prior_is_log,
       "Explicit custom log-density mode was not recorded correctly.")

# -Inf is a valid log-density outside support; +Inf, NA, and NaN are not.
valid_log <- fitdistrBayes:::.fdb_valid_log_values
assert(identical(valid_log(-Inf, 1L, "prior"), -Inf),
       "-Inf was rejected as a valid zero-density log value.")
for (bad in list(Inf, NA_real_, NaN)) {
  expect_error(valid_log(bad, 1L, "prior"), "Inf, NA, and NaN")
}

# Representable half-integers remain nonintegers regardless of magnitude.
assert(!fitdistrBayes:::.fdb_discrete(2^45 + 0.5),
       "A large representable half-integer was classified as discrete.")
assert(fitdistrBayes:::.fdb_discrete(2^45),
       "A large representable integer was rejected.")

# Constant chains at different values are a convergence failure, not success.
different_constants <- matrix(rep(0:3, each = 40), nrow = 40, ncol = 4)
same_constants <- matrix(2, nrow = 40, ncol = 4)
assert(is.infinite(fitdistrBayes:::.fdb_rhat(different_constants)),
       "Different constant chains did not produce infinite R-hat.")
assert(identical(fitdistrBayes:::.fdb_ess_matrix(different_constants), 0),
       "Different constant chains did not produce zero ESS.")
assert(identical(fitdistrBayes:::.fdb_rhat(same_constants), 1),
       "Globally identical constant chains did not use the documented R-hat convention.")
assert(fitdistrBayes:::.fdb_ess_matrix(same_constants) == 160,
       "Globally identical constant chains did not use the documented ESS convention.")

# Metropolis acceptance accounting separates warmup from post-warmup draws.
set.seed(7)
amwg <- fitdistrBayes:::.fdb_amwg(
  function(z) -sum(z^2) / 2, c(theta = 0),
  iter = 120, warmup = 40, thin = 1, chains = 2,
  control = fitdistrBayes:::.fdb_merge_control(ctrl)
)
lhs <- amwg$acceptance_all * 120
rhs <- amwg$acceptance_warmup * 40 + amwg$acceptance * 80
assert(max(abs(lhs - rhs)) < 1e-10,
       "Warmup/post-warmup acceptance counters are inconsistent.")

# Count controls reject fractions, and the removed dead option is unknown.
expect_error(
  fitdistrBayes:::.fdb_merge_control(list(entropy_exact_limit = 2.5)),
  "positive integer"
)
expect_error(
  fitdistrBayes:::.fdb_merge_control(list(max_entropy_terms = 100L)),
  "Unknown control"
)

# Custom predictive values must be finite and may be support-validated.
fit_bad_rng <- suppressWarnings(fitdistrBayes(
  c(-1, 0, 1), d_explicit_log, p_explicit_log,
  start = c(location = 0), iter = 100, warmup = 50,
  chains = 1, seed = 8,
  control = c(ctrl, list(rng = function(n, location) rep(NA_real_, n)))
))
expect_error(predict(fit_bad_rng, draws = 2, seed = 8), "finite numeric")

fit_bad_support <- suppressWarnings(fitdistrBayes(
  c(-1, 0, 1), d_explicit_log, p_explicit_log,
  start = c(location = 0), iter = 100, warmup = 50,
  chains = 1, seed = 9,
  control = c(ctrl, list(
    rng = function(n, location) rep(-1, n),
    rng_validator = function(x) x >= 0
  ))
))
expect_error(predict(fit_bad_support, draws = 2, seed = 9),
             "rng_validator")

# Disabling stored callables is explicit and survives serialization.
fit_no_callables <- fitdistrBayes(
  rexp(10), "exponential", "reference",
  iter = 100, warmup = 50, chains = 1, seed = 10,
  control = list(store_callables = FALSE)
)
assert(!fit_no_callables$capabilities$prediction &&
       !fit_no_callables$capabilities$pointwise_log_likelihood,
       "Disabled callable capabilities were reported as available.")
path <- tempfile(fileext = ".rds")
saveRDS(fit_no_callables, path)
fit_reloaded <- readRDS(path)
unlink(path)
expect_error(predict(fit_reloaded, draws = 1), "unavailable")
expect_error(log_lik(fit_reloaded, draws = 1), "disabled")

cat("All independent-review regression tests passed.\n")
