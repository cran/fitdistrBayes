library("fitdistrBayes")

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  assert(!is.null(message) && grepl(pattern, message, ignore.case = TRUE),
         sprintf("Expected error matching '%s'.", pattern))
}

ctrl <- list(rhat_threshold = 2, ess_threshold = 1,
             warn_convergence = FALSE)
fit_small <- function(x, distr, prior, fixed = NULL, ...) {
  fitdistrBayes(x, distr, prior, fixed = fixed, ...,
                iter = 300, warmup = 150, chains = 2, seed = 88,
                control = ctrl)
}

# Exact boundary samples that remain proper.
fit_zero_pois <- fit_small(rep(0, 10), "Poisson", "Jeffreys")
assert(all(fit_zero_pois$draws$lambda > 0), "All-zero Poisson failed.")
fit_zero_geom <- fit_small(rep(0, 10), "geometric", "reference")
assert(all(fit_zero_geom$draws$prob > 0 &
           fit_zero_geom$draws$prob < 1), "All-zero geometric failed.")
fit_zero_nb <- fit_small(rep(0, 10), "nbinom", "reference",
                         fixed = list(size = 2.5))
assert(all(fit_zero_nb$draws$mu > 0), "All-zero negative binomial failed.")
fit_exp <- fit_small(c(0, 0, 0.5), "exponential", "MDI")
assert(all(fit_exp$draws$rate > 0), "Exponential sample containing zeros failed.")

# Pointwise log-likelihood always keeps draws x observations dimensions,
# including the vapply simplification boundary at a one-observation sample.
fit_one_exp <- fit_small(0.5, "exponential", "reference")
ll_one <- log_lik(fit_one_exp, draws = 7, seed = 19)
assert(identical(dim(ll_one), c(7L, 1L)) && all(is.finite(ll_one)),
       "One-observation log_lik() did not return a 7 x 1 matrix.")

# Scale robustness.
fit_gamma_tiny <- fit_small(c(1e-100, 2e-100, 4e-100, 8e-100),
                            "gamma", "reference-rate")
assert(all(is.finite(fit_gamma_tiny$draws$rate)),
       "Tiny-scale Gamma data caused non-finite draws.")
fit_weibull_tiny <- fit_small(c(1e-80, 2e-80, 5e-80, 9e-80),
                              "weibull", "reference")
assert(all(is.finite(fit_weibull_tiny$draws$scale)),
       "Tiny-scale Weibull data caused non-finite draws.")
fit_lnorm_tiny <- fit_small(c(1e-100, 2e-100, 4e-100, 9e-100),
                            "log-normal", "Jeffreys")
assert(all(is.finite(fit_lnorm_tiny$draws$meanlog)),
       "Tiny-scale Lognormal data caused non-finite draws.")

# Proper heavy-tailed posteriors at n = 2 can exceed floating-point range.
# The package must stop transparently instead of silently placing atoms at
# machine limits.
x_gamma_tail <- c(1.3679237230556981e-4, 2.7959031656389898e-9)
fit_gamma_tail <- fit_small(x_gamma_tail, "gamma", "reference-rate")
assert(all(is.finite(fit_gamma_tail$draws$rate) &
           fit_gamma_tail$draws$rate > 0),
       "Representable heavy-tailed Gamma draws were not retained.")
x_weibull_tail <- c(6.8554240217917823e-4, 1.8482839074829358e-1)
expect_error(
  fitdistrBayes(
    x_weibull_tail, "weibull", "reference",
    iter = 240, warmup = 120, chains = 2, seed = 20261021,
    control = ctrl
  ),
  "no clipping"
)
fit_weibull_regular <- fit_small(
  c(0.5, 1, 2, 3), "weibull", "reference"
)
weibull_impossible <- fit_weibull_regular$.loglik(
  c(shape = 1000, scale = .Machine$double.xmin)
)
assert(any(weibull_impossible == -Inf) && !any(is.nan(weibull_impossible)),
       "Impossible-tail Weibull log likelihood was clipped instead of returning -Inf.")

# Translation invariance of the nondegeneracy check. These observations are
# close relative to their large location, but remain distinct in double
# precision and therefore form a valid nonconstant sample.
x_large_offset <- 1e12 + c(-0.002, -0.001, 0, 0.001, 0.002)
fit_normal_offset <- fit_small(x_large_offset, "normal", "reference")
assert(all(is.finite(fit_normal_offset$draws$sd)),
       "Large-offset Normal data were treated as constant.")
fit_logistic_offset <- fit_small(x_large_offset, "logistic", "reference")
assert(all(is.finite(fit_logistic_offset$draws$scale)),
       "Large-offset Logistic data were treated as constant.")
fit_t_offset <- fit_small(
  x_large_offset, "t", "reference", fixed = list(df = 6)
)
assert(all(is.finite(fit_t_offset$draws$scale)),
       "Large-offset Student-t data were treated as constant.")

# Exact multiplicity boundaries.
expect_error(
  fit_small(c(0, 0, 1, 2), "cauchy", "Jeffreys"),
  "2m < n"
)
fit_cauchy_ref <- fit_small(c(0, 0, 1, 2), "cauchy", "reference")
assert(inherits(fit_cauchy_ref, "fitdistrBayes"),
       "Cauchy reference multiplicity boundary was rejected incorrectly.")

# Missing-value policy.
fit_omit <- fitdistrBayes(
  c(-1, NA, 0, 1, 2), "normal", "reference",
  na.action = "omit", iter = 300, warmup = 150, chains = 2, seed = 9
)
assert(identical(fit_omit$omitted, 2L) && fit_omit$model$n == 4,
       "na.action='omit' failed.")

# User-defined probability parameter with two finite bounds.
dbern <- function(x, prob, log = FALSE) {
  value <- ifelse(x == 1, prob, ifelse(x == 0, 1 - prob, 0))
  if (log) log(value) else value
}
prior_jeff <- function(prob, log = FALSE) {
  value <- if (prob > 0 && prob < 1) 1 / sqrt(prob * (1 - prob)) else 0
  if (log) log(value) else value
}
fit_custom <- suppressWarnings(fitdistrBayes(
  c(1, 0, 1, 1, 0, 1, 0, 1),
  dbern, prior_jeff, start = c(prob = 0.5),
  iter = 500, warmup = 250, chains = 2, seed = 10,
  control = list(
    lower = c(prob = 0), upper = c(prob = 1),
    rhat_threshold = 2, ess_threshold = 1, warn_convergence = FALSE,
    rng = function(n, prob) rbinom(n, 1, prob)
  )
))
assert(all(fit_custom$draws$prob > 0 & fit_custom$draws$prob < 1),
       "Two-sided custom transformation failed.")
assert(all(predict(fit_custom, draws = 20, seed = 1) %in% 0:1),
       "Custom posterior prediction failed.")

# Defensive controls and starting values.
expect_error(
  fitdistrBayes(1:5, "Poisson", "Jeffreys",
                iter = 100, warmup = 50,
                control = list(rhat_threshold = 0.9)),
  "at least 1"
)
expect_error(
  fitdistrBayes(rbeta(10, 2, 3), "beta", "Jeffreys",
                start = c(shape1 = -1, shape2 = 2),
                iter = 100, warmup = 50),
  "Positive starting"
)
expect_error(
  fitdistrBayes(as.complex(1:5), "Poisson", "Jeffreys",
                iter = 100, warmup = 50),
  "numeric vector"
)
expect_error(
  fitdistrBayes(1:5, "Poisson", "Jeffreys",
                fixed = structure(list(1, 2), names = c("size", "size")),
                iter = 100, warmup = 50),
  "unique"
)
expect_error(
  fitdistrBayes(1:5, "Poisson", "Jeffreys",
                seed = 1.5, iter = 100, warmup = 50),
  "integer"
)
expect_error(
  fitdistrBayes(
    c(0, 1), dbern, prior_jeff, start = c(prob = 0.5),
    iter = 100, warmup = 50,
    control = list(lower = c(prob = NA_real_), upper = c(prob = 1))
  ),
  "real numeric"
)

cat("All numerical edge-case tests passed.\n")
