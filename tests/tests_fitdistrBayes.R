library("fitdistrBayes")

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

expect_error <- function(expr, pattern = NULL) {
  got <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  assert(!is.null(got), "Expected an error, but the expression succeeded.")
  if (!is.null(pattern)) {
    assert(grepl(pattern, got, ignore.case = TRUE),
           sprintf("Error '%s' did not match '%s'.", got, pattern))
  }
  invisible(got)
}

check_fit <- function(fit, parameters, n_save, chains) {
  assert(inherits(fit, "fitdistrBayes"), "Wrong fitted-object class.")
  assert(identical(fit$model$parameters, parameters),
         "Unexpected parameter names.")
  assert(length(fit$chains) == chains, "Wrong number of chains.")
  assert(all(vapply(fit$chains, nrow, integer(1)) == n_save),
         "Wrong saved chain length.")
  assert(all(is.finite(as.matrix(fit$draws[, parameters, drop = FALSE]))),
         "Non-finite posterior draws.")
  certified_means <- fit$summary$mean_exists %in% TRUE
  assert(all(is.finite(fit$summary$mean[certified_means])),
         "A certified posterior mean was not finite.")
  assert(all(fit$summary$q2.5 <= fit$summary$median &
             fit$summary$median <= fit$summary$q97.5),
         "Invalid posterior quantile order.")
  invisible(fit)
}

fast <- list(
  iter = 700L,
  warmup = 300L,
  chains = 4L,
  seed = 20260719,
  control = list(
    rhat_threshold = 1.20,
    ess_threshold = 20,
    warn_convergence = FALSE
  )
)
n_save <- fast$iter - fast$warmup

run_fit <- function(x, distr, prior, fixed = NULL, start = NULL,
                    control = fast$control) {
  fitdistrBayes(
    x, distr, prior, fixed = fixed, start = start,
    iter = fast$iter, warmup = fast$warmup, chains = fast$chains,
    seed = fast$seed, control = control
  )
}

set.seed(1001)
fits <- list()
fits$beta <- run_fit(rbeta(35, 2, 5), "BeTa", "Jeffreys")
check_fit(fits$beta, c("shape1", "shape2"), n_save, fast$chains)

fits$cauchy <- run_fit(rcauchy(35, 1, 2), "cauchy", "reference")
check_fit(fits$cauchy, c("location", "scale"), n_save, fast$chains)

fits$chisq <- run_fit(rchisq(30, 4), "chi-square", "reference")
check_fit(fits$chisq, "df", n_save, fast$chains)

x_exp <- rexp(40, 2.5)
fits$exp <- run_fit(x_exp, "EXPONENTIAL", "Jeffreys")
check_fit(fits$exp, "rate", n_save, fast$chains)
expected_exp <- length(x_exp) / sum(x_exp)
mcse_exp <- fits$exp$summary$mcse_mean
assert(abs(fits$exp$summary$mean[1] - expected_exp) < 6 * mcse_exp,
       "Exponential exact posterior mean check failed.")

fits$gamma <- run_fit(rgamma(40, 3, rate = 2), "gamma",
                      "reference-shape")
check_fit(fits$gamma, c("shape", "rate"), n_save, fast$chains)

fits$geometric <- run_fit(rgeom(35, 0.35), "geometric", "MDI")
check_fit(fits$geometric, "prob", n_save, fast$chains)

fits$lognormal <- run_fit(rlnorm(35, 0.5, 0.7), "log-normal",
                          "reference")
check_fit(fits$lognormal, c("meanlog", "sdlog"), n_save, fast$chains)

fits$logistic <- run_fit(rlogis(35, -1, 1.5), "logistic", "MDI")
check_fit(fits$logistic, c("location", "scale"), n_save, fast$chains)

fits$nbinom <- run_fit(rnbinom(35, size = 4, mu = 6),
                       "negative binomial", "Jeffreys",
                       fixed = list(size = 4))
check_fit(fits$nbinom, "mu", n_save, fast$chains)

fits$normal <- run_fit(rnorm(35, 2, 1.5), "normal", "MDI")
check_fit(fits$normal, c("mean", "sd"), n_save, fast$chains)

x_pois <- rpois(40, 3)
fits$poisson <- run_fit(x_pois, "Poisson", "MDI")
check_fit(fits$poisson, "lambda", n_save, fast$chains)

fits$t_fixed <- run_fit(1 + 2 * rt(40, 6), "t", "reference",
                        fixed = list(df = 6))
check_fit(fits$t_fixed, c("location", "scale"), n_save, fast$chains)

fits$t_unknown <- run_fit(rt(45, 8), "student-t",
                          "independence-jeffreys")
check_fit(fits$t_unknown, c("location", "scale", "df"),
          n_save, fast$chains)

fits$weibull <- run_fit(rweibull(40, 1.7, 2), "weibull", "reference")
check_fit(fits$weibull, c("shape", "scale"), n_save, fast$chains)
assert(is.na(fits$weibull$summary$mean[fits$weibull$summary$parameter ==
                                      "scale"]),
       "The nonexistent Weibull scale mean was not suppressed.")
assert(is.na(fits$t_unknown$summary$mean[fits$t_unknown$summary$parameter ==
                                        "df"]),
       "The nonexistent unknown-df Student-t mean was not suppressed.")
assert(isTRUE(all.equal(
  unname(coef(fits$gamma)),
  fits$gamma$summary$median,
  tolerance = 0
)), "coef() must return posterior medians.")

dcustom <- function(x, location, scale, log = FALSE) {
  dnorm(x, location, scale, log = log)
}
pcustom <- function(location, scale, log = FALSE) {
  value <- if (scale > 0) 1 / scale else 0
  if (log) log(value) else value
}
custom_control <- fast$control
custom_control$lower <- c(scale = 0)
custom_control$rng <- function(n, location, scale) rnorm(n, location, scale)
fits$custom <- suppressWarnings(
  run_fit(rnorm(30), dcustom, pcustom,
          start = c(location = 0, scale = 1),
          control = custom_control)
)
check_fit(fits$custom, c("location", "scale"), n_save, fast$chains)

# Public methods and generated quantities.
ci <- confint(fits$gamma, level = 0.90)
assert(identical(dim(ci), c(2L, 2L)), "confint() returned wrong dimensions.")
yrep <- predict(fits$gamma, draws = 20, size = 5, seed = 99)
assert(identical(dim(yrep), c(20L, 5L)) && all(yrep > 0),
       "Posterior prediction failed.")
ll <- log_lik(fits$gamma, draws = 25, seed = 99)
assert(identical(dim(ll), c(25L, 40L)) && all(is.finite(ll)),
       "Pointwise log-likelihood failed.")

# Reproducibility under an explicit seed.
x_rng <- c(0.2, 0.4, 0.7, 1.1, 1.5, 2.0)
fit_rng_a <- run_fit(x_rng, "exponential", "reference")
fit_rng_b <- run_fit(x_rng, "exponential", "reference")
assert(
  identical(fit_rng_a$draws, fit_rng_b$draws),
  "An explicit seed did not reproduce the posterior draws."
)

# Mathematical and input failure tests.
expect_error(run_fit(c(0.2, 0.8), "beta", "MDI"), "improper")
expect_error(run_fit(c(1), "gamma", "Jeffreys"), "n >= 2")
expect_error(run_fit(c(1, 1, 1), "gamma", "Jeffreys"), "nonconstant")
expect_error(run_fit(c(1, 2, 3), "gamma", "reference"), "reference-shape")
expect_error(run_fit(c(1, 2, 3), "weibull", "MDI"), "improper")
expect_error(run_fit(c(1, 2, 3), "negative binomial", "Jeffreys"),
             "fixed")
expect_error(run_fit(c(-1, 0, 2), "Poisson", "Jeffreys"),
             "nonnegative integers")
expect_error(run_fit(c(0, 0), "exponential", "Jeffreys"), "sum")
expect_error(run_fit(c(1, 2, 3), "t", "Jeffreys"), "improper")
expect_error(run_fit(c(1, 2, 3), "lognormal", "MDI"), "improper")
expect_error(
  fitdistrBayes(c(1, NA, 2), "normal", "reference",
                iter = 100, warmup = 50, chains = 2),
  "missing"
)

cat(sprintf(
  "All fitdistrBayes prototype tests passed (%d fitted models; %d failure checks).\n",
  length(fits), 11L
))
