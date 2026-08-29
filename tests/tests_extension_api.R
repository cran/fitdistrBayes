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

control_fast <- list(
  rhat_threshold = 1.25, ess_threshold = 20,
  warn_convergence = FALSE
)

# A new distribution with a proper non-objective prior and automatic MCMC.
d_laplace <- function(x, location, scale, log = FALSE) {
  value <- -log(2 * scale) - abs(x - location) / scale
  if (log) value else exp(value)
}
p_laplace <- function(location, scale, log = FALSE) {
  value <- dnorm(location, 0, 5, log = TRUE) +
    dlnorm(scale, 0, 0.75, log = TRUE)
  if (log) value else exp(value)
}
r_laplace <- function(n, location, scale) {
  location + scale * ifelse(runif(n) < 0.5, -1, 1) * rexp(n)
}
laplace_moments <- data.frame(
  parameter = c("location", "scale"),
  mean_exists = c(TRUE, TRUE),
  variance_exists = c(TRUE, TRUE),
  note = rep("finite under the stated proper Normal--Lognormal prior", 2L),
  stringsAsFactors = FALSE
)
laplace_model <- fitdistrBayes_model(
  density = d_laplace,
  prior = p_laplace,
  start = c(location = 0, scale = 1),
  lower = c(scale = 0),
  name = "Laplace",
  engine = "adaptive_metropolis",
  propriety = "proper Normal--Lognormal prior and a nonconstant sample",
  moments = laplace_moments,
  rng = r_laplace,
  validate = function(x) length(x) >= 2L && diff(range(x)) > 0,
  prior_label = "Normal--Lognormal",
  prior_kernel = "Normal(location; 0, 5^2) Lognormal(scale; 0, 0.75^2)"
)
set.seed(1)
x_laplace <- r_laplace(35, 1, 1.5)
fit_laplace <- fitdistrBayes(
  x_laplace, laplace_model,
  iter = 500, warmup = 200, chains = 2, seed = 2,
  control = control_fast
)
assert(inherits(fit_laplace, "fitdistrBayes"),
       "The model specification did not produce a fitted object.")
assert(identical(fit_laplace$engine$algorithm,
                 "generic adaptive Metropolis"),
       "The requested adaptive Metropolis engine was not used.")
assert(isTRUE(fit_laplace$model$support_checked),
       "The user data validator was not recorded.")
assert(all(fit_laplace$summary$mean_exists) &&
       all(fit_laplace$summary$variance_exists),
       "User-supplied moment information was not propagated.")
assert(all(dim(predict(fit_laplace, draws = 5, size = 4, seed = 3)) ==
           c(5, 4)), "The extension predictive generator failed.")

# A one-parameter extension using the package's slice sampler.
d_exponential <- function(x, rate, log = FALSE) {
  dexp(x, rate = rate, log = log)
}
p_rate <- function(rate, log = FALSE) {
  dgamma(rate, shape = 2, rate = 1, log = log)
}
rate_model <- fitdistrBayes_model(
  d_exponential, p_rate, start = c(rate = 1), lower = 0,
  name = "Exponential with Gamma prior", engine = "slice",
  propriety = TRUE,
  moments = data.frame(
    parameter = "rate", mean_exists = TRUE, variance_exists = TRUE,
    note = "Gamma posterior", stringsAsFactors = FALSE
  ),
  validate = function(x) all(x >= 0)
)
fit_rate <- fitdistrBayes(
  rexp(30, 2), rate_model,
  iter = 400, warmup = 100, chains = 2, seed = 4,
  control = control_fast
)
assert(identical(fit_rate$engine$algorithm,
                 "generic univariate slice sampling"),
       "The requested slice engine was not used.")

# A proper informative prior with exact conjugate posterior simulation.
d_poisson <- function(x, lambda, log = FALSE) {
  dpois(x, lambda = lambda, log = log)
}
p_lambda <- function(lambda, log = FALSE) {
  dgamma(lambda, shape = 2, rate = 1, log = log)
}
sample_poisson_gamma <- function(x, n_save, chains, ...) {
  lapply(seq_len(chains), function(chain) {
    answer <- matrix(
      rgamma(n_save, shape = sum(x) + 2, rate = length(x) + 1),
      ncol = 1L
    )
    colnames(answer) <- "lambda"
    answer
  })
}
poisson_gamma <- fitdistrBayes_model(
  d_poisson, p_lambda, start = c(lambda = 1), lower = 0,
  name = "Poisson--Gamma", engine = "custom",
  sampler = sample_poisson_gamma, independent = TRUE,
  engine_label = "exact Gamma posterior",
  propriety = function(x, fixed) list(
    proper = length(x) >= 1L && all(x >= 0) && all(x == floor(x)),
    message = "the Gamma posterior has positive shape and rate"
  ),
  moments = data.frame(
    parameter = "lambda", mean_exists = TRUE, variance_exists = TRUE,
    note = "all positive moments of the Gamma posterior are finite",
    stringsAsFactors = FALSE
  ),
  validate = function(x) all(x >= 0) && all(x == floor(x)),
  prior_label = "Gamma(2, 1)"
)
x_count <- c(2, 4, 3, 5, 1, 3)
fit_exact <- fitdistrBayes(
  x_count, poisson_gamma,
  iter = 220, warmup = 20, chains = 2, seed = 5
)
target <- (sum(x_count) + 2) / (length(x_count) + 1)
assert(fit_exact$diagnostics$exact_or_independent,
       "Independent exact simulation was not recorded.")
assert(identical(fit_exact$engine$algorithm, "exact Gamma posterior"),
       "The custom engine label was not retained.")
assert(abs(fit_exact$summary$mean - target) <
         7 * fit_exact$summary$mcse_mean,
       "Exact custom draws disagree with the analytical posterior mean.")

# Invalid user declarations must stop before sampling.
bad_propriety <- fitdistrBayes_model(
  d_exponential, p_rate, start = c(rate = 1), lower = 0,
  engine = "slice", propriety = FALSE
)
expect_error(
  fitdistrBayes(rexp(5), bad_propriety,
                iter = 100, warmup = 50, chains = 1),
  "propriety"
)
expect_error(
  fitdistrBayes(rep(1, 3), laplace_model,
                iter = 100, warmup = 50, chains = 1),
  "validator"
)
expect_error(
  fitdistrBayes(x_count, poisson_gamma, prior = "jeffreys",
                iter = 100, warmup = 50, chains = 1),
  "Do not supply"
)

cat("All extension-interface tests passed.\n")
