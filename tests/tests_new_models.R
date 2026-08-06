library("fitdistrBayes")

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

# Each newly registered density integrates to one under the package
# parameterization.
gumbel_density <- function(x, location, scale) {
  z <- (x - location) / scale
  exp(-z - exp(-z)) / scale
}
frechet_density <- function(x, shape, scale) {
  scale * shape * x^(-shape - 1) * exp(-scale * x^(-shape))
}
lomax_density <- function(x, shape, scale) {
  shape / scale * (1 + x / scale)^(-shape - 1)
}
nakagami_density <- function(x, shape, spread) {
  2 / gamma(shape) * (shape / spread)^shape * x^(2 * shape - 1) *
    exp(-shape * x^2 / spread)
}
el_density <- function(x, theta, rate) {
  -rate * (1 - theta) * exp(-rate * x) /
    (log(theta) * (1 - (1 - theta) * exp(-rate * x)))
}
rician_density <- function(x, noncentrality, scale) {
  argument <- noncentrality * x / scale^2
  x / scale^2 * besselI(argument, 0, expon.scaled = TRUE) *
    exp(argument - (x^2 + noncentrality^2) / (2 * scale^2))
}

checks <- c(
  integrate(gumbel_density, -Inf, Inf, location = 1, scale = 2)$value,
  integrate(frechet_density, 0, Inf, shape = 2.5, scale = 4)$value,
  integrate(lomax_density, 0, Inf, shape = 3, scale = 2)$value,
  integrate(nakagami_density, 0, Inf, shape = 2.5, spread = 4)$value,
  integrate(el_density, 0, Inf, theta = 0.4, rate = 1.3)$value,
  integrate(rician_density, 0, Inf, noncentrality = 5, scale = 2)$value
)
stopifnot(max(abs(checks - 1)) < 2e-7)

# The accelerated Rician Jeffreys factor is checked against the defining
# integral rather than against the implementation that it replaced.
log_q_factory <- getFromNamespace(".fdb_rician_logq_factory", "fitdistrBayes")
log_q <- log_q_factory()
q_direct <- function(rho) {
  integrand <- function(y) {
    i0e <- besselI(y, 0, expon.scaled = TRUE)
    i1e <- besselI(y, 1, expon.scaled = TRUE)
    y^3 / rho^2 * exp(-(y - rho)^2 / (2 * rho)) * i1e^2 / i0e
  }
  integral <- integrate(integrand, 0, Inf, rel.tol = 1e-10,
                        subdivisions = 2000L)$value
  psi <- integral - rho
  (rho + 1) * psi - rho
}
rho_grid <- c(1e-3, 0.1, 1, 10, 30, 100)
q_fast <- exp(vapply(log(rho_grid), log_q, numeric(1)))
q_truth <- vapply(rho_grid, q_direct, numeric(1))
stopifnot(max(abs(q_fast - q_truth) / q_truth) < 7e-4)

# Automatic starts expose the classical construction used by each model.
ctrl <- list(rhat_threshold = 5, ess_threshold = 1,
             warn_convergence = FALSE)
set.seed(20260719)
x_gumbel <- 1 - 2 * log(-log(runif(30)))
x_frechet <- (4 / rexp(30))^(1 / 2.5)
x_lomax <- 3 * (runif(30)^(-1 / 3) - 1)
x_nakagami <- sqrt(rgamma(30, 2.5, rate = 2.5 / 4))
x_rician <- sqrt(rnorm(30, 5, 2)^2 + rnorm(30, 0, 2)^2)
theta_el <- 0.4
rate_el <- 1.3
u_el <- runif(30)
x_el <- -(log(-expm1((1 - u_el) * log(theta_el))) -
            log1p(-theta_el)) / rate_el

data_list <- list(
  gumbel = x_gumbel, frechet = x_frechet, lomax = x_lomax,
  nakagami = x_nakagami, rician = x_rician,
  `exponential-logarithmic` = x_el
)
prior_list <- c(
  gumbel = "jeffreys", frechet = "jeffreys", lomax = "jeffreys",
  nakagami = "reference", rician = "jeffreys",
  `exponential-logarithmic` = "mdi"
)
method_patterns <- c(
  gumbel = "method of moments", frechet = "L-moments",
  lomax = "L-moments", nakagami = "method of moments",
  rician = "method of moments",
  `exponential-logarithmic` = "moments"
)

for (model in names(data_list)) {
  fit <- fitdistrBayes(
    data_list[[model]], model, prior_list[[model]],
    iter = 240, warmup = 120, chains = 2, seed = 200 + match(model, names(data_list)),
    control = ctrl
  )
  stopifnot(
    identical(fit$initialization$source, "automatic"),
    grepl(method_patterns[[model]], fit$initialization$method,
          ignore.case = TRUE),
    all(is.finite(as.matrix(fit$draws[, fit$model$parameters, drop = FALSE]))),
    all(is.finite(log_lik(fit, draws = 8))),
    all(is.finite(predict(fit, draws = 8, seed = 99)))
  )
}

# Invalid or unproved routes fail before sampling.
expect_error(
  fitdistrBayes(x_lomax, "lomax", "reference",
                iter = 40, warmup = 20, chains = 1),
  "improper"
)
expect_error(
  fitdistrBayes(x_nakagami, "nakagami", "mdi",
                iter = 40, warmup = 20, chains = 1),
  "improper"
)
expect_error(
  fitdistrBayes(x_rician, "rician", "reference",
                iter = 40, warmup = 20, chains = 1),
  "not available"
)
expect_error(
  fitdistrBayes(x_el[1:2], "el", "jeffreys",
                iter = 40, warmup = 20, chains = 1),
  "n > 2"
)

cat("New-model analytic, numerical, initialization, and route tests passed.\n")
