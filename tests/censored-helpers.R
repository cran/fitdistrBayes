library(fitdistrBayes)
E <- asNamespace("fitdistrBayes")
control_quiet <- list(warn_convergence = FALSE)
fixtures <- list(
  beta = c(shape1 = 2, shape2 = 5),
  cauchy = c(location = 1, scale = 0.7),
  `chi-squared` = c(df = 4),
  exponential = c(rate = 0.7),
  gamma = c(shape = 2.5, rate = 1.3),
  geometric = c(prob = 0.3),
  lognormal = c(meanlog = 0.3, sdlog = 0.7),
  logistic = c(location = 1, scale = 0.7),
  `negative binomial` = c(mu = 5),
  normal = c(mean = 1, sd = 0.7),
  Poisson = c(lambda = 5),
  t = c(location = 1, scale = 0.7, df = 5),
  weibull = c(shape = 1.7, scale = 2),
  frechet = c(shape = 3, scale = 2),
  gumbel = c(location = 1, scale = 0.7),
  lomax = c(shape = 3, scale = 2),
  nakagami = c(shape = 2, spread = 3),
  `exponential-logarithmic` = c(theta = 0.5, rate = 1),
  rician = c(noncentrality = 2, scale = 0.7),
  `weighted lindley` = c(lambda = 1.5, phi = 2))

get_fixed <- function(model, prior) {
  if (model == "negative binomial") return(list(size = 2))
  if (model == "t" && prior != "independence-jeffreys") return(list(df = 5))
  NULL
}

simulate_model <- function(model, n = 100, seed = 11, theta = fixtures[[model]]) {
  set.seed(seed)
  prior <- E$fitdistrBayes_routes(model)$prior[1]
  fixed <- get_fixed(model, prior)
  dummy <- if (model == "beta") seq(0.1, 0.9, length.out = 20) else seq_len(20)
  built <- E$.fdb_build_builtin(dummy, model, prior, fixed)
  built$rng(theta, n)
}

censor_data <- function(y, fraction = 0.35) {
  limit <- as.numeric(stats::quantile(y, 1 - fraction, type = 1))
  list(x = pmin(y, limit), status = as.integer(y <= limit), limit = limit)
}

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
expect_error <- function(expr, pattern = NULL) {
  msg <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  assert(!is.null(msg), "Expected an error but the call succeeded.")
  if (!is.null(pattern)) assert(grepl(pattern, msg, ignore.case = TRUE), paste("Unexpected error:", msg))
  invisible(msg)
}
