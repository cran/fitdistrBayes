library("fitdistrBayes")

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

close_to <- function(x, y, tolerance = 1e-10) {
  isTRUE(all.equal(unname(x), unname(y), tolerance = tolerance))
}

ctrl <- list(
  rhat_threshold = 2,
  ess_threshold = 1,
  warn_convergence = FALSE
)

fit_start <- function(x, distr, prior, fixed = NULL, start = NULL) {
  fitdistrBayes(
    x, distr, prior, fixed = fixed, start = start,
    iter = 300, warmup = 150, chains = 2, seed = 1701,
    control = ctrl
  )
}

# Beta: ordinary empirical moments (with divisor n).
x_beta <- c(0.08, 0.15, 0.22, 0.31, 0.48, 0.62)
beta_mean <- mean(x_beta)
beta_m2 <- mean((x_beta - beta_mean)^2)
beta_common <- beta_mean * (1 - beta_mean) / beta_m2 - 1
beta_expected <- c(
  shape1 = beta_mean * beta_common,
  shape2 = (1 - beta_mean) * beta_common
)
fit_beta <- fit_start(x_beta, "beta", "Jeffreys")
assert(close_to(fit_beta$initialization$automatic, beta_expected),
       "Beta moment initialization is incorrect.")
assert(identical(fit_beta$initialization$source, "automatic"),
       "Beta automatic initialization was not recorded.")

# User values replace the automatic center and remain visible in the object.
beta_manual <- c(shape1 = 3, shape2 = 7)
fit_beta_manual <- fit_start(
  x_beta, "beta", "Jeffreys", start = beta_manual
)
assert(close_to(fit_beta_manual$initialization$center, beta_manual),
       "User-supplied Beta start was not honored.")
assert(identical(fit_beta_manual$initialization$source, "user-supplied"),
       "User-supplied initialization source was not recorded.")

# Cauchy: median and half interquartile range.
x_cauchy <- c(-5, -2, -1, 0, 0.5, 2, 7)
q_cauchy <- quantile(
  x_cauchy, c(0.25, 0.5, 0.75), names = FALSE, type = 8
)
cauchy_expected <- c(
  location = q_cauchy[2],
  scale = (q_cauchy[3] - q_cauchy[1]) / 2
)
fit_cauchy <- fit_start(x_cauchy, "cauchy", "reference")
assert(close_to(fit_cauchy$initialization$automatic, cauchy_expected),
       "Cauchy quantile initialization is incorrect.")

# Chi-squared and Gamma: ordinary moments.
x_chisq <- c(0.4, 1.1, 2.0, 3.7, 5.2)
fit_chisq <- fit_start(x_chisq, "chi-squared", "reference")
assert(close_to(fit_chisq$initialization$automatic, c(df = mean(x_chisq))),
       "Chi-squared moment initialization is incorrect.")

x_gamma <- c(0.4, 0.7, 1.2, 1.9, 3.1, 5.0)
gamma_mean <- mean(x_gamma)
gamma_m2 <- mean((x_gamma - gamma_mean)^2)
gamma_expected <- c(
  shape = gamma_mean^2 / gamma_m2,
  rate = gamma_mean / gamma_m2
)
fit_gamma <- fit_start(x_gamma, "gamma", "reference-shape")
assert(close_to(fit_gamma$initialization$automatic, gamma_expected),
       "Gamma moment initialization is incorrect.")
assert(identical(fit_gamma$initialization$sampled_parameters, "shape"),
       "Gamma conditional initialization metadata is incorrect.")

# Discrete MDI routes: classical means in the interior.
x_geom <- c(0, 1, 0, 3, 2, 1, 4)
fit_geom <- fit_start(x_geom, "geometric", "MDI")
assert(close_to(
  fit_geom$initialization$automatic,
  c(prob = 1 / (1 + mean(x_geom)))
), "Geometric moment initialization is incorrect.")

x_nb <- c(0, 1, 2, 3, 4, 6, 8)
fit_nb <- fit_start(
  x_nb, "negative binomial", "MDI", fixed = list(size = 4)
)
assert(close_to(fit_nb$initialization$automatic, c(mu = mean(x_nb))),
       "Negative-binomial moment initialization is incorrect.")

x_pois <- c(0, 1, 2, 3, 4, 5)
fit_pois <- fit_start(x_pois, "Poisson", "MDI")
assert(close_to(
  fit_pois$initialization$automatic, c(lambda = mean(x_pois))
), "Poisson moment initialization is incorrect.")

# All-zero samples use an explicit interior continuity correction.
fit_pois_zero <- fit_start(rep(0, 8), "Poisson", "MDI")
assert(close_to(
  fit_pois_zero$initialization$automatic, c(lambda = 0.5 / 8)
), "All-zero Poisson initialization is not interior.")

# Logistic: mean and variance identities.
x_logistic <- c(-3, -1, 0, 0.5, 1.2, 2.7, 4)
logistic_location <- mean(x_logistic)
logistic_scale <- sqrt(
  3 * mean((x_logistic - logistic_location)^2)
) / pi
fit_logistic <- fit_start(x_logistic, "logistic", "reference")
assert(close_to(
  fit_logistic$initialization$automatic,
  c(location = logistic_location, scale = logistic_scale)
), "Logistic moment initialization is incorrect.")

# Student-t with fixed df > 2: mean and variance identities.
x_t <- c(-2, -1, 0, 0.5, 1, 2, 4, 6)
df_fixed <- 6
t_location <- mean(x_t)
t_scale <- sqrt(
  mean((x_t - t_location)^2) * (df_fixed - 2) / df_fixed
)
fit_t <- fit_start(
  x_t, "t", "reference", fixed = list(df = df_fixed)
)
assert(close_to(
  fit_t$initialization$automatic,
  c(location = t_location, scale = t_scale)
), "Fixed-df Student-t moment initialization is incorrect.")

# Weibull: closed-form L-moment estimator.
x_weibull <- c(0.3, 0.6, 0.9, 1.4, 2.1, 3.0, 4.8)
n_weibull <- length(x_weibull)
ordered <- sort(x_weibull)
l1 <- mean(ordered)
b1 <- sum(
  ((seq_len(n_weibull) - 1) / (n_weibull - 1)) * ordered
) / n_weibull
l2 <- 2 * b1 - l1
tau <- l2 / l1
weibull_shape <- -log(2) / log1p(-tau)
weibull_scale <- l1 / gamma(1 + 1 / weibull_shape)
fit_weibull <- fit_start(x_weibull, "weibull", "reference")
assert(close_to(
  fit_weibull$initialization$automatic,
  c(shape = weibull_shape, scale = weibull_scale)
), "Weibull L-moment initialization is incorrect.")
assert(grepl("L-moments", fit_weibull$initialization$method, fixed = TRUE),
       "Weibull initialization method was not recorded.")
assert(identical(fit_weibull$initialization$sampled_parameters, "shape"),
       "Weibull conditional initialization metadata is incorrect.")

# Exact posterior simulation correctly reports that no start is required.
fit_exp <- fit_start(c(0.2, 0.5, 0.9, 1.4), "exponential", "Jeffreys")
assert(identical(fit_exp$initialization$source, "not applicable") &&
       is.null(fit_exp$initialization$center),
       "Exact posterior simulation should not report a starting point.")

cat("All classical-initialization tests passed.\n")
