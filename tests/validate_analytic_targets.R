library("fitdistrBayes")

assert_mc <- function(estimate, target, mcse, label, multiplier = 8) {
  error <- abs(estimate - target)
  tolerance <- multiplier * mcse + 1e-10 * max(1, abs(target))
  if (!is.finite(error) || error > tolerance) {
    stop(sprintf(
      "%s failed: estimate=%g, target=%g, error=%g, tolerance=%g.",
      label, estimate, target, error, tolerance
    ), call. = FALSE)
  }
}

settings <- list(iter = 6000L, warmup = 1000L, chains = 4L)

# Exponential: rate | x ~ Gamma(n, sum(x)).
x <- c(0.2, 0.4, 0.8, 1.1, 1.6, 2.2)
fit <- fitdistrBayes(
  x, "exponential", "Jeffreys", seed = 1,
  iter = settings$iter, warmup = settings$warmup,
  chains = settings$chains
)
assert_mc(fit$summary$mean, length(x) / sum(x),
          fit$summary$mcse_mean, "Exponential posterior mean")

# Poisson: lambda | x ~ Gamma(sum(x)+1/2, n).
x <- c(0, 1, 2, 3, 4, 2, 1, 5)
fit <- fitdistrBayes(
  x, "Poisson", "Jeffreys", seed = 2,
  iter = settings$iter, warmup = settings$warmup,
  chains = settings$chains
)
target <- (sum(x) + 0.5) / length(x)
assert_mc(fit$summary$mean, target, fit$summary$mcse_mean,
          "Poisson posterior mean")

# Geometric: p | x ~ Beta(n, sum(x)+1/2).
x <- c(0, 1, 0, 3, 2, 1, 4, 0)
fit <- fitdistrBayes(
  x, "geometric", "reference", seed = 3,
  iter = settings$iter, warmup = settings$warmup,
  chains = settings$chains
)
target <- length(x) / (length(x) + sum(x) + 0.5)
assert_mc(fit$summary$mean, target, fit$summary$mcse_mean,
          "Geometric posterior mean")

# Fixed-size negative binomial: p | x ~ Beta(nr, sum(x)+1/2),
# mu = r(1-p)/p, E(mu)=r*b/(a-1).
x <- c(0, 1, 3, 2, 4, 5, 1, 2)
r <- 4
a <- length(x) * r
b <- sum(x) + 0.5
fit <- fitdistrBayes(
  x, "negative binomial", "Jeffreys", fixed = list(size = r),
  seed = 4, iter = settings$iter, warmup = settings$warmup,
  chains = settings$chains
)
target <- r * b / (a - 1)
assert_mc(fit$summary$mean, target, fit$summary$mcse_mean,
          "Negative-binomial posterior mean")

# Normal/reference: E(mu|x)=xbar and sigma^2 has an inverse-Gamma posterior.
x <- c(-1.2, -0.4, 0.1, 0.7, 1.3, 2.1, 2.4)
fit <- fitdistrBayes(
  x, "normal", "reference", seed = 5,
  iter = settings$iter, warmup = settings$warmup,
  chains = settings$chains
)
row_mu <- fit$summary[fit$summary$parameter == "mean", ]
assert_mc(row_mu$mean, mean(x), row_mu$mcse_mean,
          "Normal posterior location mean")
shape_s2 <- (length(x) - 1) / 2
Q <- sum((x - mean(x))^2)
target_s2 <- (Q / 2) / (shape_s2 - 1)
s2 <- fit$draws$sd^2
mcse_s2 <- sd(s2) / sqrt(nrow(fit$draws))
assert_mc(mean(s2), target_s2, mcse_s2,
          "Normal posterior variance mean")

# Gamma marginal/conditional identity:
# E(rate | x) = n E(shape | x) / sum(x).
x <- c(0.4, 0.7, 0.9, 1.1, 1.5, 2.0, 2.6, 3.2)
fit <- fitdistrBayes(
  x, "gamma", "reference-shape", seed = 6,
  iter = settings$iter, warmup = settings$warmup,
  chains = settings$chains
)
target_rate <- length(x) * mean(fit$draws$shape) / sum(x)
row_rate <- fit$summary[fit$summary$parameter == "rate", ]
assert_mc(row_rate$mean, target_rate, row_rate$mcse_mean,
          "Gamma conditional-rate identity")

# Weibull conditional identity: z=scale^(-shape) and, conditionally on shape,
# z*A(shape) follows Gamma(n, 1).  Its unconditional mean must be n.
x <- c(0.5, 0.8, 1.0, 1.4, 1.9, 2.5, 3.1, 4.0)
fit <- fitdistrBayes(
  x, "weibull", "reference", seed = 7,
  iter = settings$iter, warmup = settings$warmup,
  chains = settings$chains
)
shape <- fit$draws$shape
scale <- fit$draws$scale
conditional_gamma <- vapply(seq_along(shape), function(i) {
  sum(x^shape[i]) * scale[i]^(-shape[i])
}, numeric(1))
mcse_identity <- sd(conditional_gamma) / sqrt(length(conditional_gamma))
assert_mc(mean(conditional_gamma), length(x), mcse_identity,
          "Weibull conditional-Gamma identity")

cat("All analytical posterior-target validations passed.\n")
