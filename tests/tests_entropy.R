library("fitdistrBayes")

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

nb_entropy <- getFromNamespace(".fdb_nbinom_entropy", "fitdistrBayes")
merge_control <- getFromNamespace(".fdb_merge_control", "fitdistrBayes")
control <- merge_control(list())

direct_entropy <- function(size, mu, tail = 1e-12) {
  upper <- qnbinom(1 - tail, size = size, mu = mu)
  k <- 0:as.integer(upper)
  log_probability <- dnbinom(k, size = size, mu = mu, log = TRUE)
  -sum(exp(log_probability) * log_probability)
}

# The integral implementation must agree with direct summation over regular,
# overdispersed, and nearly Poisson negative-binomial laws.
grid <- expand.grid(
  size = c(0.2, 0.5, 1, 2, 20),
  mu = c(0.01, 0.1, 1, 10, 100)
)
for (i in seq_len(nrow(grid))) {
  size <- grid$size[i]
  mu <- grid$mu[i]
  integral_value <- nb_entropy(size, mu, control)
  direct_value <- direct_entropy(size, mu)
  assert(
    isTRUE(all.equal(integral_value, direct_value, tolerance = 1e-7)),
    sprintf("Negative-binomial entropy mismatch at size=%g, mu=%g.",
            size, mu)
  )
}

# Extremely dispersed and concentrated cases remain finite without a sum whose
# length grows with the mean.
extreme <- rbind(
  c(size = 1e-6, mu = 1e6),
  c(size = 0.2, mu = 1e6),
  c(size = 1e4, mu = 1e6)
)
extreme_values <- apply(
  extreme, 1, function(z) nb_entropy(z[["size"]], z[["mu"]], control)
)
assert(all(is.finite(extreme_values) & extreme_values >= 0),
       "Extreme negative-binomial entropies were not finite.")

# Regression for the all-zero, small-size MDI case that previously required
# repeated probability sums of up to 200,000 terms.
elapsed <- system.time({
  fit <- fitdistrBayes(
    0, "negative binomial", "MDI", fixed = list(size = 0.2),
    iter = 100, warmup = 50, chains = 2, seed = 81,
    control = list(rhat_threshold = 10, ess_threshold = 1,
                   warn_convergence = FALSE)
  )
})[["elapsed"]]
assert(all(is.finite(fit$draws$mu) & fit$draws$mu > 0),
       "Small-size negative-binomial MDI fit produced invalid draws.")
assert(elapsed < 10,
       sprintf("Small-size negative-binomial MDI fit is too slow: %.2fs.",
               elapsed))

cat(sprintf(
  "Negative-binomial entropy tests passed; performance regression: %.3fs.\n",
  elapsed
))
