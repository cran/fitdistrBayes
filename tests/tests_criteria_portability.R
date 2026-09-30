library(fitdistrBayes)
local({
  E <- asNamespace("fitdistrBayes")
  close <- function(x, y, tol = 1e-12) {
    stopifnot(isTRUE(all.equal(x, y, tolerance = tol, check.attributes = FALSE)))
  }

  # Constant columns must have zero variance even without long-double sums.
  for (v in c(0, 1e-300, -1e-300, 1e307, -1e307,
              .Machine$double.xmax, -.Machine$double.xmax)) {
    stopifnot(identical(E$.fdb_ic_variance(rep(v, 500)), 0),
              identical(E$.fdb_ic_se(rep(v, 500)), 0))
  }
  for (factor in c(1e-150, 1, 1e150)) {
    z <- c(-2, -1, 0, 1, 2) * factor
    close(E$.fdb_ic_variance(z) / (2.5 * factor^2), 1)
  }
  for (factor in c(1e-200, 1e200)) {
    close(E$.fdb_ic_se(c(1, 2, 3) * factor) / (sqrt(3) * factor), 1)
  }
  # Large offsets with much smaller (but representable) deviations.
  z <- 1e150 + c(-2, -1, 0, 1, 2) * 1e135
  centered <- z - z[1L]
  expected <- var(centered / 1e135) * 1e270
  close(E$.fdb_ic_variance(z) / expected, 1)
  stopifnot(is.infinite(E$.fdb_ic_variance(c(-1e200, 1e200))),
            is.infinite(E$.fdb_ic_variance(c(-1e308, 1e308))),
            is.na(E$.fdb_ic_se(1)), is.na(E$.fdb_ic_se(c(1, Inf))))

  # Emulate a variance routine restricted to double-precision accumulation.
  # The copy has a private environment; neither the namespace nor GlobalEnv
  # is modified. Also assert that every input reaching var() is bounded.
  emulation <- new.env(parent = E)
  emulation$bounded_var <- function(x) {
    stopifnot(all(is.finite(x)), all(abs(x) <= 1))
    accumulate <- function(y) Reduce(`+`, y, init = 0)
    avg <- accumulate(x) / length(x)
    centered <- x - avg
    corrected <- accumulate(centered^2) - accumulate(centered)^2 / length(x)
    max(0, corrected) / (length(x) - 1L)
  }
  replace_var <- function(expr) {
    if (identical(expr, quote(stats::var))) return(quote(bounded_var))
    if (is.call(expr)) return(as.call(lapply(as.list(expr), replace_var)))
    expr
  }
  for (nm in c(".fdb_ic_sd", ".fdb_ic_variance", ".fdb_ic_se",
               ".fdb_compute_criteria")) {
    f <- E[[nm]]
    if (nm == ".fdb_ic_sd") body(f) <- replace_var(body(f))
    environment(f) <- emulation
    emulation[[nm]] <- f
  }

  set.seed(1909)
  fit <- fitdistrBayes(1:20, "exponential", "reference",
    iter = 500, warmup = 250, chains = 2, seed = 9)
  original_chains <- fit$chains
  state <- .Random.seed
  normal <- WAIC(fit)
  emulated <- emulation$.fdb_compute_criteria(fit, "waic")
  close(normal$estimates$estimate, emulated$estimates$estimate)
  close(normal$details$waic$pointwise, emulated$details$waic$pointwise)

  unavailable <- function(result, kind) {
    stopifnot(inherits(result, "fitdistrBayes_criteria"),
      !result$estimates$available, !result$estimates$reliable,
      all(is.na(result$estimates[c("estimate", "se", "p_eff", "elpd")])),
      grepl(kind, result$estimates$reason, fixed = TRUE))
  }
  big <- fit
  big$.loglik <- function(theta) rep(-1e307, 20)
  for (compute in list(function(x) WAIC(x), function(x)
                      emulation$.fdb_compute_criteria(x, "waic"))) {
    result <- suppressWarnings(compute(big))
    unavailable(result, "aggregate")
    stopifnot(all(result$details$waic$pointwise[, "p_waic"] == 0),
              all(is.finite(result$details$waic$pointwise)))
  }
  big$.loglik <- function(theta) rep(-1e308, 20)
  result <- suppressWarnings(WAIC(big))
  unavailable(result, "pointwise")
  stopifnot(identical(unname(result$details$waic$unrepresentable_observations), 1:20),
            all(is.na(result$details$waic$pointwise[, "waic"])))

  # True variance overflow is not confused with a harmless constant offset.
  cutoff <- median(fit$draws$rate)
  big$.loglik <- function(theta) rep(if (theta[["rate"]] < cutoff) -1e200 else 1e200, 20)
  result <- suppressWarnings(WAIC(big))
  unavailable(result, "pointwise")
  stopifnot(all(is.na(result$details$waic$pointwise)))
  stopifnot(identical(state, .Random.seed), identical(original_chains, fit$chains))

  # Agreement with the independently maintained implementation on ordinary
  # complete and right-censored observed-data likelihoods.
  cens <- fitcensBayes(pmin(1:20, 12), as.integer(1:20 <= 12),
    "exponential", "reference", iter = 500, warmup = 250, chains = 2, seed = 10)
  for (object in list(fit, cens)) {
    actual <- suppressWarnings(WAIC(object))
    ll <- log_lik(object)
    penalty <- apply(ll, 2L, var)
    maxima <- apply(ll, 2L, max)
    lppd <- maxima + log(colMeans(exp(sweep(ll, 2L, maxima, "-"))))
    close(actual$estimates$estimate, -2 * sum(lppd - penalty))
    close(actual$estimates$p_eff, sum(penalty))
    if (requireNamespace("loo", quietly = TRUE)) {
      reference <- suppressWarnings(loo::waic(ll))
      close(actual$details$waic$pointwise, reference$pointwise)
      close(actual$estimates$se, reference$estimates["waic", "SE"])
    }
  }
  cat("PASS: WAIC portability, double-only arithmetic, genuine overflow, and reference agreement.\n")
})
