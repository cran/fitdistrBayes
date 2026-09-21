library(fitdistrBayes)
local({
  E <- asNamespace("fitdistrBayes")
  err <- function(expr, pattern) {
    msg <- tryCatch({ force(expr); "NO ERROR" }, error = conditionMessage)
    stopifnot(grepl(pattern, msg, ignore.case = TRUE))
  }
  set.seed(321)
  x <- rexp(60, 0.7)
  a <- fitdistrBayes(x, "exponential", "reference", seed = 42,
                    iter = 1200, warmup = 600, chains = 2)
  state <- .Random.seed
  stopifnot(is.null(a$criteria))
  b <- suppressWarnings(fitdistrBayes(x, "exponential", "reference", seed = 42,
       iter = 1200, warmup = 600, chains = 2, criteria = c("waic", "dic")))
  stopifnot(identical(a$chains, b$chains), identical(a$summary, b$summary),
            identical(state, .Random.seed))
  scores <- suppressWarnings(criteria(a, c("waic", "dic")))
  stopifnot(identical(state, .Random.seed), all(scores$estimates$available),
            identical(scores, b$criteria), nrow(as.data.frame(scores)) == 2L)
  ll <- log_lik(a)
  lambda <- a$draws$rate
  dic <- 2 * mean(-2 * rowSums(ll)) + 2 * sum(dexp(x, mean(lambda), log = TRUE))
  stopifnot(abs(dic - DIC(a)$estimates$estimate) < 1e-10)
  # Analytic exponential E[D] and D(E[lambda]); modest Monte Carlo tolerance.
  shape <- length(x); rate <- sum(x)
  analytic <- 2 * (-2 * (shape * (digamma(shape) - log(rate)) - shape)) +
    2 * sum(dexp(x, shape / rate, log = TRUE))
  stopifnot(abs(dic - analytic) < 0.5)
  # One likelihood evaluation per draw, shared between criteria, plus one for DIC.
  counted <- a; count <- 0L; callback <- a$.loglik
  counted$.loglik <- function(theta) { count <<- count + 1L; callback(theta) }
  suppressWarnings(criteria(counted, c("waic", "dic")))
  stopifnot(count == nrow(a$draws) + 1L)
  # Censoring must be evaluated as log S, not log f of an imputed lifetime.
  status <- as.integer(x <= 1.8); observed <- pmin(x, 1.8)
  c0 <- fitcensBayes(observed, status, "exponential", "reference", seed = 43,
                    iter = 1200, warmup = 600, chains = 2)
  state <- .Random.seed
  c1 <- suppressWarnings(fitcensBayes(observed, status, "exponential", "reference", seed = 43,
       iter = 1200, warmup = 600, chains = 2, criteria = c("waic", "dic")))
  stopifnot(identical(c0$chains, c1$chains), identical(state, .Random.seed))
  lc <- log_lik(c0)
  expected <- outer(log(c0$draws$rate), status) - outer(c0$draws$rate, observed)
  stopifnot(max(abs(lc - expected)) < 1e-12)
  # Invalid options fail before any sampling/RNG change.
  state <- .Random.seed
  for (invalid in list(NA, c(TRUE, FALSE), "AIC", list("waic"), character()))
    err(fitdistrBayes(x, "exponential", "reference", criteria = invalid), "criter|Use FALSE")
  err(fitdistrBayes(x, "exponential", "reference", criteria = "dic",
                  control = list(store_callables = FALSE)), "store_callables")
  stopifnot(identical(state, .Random.seed))
  # Missing/infinite means must not be replaced by finite empirical averages.
  badmean <- a; badmean$moment_status$mean_exists <- FALSE
  badmean$.loglik <- function(theta) stop("Must not be evaluated")
  dc <- suppressWarnings(DIC(badmean))
  stopifnot(is.na(dc$estimates$estimate), !dc$estimates$available)
  badll <- a; badll$.loglik <- function(theta) rep(-Inf, length(x))
  err(WAIC(badll), "finite")
  err(compare_models(a, c0, criterion = "waic"), "identical retained data")
  reordered <- a; reordered$data <- rev(x)
  err(compare_models(a, reordered, criterion = "waic"), "identical retained data")
  changed <- a; changed$model$name <- "Poisson"
  err(compare_models(a, changed, criterion = "waic"), "likelihood measure")
  omitted <- a; omitted$omitted <- 3L
  err(compare_models(a, omitted, criterion = "waic"), "omitted rows")
  cc <- suppressWarnings(compare_models(A = scores, B = scores, criterion = "waic"))
  stopifnot(all(cc$table$delta == 0), all(cc$table$se_delta == 0))
  # Guard corner cases, including ties and held-out single events.
  guard <- function(model, z, prior = "reference", status = rep(1L, length(z)), fixed = list()) {
    o <- list(model = list(name = model, fixed = fixed), data = z, status = status,
              prior = list(key = prior))
    class(o) <- "fitcensBayes"
    E$.fdb_loo_guard(o)$certified
  }
  stopifnot(identical(guard("normal", c(1, 1, 2)), c(TRUE, TRUE, FALSE)),
    !any(guard("normal", c(1, 2))), all(guard("normal", c(1, 2, 3))),
    identical(guard("exponential", c(1, 2, 3), status = c(1, 0, 0)), c(FALSE, TRUE, TRUE)),
    all(guard("exponential", c(1, 2, 3), "mdi", c(0, 0, 0))),
    all(guard("Poisson", c(0, 0, 0))),
    !any(guard("rician", c(1, 2, 3))),
    all(guard("rician", c(1, 2, 3, 4))),
    !any(guard("t", c(1, 2), "independence-jeffreys")),
    all(guard("t", c(1, 2, 3), "independence-jeffreys")),
    !any(guard("t", c(1, 2, 3), "jeffreys", fixed = list(df = 1))))
  # Compare native loo implementations when available; no optional dependency
  # is necessary for the preceding WAIC and DIC tests.
  if (requireNamespace("loo", quietly = TRUE)) {
    allcrit <- suppressWarnings(criteria(a))
    native_waic <- suppressWarnings(loo::waic(ll))
    native_loo <- suppressWarnings(loo::loo(ll, r_eff = 1, cores = 1))
    stopifnot(abs(allcrit$estimates$estimate[1] - native_waic$estimates["waic", "Estimate"]) < 1e-10,
      abs(allcrit$estimates$estimate[2] - native_loo$estimates["looic", "Estimate"]) < 1e-10)
    # Exact LOO predictive densities from conjugate exponential posteriors.
    exact_loo <- lgamma(shape) - lgamma(shape - 1) +
      (shape - 1) * log(rate - x) - shape * log(rate)
    stopifnot(abs(sum(exact_loo) - allcrit$details$looic$estimates["elpd_loo", "Estimate"]) < 0.5)
    cs <- suppressWarnings(LOOIC(c0))
    sh <- sum(status); rt <- sum(observed)
    exact_c <- lgamma(sh) - lgamma(sh - status) +
      (sh - status) * log(rt - observed) - sh * log(rt)
    stopifnot(abs(sum(exact_c) - cs$estimates$elpd) < 0.5)
    tiny <- fitdistrBayes(c(1, 2), "normal", "reference", iter = 100, warmup = 50,
                         chains = 2, seed = 15)
    stopifnot(!suppressWarnings(LOOIC(tiny))$estimates$available)
    mcmc <- fitdistrBayes(x, "gamma", "reference-shape", iter = 1000, warmup = 500,
      chains = 2, seed = 91, control = list(warn_convergence = FALSE))
    mc <- suppressWarnings(LOOIC(mcmc)); ml <- log_lik(mcmc)
    re <- loo::relative_eff(exp(sweep(ml, 2, apply(ml, 2, max), "-")),
                            chain_id = mcmc$draws$.chain, cores = 1)
    stopifnot(identical(unname(mc$details$loo_r_eff), unname(re)))
    mdi <- fitdistrBayes(x, "exponential", "mdi", seed = 24,
                         iter = 1200, warmup = 600, chains = 2)
    comp <- suppressWarnings(compare_models(A = a, B = mdi))
    pw <- lapply(comp$results, function(z) z$details$looic$pointwise[, "looic"])
    stopifnot(abs(comp$table$se_delta[2] - sqrt(length(x) * var(pw[[1]] - pw[[2]]))) < 1e-10)
  }
  cat("PASS: criteria contracts, likelihood identity, RNG/chains, analytic targets, diagnostics, and comparison checks.\n")
})
