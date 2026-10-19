source("censored-helpers.R", local = TRUE)
set.seed(90515)
counter <- new.env(parent = emptyenv())
counter$checks <- 0L
check <- function(value, message) {
  assert(value, message)
  counter$checks <- counter$checks + 1L
}
routes <- fitcensBayes_models()
check(nrow(routes) == 59 && length(unique(routes$model)) == 21, "Catalog changed.")
check(identical(names(formals(fitcensBayes))[1:4], c("x", "status", "distr", "prior")), "Argument order.")

# Every density and prior agrees with the frozen original at multiple points.
prior_checks <- 0L
for (r in seq_len(nrow(routes))) {
  model <- E$.fdb_model_name(routes$model[r]); prior <- routes$prior[r]
  fixed <- get_fixed(model, prior)
  y <- simulate_model(model, 70, seed = 10 + r)
  if (model == "t" && prior == "independence-jeffreys") fixed <- NULL
  built <- E$.fdb_build_builtin(y, model, prior, fixed)
  th <- fixtures[[model]][built$parameters]
  ours <- E$.cens_density(model, y, th, fixed)
  check(isTRUE(all.equal(unname(ours), unname(built$loglik(th)), tolerance = 1e-10)),
        paste("Density mismatch", model, prior))
  env <- environment(built$loglik)
  if (exists("logpost", env, inherits = FALSE) && exists("init", env, inherits = FALSE) &&
      exists("transform", env, inherits = FALSE) &&
      !model %in% c("chi-squared", "geometric", "weibull", "frechet", "nakagami")) {
    old <- get("logpost", env); transform <- get("transform", env); u <- unname(get("init", env))
    if (model == "exponential-logarithmic") u <- c(log(-log(0.5)), 0)
    logjac <- function(v) {
      h <- 1e-5
      J <- sapply(seq_along(v), function(j) {
        l <- v; r <- v; l[j] <- l[j] - h; r[j] <- r[j] + h
        (transform(r) - transform(l)) / (2 * h)
      })
      as.numeric(determinant(matrix(J, length(v)), logarithm = TRUE)$modulus)
    }
    rq <- if (model == "rician") E$.fdb_rician_logq_factory() else NULL
    newer <- function(v) {
      th <- transform(v); names(th) <- built$parameters
      sum(E$.cens_density(model, y, th, fixed)) +
        E$.cens_prior(model, prior, th, fixed, E$.cens_control(list()), rq) + logjac(v)
    }
    v <- u + seq_along(u) * 0.04
    check(abs((newer(v) - newer(u)) - (old(v) - old(u))) < 1e-5,
          paste("Prior/coordinate mismatch", model, prior))
    prior_checks <- prior_checks + 1L
  }
}

# Survival functions independently checked by differentiation/integration.
for (model in names(fixtures)) {
  th <- fixtures[[model]]; fixed <- get_fixed(model, if (model == "t") "independence-jeffreys" else "jeffreys")
  y <- simulate_model(model, 200, 711)
  grid <- as.numeric(stats::quantile(y, c(0.2, 0.5, 0.8)))
  s <- E$.cens_survival(model, grid, th, fixed)
  check(all(is.finite(s)) && all(diff(s) <= 1e-12) && all(s <= 0), paste("Survival", model))
  if (model %in% c("geometric", "Poisson", "negative binomial")) {
    p <- exp(E$.cens_survival(model, grid - 1, th, fixed)) - exp(s)
    check(max(abs(p - exp(E$.cens_density(model, grid, th, fixed)))) < 1e-10, paste("Discrete tail", model))
  } else {
    h <- 1e-5 * max(0.01, stats::sd(y))
    delta <- (exp(E$.cens_survival(model, grid - h, th, fixed)) -
              exp(E$.cens_survival(model, grid + h, th, fixed))) / (2 * h)
    d <- exp(E$.cens_density(model, grid, th, fixed))
    check(max(abs(delta / d - 1)) < 2e-4, paste("Tail derivative", model))
    integ <- stats::integrate(function(z) exp(E$.cens_density(model, z, th, fixed)),
      grid[2], if (model == "beta") 1 else Inf, rel.tol = 1e-7)$value
    check(abs(integ - exp(s[2])) < 1e-6, paste("Tail integral", model))
  }
  nn <- if (model == "rician") 150L else 2000L
  lower <- rep(grid[2], nn)
  z <- E$.cens_truncated(model, lower, th, fixed)
  check(all(z > lower) && all(is.finite(z)), paste("Imputation bounds", model))
  if (!model %in% c("geometric", "Poisson", "negative binomial")) {
    uniforms <- exp(E$.cens_survival(model, z, th, fixed) - E$.cens_survival(model, lower, th, fixed))
    check(abs(mean(uniforms) - 0.5) < 4 / sqrt(12 * nn), paste("Truncation distribution", model))
  }
}

# API validation, paired missingness, all censored exceptions, and ties.
bad <- list(
  quote(fitcensBayes(1:3, c(1, 0), "gamma", "jeffreys")),
  quote(fitcensBayes(1:3, c(1, 2, 0), "gamma", "jeffreys")),
  quote(fitcensBayes(c(1, Inf), c(1, 0), "exponential", "reference")),
  quote(fitcensBayes(c(1, NA), c(1, 0), "exponential", "reference")),
  quote(fitcensBayes(1:3, c(0, 0, 0), "gamma", "jeffreys")),
  quote(fitcensBayes(1:3, c(1, 0, 0), "gamma", "jeffreys")),
  quote(fitcensBayes(c(1, 1, 3), c(1, 1, 0), "gamma", "jeffreys")),
  quote(fitcensBayes(1:4, c(1, 1, 0, 0), "gamma", "mdi")),
  quote(fitcensBayes(1:4, c(1, 1, 0, 0), "weibull", "mdi")),
  quote(fitcensBayes(1:4, c(1, 1, 0, 0), "lognormal", "mdi")),
  quote(fitcensBayes(c(0.5, 2), c(1, 0), "Poisson", "jeffreys")),
  quote(fitcensBayes(c(0.2, 1), c(1, 0), "beta", "jeffreys")),
  quote(fitcensBayes(1:4, c(1, 1, 0, 0), "negative binomial", "jeffreys")),
  quote(fitcensBayes(c(1, 1, 2, 4), c(1, 1, 1, 0), "t", "independence-jeffreys")),
  quote(fitcensBayes(1:3, c(1, 1, 0), "normal", "reference", start = c(sd = -1))),
  quote(fitcensBayes(1:3, c(1, 1, 0), "normal", "reference", chains = 1)),
  quote(fitcensBayes(c(1e308, 1e308), c(1, 0), "geometric", "reference")),
  quote(fitcensBayes(1:3, c(1, 1, 0), "normal", "reference", control = list(typo = 1)))
)
for (expr in bad) { expect_error(eval(expr)); counter$checks <- counter$checks + 1L }
fit <- fitcensBayes(c(1, NA, 3, 4), c(1, 0, 0, 1), "exponential", "reference",
                   na.action = "omit", iter = 200, chains = 2, seed = 1)
check(identical(fit$omitted, 2L) && identical(fit$status, c(1L, 0L, 1L)), "Paired NA omission.")
fit0 <- fitcensBayes(c(1, 2), c(0, 0), "exponential", "mdi", iter = 200, chains = 2, seed = 1)
check(fit0$censoring$events == 0 && fit0$diagnostics$converged, "All censored MDI exponential.")
fitz <- fitcensBayes(c(0, 2), c(1, 0), "exponential", "reference", iter = 200, chains = 2, seed = 1)
check(fitz$censoring$events == 1, "Zero event plus positive censored exposure.")
fit2 <- fitcensBayes(c(1, NA, 3, 4), c(1, 0, 0, 1), "exponential", "reference",
                    na.action = "omit", iter = 200, chains = 2, seed = 1)
check(identical(fit$chains, fit2$chains), "Reproducibility.")
ll <- log_lik_cens(fit, draws = 10, seed = 2)
check(identical(dim(ll), c(10L, 3L)) && all(is.finite(ll)), "Loglik shape.")
imp <- predict(fit, type = "impute", draws = 30, seed = 2)
check(all(imp > 3), "Predictive imputation.")
surv <- predict(fit, type = "survival", times = c(0, 1, 2), draws = 20, seed = 2)
check(all(surv[, 1] == 1) && all(surv[, 2] >= surv[, 3]), "Survival prediction.")
check(identical(dim(confint(fit)), c(1L, 2L)), "Credible intervals.")
check(inherits(summary(fit), "summary.fitcensBayes"), "Summary class.")

# Raw pointwise values must use the OBSERVED likelihood, not imputed data.
theta <- c(rate = 0.7)
check(max(abs(fit$.loglik(theta) - c(log(0.7) - 0.7, -2.1, log(0.7) - 2.8))) < 1e-12,
      "Observed likelihood mismatch.")

# Transforms and Jacobians, including gamma/mean, beta/precision, WL orthogonality.
for (model in names(fixtures)) {
  th <- fixtures[[model]]
  trans <- E$.cens_transform(model, th); u <- trans$to(th)
  check(max(abs(trans$from(u)$theta - th)) < 1e-9, paste("Round-trip", model))
  h <- 1e-5
  J <- sapply(seq_along(u), function(j) {
    a <- b <- u; a[j] <- a[j] - h; b[j] <- b[j] + h
    (trans$from(b)$theta - trans$from(a)$theta) / (2 * h)
  })
  lj <- as.numeric(determinant(matrix(J, length(u)), logarithm = TRUE)$modulus)
  check(abs(lj - trans$from(u)$jacobian) < 1e-6, paste("Jacobian", model))
}

# A frozen chain must never receive a reassuring convergence flag.
frozen <- replicate(4, matrix(1, 100, 1, dimnames = list(NULL, "p")), simplify = FALSE)
diag <- E$.cens_summarize(frozen, FALSE)
check(is.infinite(diag$rhat) && diag$ess_bulk == 0, "Frozen chains falsely passed.")
set.seed(515)
raw <- replicate(4, matrix(rnorm(500), 500, 1, dimnames = list(NULL, "p")), simplify = FALSE)
small <- lapply(raw, function(z) z * 1e-9)
large <- lapply(raw, function(z) z * 1e9)
check(abs(E$.cens_summarize(small,FALSE)$ess_mean - E$.cens_summarize(large,FALSE)$ess_mean) < 1e-8,
      "Mean ESS changed with measurement units.")
arguments <- c(500, 999, 1000, 1001, 2000, 10000)
check(max(abs(exp(E$.cens_log_i0_scaled(log(arguments))) /
  besselI(arguments, 0, expon.scaled = TRUE) - 1)) < 1e-12,
  "Scaled Bessel asymptotic continuity/accuracy.")
for (a in c(10, 400, 1000)) {
  th <- c(noncentrality = a, scale = 1)
  ss <- E$.cens_survival("rician", a+c(-3,0,3,10), th, NULL)
  check(all(is.finite(ss)) && all(diff(ss) < 0), "Large-ncp Rician survival.")
}
cat(sprintf("PASS: %d deterministic/numerical contract assertions; %d transformed original-prior comparisons.\n", counter$checks, prior_checks))
