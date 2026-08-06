library("fitdistrBayes")

ctrl <- list(
  rhat_threshold = 10, ess_threshold = 1,
  warn_convergence = FALSE, proposal_scale = 0.3
)

generate <- list(
  gumbel = function(n, p) p[1] - p[2] * log(-log(runif(n))),
  frechet = function(n, p) (p[2] / rexp(n))^(1 / p[1]),
  lomax = function(n, p) p[2] * expm1(-log(runif(n)) / p[1]),
  nakagami = function(n, p) sqrt(rgamma(n, p[1], rate = p[1] / p[2])),
  el = function(n, p) {
    u <- runif(n)
    -(log(-expm1((1 - u) * log(p[1]))) - log1p(-p[1])) / p[2]
  },
  rician = function(n, p) {
    sqrt(rnorm(n, p[1], p[2])^2 + rnorm(n, 0, p[2])^2)
  }
)

parameters <- list(
  gumbel = list(c(-10, 0.05), c(0, 1), c(1000, 20)),
  frechet = list(c(0.7, 0.1), c(2.5, 4), c(8, 200)),
  lomax = list(c(0.7, 0.1), c(3, 3), c(10, 100)),
  nakagami = list(c(0.2, 0.1), c(0.5, 4), c(5, 100)),
  el = list(c(0.01, 0.1), c(0.4, 1.3), c(0.95, 10)),
  rician = list(c(0.1, 2), c(5, 2), c(50, 0.5))
)

priors <- list(
  gumbel = c("jeffreys", "reference", "mdi"),
  frechet = c("jeffreys", "reference"),
  lomax = "jeffreys",
  nakagami = c("jeffreys", "reference"),
  el = c("jeffreys", "mdi", "reference-theta", "reference-rate"),
  rician = "jeffreys"
)

sizes <- list(
  gumbel = c(2L, 8L, 60L),
  frechet = c(2L, 8L, 60L),
  lomax = c(1L, 8L, 60L),
  nakagami = c(2L, 8L, 60L),
  el = c(1L, 8L, 60L),
  rician = c(3L, 8L, 60L)
)

public_name <- c(
  gumbel = "gumbel", frechet = "frechet", lomax = "lomax",
  nakagami = "nakagami-m", el = "exponential-logarithmic",
  rician = "rician"
)

set.seed(20260720)
records <- list()
position <- 0L
for (model in names(generate)) {
  for (scenario in seq_along(parameters[[model]])) {
    for (n in sizes[[model]]) {
      x <- generate[[model]](n, parameters[[model]][[scenario]])
      # A continuous simulation can be exactly constant only through floating
      # point collapse at an extreme parameter setting.
      if (model %in% c("gumbel", "frechet", "rician") &&
          length(unique(x)) < 2L) next
      for (prior in priors[[model]]) {
        if (model == "el" && n <= 2L && prior != "mdi") next
        position <- position + 1L
        fit_result <- tryCatch(
          list(
            value = fitdistrBayes(
              x, public_name[[model]], prior,
              iter = 240, warmup = 120, chains = 2,
              seed = 5000 + position, control = ctrl
            ),
            error = NULL
          ),
          error = function(e) {
            list(value = NULL, error = conditionMessage(e))
          }
        )
        fit <- fit_result$value
        fit_error <- fit_result$error
        if (is.null(fit)) {
          if (!grepl("no clipping|could not be represented|floating-point range",
                     fit_error, ignore.case = TRUE)) {
            stop(sprintf(
              "Stress failure for %s / %s / scenario %d / n=%d: %s",
              model, prior, scenario, n, fit_error
            ), call. = FALSE)
          }
          records[[position]] <- data.frame(
            model = model, prior = prior, scenario = scenario, n = n,
            max_rhat = NA_real_, min_ess = NA_real_,
            status = "explicit numerical-range stop"
          )
          next
        }
        draws <- as.matrix(fit$draws[, fit$model$parameters, drop = FALSE])
        predictive <- predict(fit, draws = 12, seed = 9000 + position)
        likelihood <- log_lik(fit, draws = 6)
        positive_parameters <- setdiff(colnames(draws), "location")
        bad_positive <- length(positive_parameters) &&
          any(draws[, positive_parameters, drop = FALSE] <= 0)
        if (any(!is.finite(draws)) || bad_positive ||
            any(!is.finite(predictive)) || any(!is.finite(likelihood))) {
          stop(sprintf(
            "Stress failure for %s / %s / scenario %d / n=%d.",
            model, prior, scenario, n
          ), call. = FALSE)
        }
        if (model == "el" && any(draws[, "theta"] >= 1)) {
          stop("EL theta draw reached the closed boundary.", call. = FALSE)
        }
        records[[position]] <- data.frame(
          model = model, prior = prior, scenario = scenario, n = n,
          max_rhat = fit$diagnostics$max_rhat,
          min_ess = fit$diagnostics$min_ess_bulk,
          status = "fitted"
        )
      }
    }
  }
}

records <- do.call(rbind, records)
cat(sprintf(
  "New-model stress validation completed: %d fits across wide parameter and sample-size regimes.\n",
  nrow(records)
))
print(aggregate(cbind(max_rhat, min_ess) ~ model, records, max),
      row.names = FALSE)
print(table(records$status))
