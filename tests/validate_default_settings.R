library("fitdistrBayes")

set.seed(19072026)
datasets <- list(
  beta = rbeta(50, 2.5, 6),
  cauchy = rcauchy(50, 1, 1.5),
  gamma = rgamma(50, 3, rate = 2),
  logistic = rlogis(50, -0.5, 1.2),
  t_fixed = 1 + 1.5 * rt(50, 7),
  t_unknown = rt(60, 9),
  weibull = rweibull(50, 1.8, 2.2),
  cauchy_small = rcauchy(50, 0, 0.01),
  logistic_offset = 1e8 + rlogis(50, 0, 1e-3),
  t_offset = 1e8 + 1e-3 * rt(50, 30)
)

fits <- list(
  beta = fitdistrBayes(datasets$beta, "beta", "Jeffreys", seed = 11),
  cauchy = fitdistrBayes(datasets$cauchy, "cauchy", "reference", seed = 12),
  gamma = fitdistrBayes(datasets$gamma, "gamma", "reference-shape", seed = 13),
  logistic = fitdistrBayes(datasets$logistic, "logistic", "reference", seed = 14),
  t_fixed = fitdistrBayes(datasets$t_fixed, "t", "reference",
                          fixed = list(df = 7), seed = 15),
  t_unknown = fitdistrBayes(datasets$t_unknown, "t",
                            "independence-jeffreys", seed = 16),
  weibull = fitdistrBayes(datasets$weibull, "weibull", "reference", seed = 17),
  cauchy_small = fitdistrBayes(datasets$cauchy_small, "cauchy",
                               "Jeffreys", seed = 18),
  logistic_offset = fitdistrBayes(datasets$logistic_offset, "logistic",
                                  "MDI", seed = 19),
  t_offset = fitdistrBayes(datasets$t_offset, "t", "reference",
                           fixed = list(df = 30), seed = 20)
)

result <- do.call(rbind, lapply(names(fits), function(name) {
  fit <- fits[[name]]
  data.frame(
    model = name,
    algorithm = fit$engine$algorithm,
    converged = fit$diagnostics$converged,
    max_rhat = fit$diagnostics$max_rhat,
    min_ess_bulk = fit$diagnostics$min_ess_bulk,
    min_ess_tail = fit$diagnostics$min_ess_tail,
    mean_acceptance = if (is.null(fit$diagnostics$acceptance)) NA_real_
      else mean(fit$diagnostics$acceptance$overall),
    row.names = NULL
  )
}))

print(result, row.names = FALSE, digits = 4)
if (!all(result$converged)) {
  stop("At least one default-setting validation fit missed its convergence targets.",
       call. = FALSE)
}
cat("All default-setting validation fits satisfied R-hat and ESS targets.\n")
