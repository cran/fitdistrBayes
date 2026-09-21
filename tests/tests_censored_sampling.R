source("censored-helpers.R", local = TRUE)

# Short deterministic smoke runs check API/numerics, not MCMC convergence.
routes <- fitcensBayes_models()
for (r in seq_len(nrow(routes))) {
  model <- E$.fdb_model_name(routes$model[r])
  prior <- routes$prior[r]
  y <- simulate_model(model, n = 60, seed = 800 + r)
  d <- censor_data(y, 0.35)
  fit <- fitcensBayes(d$x, d$status, model, prior,
    fixed = get_fixed(model, prior), iter = 80, warmup = 40, chains = 2,
    seed = 900 + r, control = control_quiet)
  stopifnot(inherits(fit, "fitcensBayes"),
    !inherits(fit, "fitdistrBayes"),
    nrow(as.data.frame(fit)) == 80L,
    identical(dim(log_lik(fit, draws = 3)), c(3L, 60L)),
    all(is.finite(log_lik_cens(fit, draws = 3))))
  imp <- predict(fit, type = "impute", draws = 2, seed = 42)
  stopifnot(all(sweep(imp, 2, d$x[d$status == 0], "-") > 0))
}
cat("PASS: all 56 censored routes and registered extraction/prediction methods.\n")

# Each model also exercises its augmentation route; no convergence claim.
for (model in names(fixtures)) {
  prior <- E$fitdistrBayes_routes(model)$prior[1]
  y <- simulate_model(model, n = 25, seed = 80)
  d <- censor_data(y, 0.2)
  fit <- fitcensBayes(d$x, d$status, model, prior,
    fixed = get_fixed(model, prior), method = "augmentation",
    iter = 40, warmup = 20, chains = 2, seed = 90,
    control = control_quiet)
  stopifnot(all(is.finite(as.matrix(fit$draws[, fit$model$parameters, drop = FALSE]))))
}
cat("PASS: all 20 augmentation model routes.\n")

# Independent exact-posterior checks.
x <- c(0.3, 0.7, 1.2, 1.5, 2)
status <- c(1, 1, 0, 1, 0)
for (prior in c("jeffreys", "reference", "mdi")) {
  fit <- fitcensBayes(x, status, "exponential", prior,
    iter = 6000, warmup = 1000, chains = 2, seed = 700)
  a <- sum(status) + if (prior == "mdi") 2 else 0
  b <- sum(x)
  stopifnot(abs(mean(fit$draws$rate) - a / b) <
    5 * sqrt(a / (b^2 * nrow(fit$draws))))
}

# Serialization must preserve callbacks and public S3 dispatch.
path <- tempfile(fileext = ".rds")
saveRDS(fit, path)
restored <- readRDS(path)
unlink(path)
stopifnot(identical(log_lik(fit), log_lik(restored)),
  identical(predict(fit, draws = 3, seed = 88),
            predict(restored, draws = 3, seed = 88)),
  identical(coef(fit), coef(restored)),
  identical(confint(fit), confint(restored)))
cat("PASS: exact posteriors and serialized fitted object.\n")
