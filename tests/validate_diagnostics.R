library("fitdistrBayes")

if (!requireNamespace("posterior", quietly = TRUE)) {
  cat("Package 'posterior' is unavailable; independent diagnostic comparison skipped.\n")
  quit(save = "no", status = 0)
}

set.seed(1926)
n <- 1500
m <- 4
rho <- 0.8
chains <- lapply(seq_len(m), function(ch) {
  z <- numeric(n)
  z[1] <- rnorm(1, sd = 2)
  for (i in 2:n) z[i] <- rho * z[i - 1] + rnorm(1)
  matrix(z + c(-0.05, 0, 0.03, 0.01)[ch],
         ncol = 1, dimnames = list(NULL, "theta"))
})

ours <- fitdistrBayes:::.fdb_summarize(chains, independent = FALSE)
array_draws <- array(
  unlist(lapply(chains, function(z) z[, 1])),
  dim = c(n, m, 1),
  dimnames = list(iteration = NULL, chain = NULL, variable = "theta")
)
reference <- posterior::summarise_draws(
  posterior::as_draws_array(array_draws),
  rhat = posterior::rhat,
  ess_mean = posterior::ess_mean,
  ess_bulk = posterior::ess_bulk,
  ess_tail = posterior::ess_tail,
  mcse_mean = posterior::mcse_mean
)

if (abs(ours$rhat - reference$rhat) > 0.01) {
  stop(sprintf("R-hat mismatch: prototype=%g, posterior=%g.",
               ours$rhat, reference$rhat), call. = FALSE)
}
relative_bulk <- abs(ours$ess_bulk - reference$ess_bulk) /
  reference$ess_bulk
relative_tail <- abs(ours$ess_tail - reference$ess_tail) /
  reference$ess_tail
if (relative_bulk > 0.20 || relative_tail > 0.20) {
  stop(sprintf(
    "ESS mismatch: prototype bulk/tail=%g/%g; posterior=%g/%g.",
    ours$ess_bulk, ours$ess_tail,
    reference$ess_bulk, reference$ess_tail
  ), call. = FALSE)
}
relative_mean <- abs(ours$ess_mean - reference$ess_mean) /
  reference$ess_mean
relative_mcse <- abs(ours$mcse_mean - reference$mcse_mean) /
  reference$mcse_mean
if (relative_mean > 0.20 || relative_mcse > 0.20) {
  stop(sprintf(
    "Mean ESS/MCSE mismatch: prototype=%g/%g; posterior=%g/%g.",
    ours$ess_mean, ours$mcse_mean,
    reference$ess_mean, reference$mcse_mean
  ), call. = FALSE)
}

cat(sprintf(
  paste0("Diagnostic validation passed. R-hat prototype/reference: %.4f/%.4f; ",
         "mean ESS: %.1f/%.1f; bulk ESS: %.1f/%.1f; ",
         "tail ESS: %.1f/%.1f; MCSE(mean): %.4f/%.4f.\n"),
  ours$rhat, reference$rhat,
  ours$ess_mean, reference$ess_mean,
  ours$ess_bulk, reference$ess_bulk,
  ours$ess_tail, reference$ess_tail,
  ours$mcse_mean, reference$mcse_mean
))
