# Public methods for the independent fitcensBayes class.
fitcensBayes_models <- function() .cens_models()
print.fitcensBayes <- function(x, digits = 4L, ...) .cens_print(x, digits, ...)
summary.fitcensBayes <- function(object, ...) {
  ans <- object[c("model", "prior", "summary", "moment_status", "diagnostics",
                  "initialization", "engine", "censoring")]
  class(ans) <- "summary.fitcensBayes"
  ans
}
print.summary.fitcensBayes <- function(x, digits = 4L, ...) {
  .cens_print(x, digits, ...)
  cat("\nDiagnostics:\n", paste(x$diagnostics$messages, collapse = "\n"), "\n")
  cat("\nPosterior moment assessment:\n")
  print(x$moment_status, row.names = FALSE)
  invisible(x)
}
coef.fitcensBayes <- function(object, ...) object$estimates
confint.fitcensBayes <- function(object, parm = object$model$parameters,
                               level = 0.95, ...) {
  confint.fitdistrBayes(object, parm, level, ...)
}
as.data.frame.fitcensBayes <- function(x, row.names = NULL, optional = FALSE, ...) x$draws
plot.fitcensBayes <- function(x, type = c("trace", "density", "acf", "pairs"),
                             pars = x$model$parameters, ...) {
  .cens_plot(x, match.arg(type), pars, ...)
}
predict.fitcensBayes <- function(object, type = c("response", "survival", "impute"),
                                times = NULL, draws = 1000L, size = 1L,
                                seed = NULL, ...) {
  if (length(list(...))) stop("Unused arguments in predict.fitcensBayes().")
  .cens_predict(object, match.arg(type), times, draws, size, seed)
}
log_lik.fitcensBayes <- function(object, draws = NULL, seed = NULL, ...) {
  if (length(list(...))) stop("Unused arguments in log_lik.fitcensBayes().")
  log_lik.fitdistrBayes(object, draws, seed)
}
log_lik_cens <- function(object, draws = NULL, seed = NULL) {
  if (!inherits(object, "fitcensBayes")) stop("Expected a fitcensBayes object.")
  log_lik.fitcensBayes(object, draws, seed)
}
