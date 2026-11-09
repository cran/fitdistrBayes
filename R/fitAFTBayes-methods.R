print.fitAFTBayes <- function(x, digits = 4L, ...) {
  cat("\nBayesian ", x$model$name, " regression\n", sep = "")
  cat("Prior:", x$prior$key, "| Observations:", x$model$n, "| Failures:", x$censoring$events, "\n")
  cat("Diagnostics:", if (x$diagnostics$converged) "targets met" else "attention required", "\n\n")
  print(x$summary[, c("parameter", "median", "q2.5", "q97.5", "rhat", "ess_bulk", "ess_tail")],
    row.names = FALSE, digits = digits)
  cat("\ncoef() returns regression medians; use type='time_ratio' for exp(beta).\n")
  if (!is.null(x$quantile_precision)) {
    ratio <- x$diagnostics$max_quantile_mcse_ratio
    cat("Maximum approximate endpoint MCSE / 95% interval width:",
        if (is.finite(ratio)) format(ratio, digits = digits) else "not reliably estimated", "\n")
    cat("Inspect $quantile_precision for interval accuracy; R-hat/ESS targets alone do not set precision.\n")
  }
  invisible(x)
}
summary.fitAFTBayes <- function(object, ...) {
  ans <- object[c("call", "prior", "model", "summary", "moment_status", "quantile_precision", "diagnostics", "censoring", "propriety")]
  class(ans) <- "summary.fitAFTBayes"; ans
}
print.summary.fitAFTBayes <- function(x, digits = 4L, ...) {
  print.fitAFTBayes(x, digits = digits)
  cat("\nPosterior moments (NA means not certified):\n")
  print(x$moment_status, row.names = FALSE)
  if (!is.null(x$quantile_precision)) {
    cat("\nApproximate Monte Carlo standard errors for posterior quantiles:\n")
    print(x$quantile_precision, row.names = FALSE, digits = digits)
  }
  cat("\n", x$diagnostics$messages, "\n", sep = "")
  invisible(x)
}
coef.fitAFTBayes <- function(object, type = c("regression", "all", "time_ratio"), ...) {
  type <- match.arg(type)
  if (type == "all") return(object$estimates)
  ans <- stats::setNames(object$estimates[seq_len(object$model$p)], colnames(object$design))
  if (type == "time_ratio") {
    ans <- exp(ans[names(ans) != "(Intercept)"])
    if (any(!is.finite(ans) | ans == 0)) stop("Time ratio is outside the representable range; use regression coefficients.")
  }
  ans
}
confint.fitAFTBayes <- function(object, parm = NULL, level = 0.95, type = c("regression", "all", "time_ratio"), ...) {
  type <- match.arg(type)
  if (!is.numeric(level) || is.complex(level) || length(level)!=1L ||
      !is.finite(level) || level<=0 || level>=1) stop("'level' must lie in (0,1).")
  available <- if (type == "all") object$model$parameters else colnames(object$design)
  if (type == "time_ratio") available <- setdiff(available, "(Intercept)")
  if (is.null(parm)) parm <- available
  if (is.numeric(parm)) {
    if (any(!is.finite(parm)) || any(parm != floor(parm)) || any(parm < 1 | parm > length(available))) stop("Invalid parameter indices.")
    parm <- available[parm]
  }
  if (!is.character(parm) || anyNA(parm) || any(!parm %in% available)) stop("Unknown AFT parameter.")
  if (!length(parm)) return(matrix(numeric(), 0L, 2L, dimnames = list(character(), c("lower", "upper"))))
  internal <- if (type == "all") parm else paste0("beta:", parm)
  ans <- confint.fitdistrBayes(object, parm = internal, level = level)
  rownames(ans) <- parm
  if (type == "time_ratio") {
    ans <- exp(ans)
    if (any(!is.finite(ans) | ans == 0)) stop("Time-ratio interval is outside the representable range; use regression intervals.")
  }
  ans
}
as.data.frame.fitAFTBayes <- function(x, row.names = NULL, optional = FALSE, ...) x$draws
log_lik.fitAFTBayes <- function(object, draws = NULL, seed = NULL, ...) {
  if (length(list(...))) stop("Unused arguments in log_lik.fitAFTBayes().")
  log_lik.fitdistrBayes(object, draws = draws, seed = seed)
}
nobs.fitAFTBayes <- function(object, ...) object$model$n
formula.fitAFTBayes <- function(x, ...) x$formula
model.matrix.fitAFTBayes <- function(object, ...) object$design
plot.fitAFTBayes <- function(x, type = c("trace", "density", "acf", "pairs"), pars = x$model$parameters, ...) {
  .cens_plot(x, match.arg(type), pars, ...)
}

.aft_newdesign <- function(object, newdata) {
  if (is.null(newdata)) return(object$design)
  tt <- stats::delete.response(object$terms)
  mf <- stats::model.frame(tt, data = newdata, na.action = stats::na.fail, xlev = object$xlevels)
  X <- stats::model.matrix(tt, mf, contrasts.arg = object$contrasts)
  if (!identical(colnames(X), colnames(object$design)) || !nrow(X) || any(!is.finite(X)))
    stop("Prediction design must be finite and have the fitted design columns.")
  X
}

predict.fitAFTBayes <- function(object, newdata = NULL, type = c("survival", "linear"),
                                times = NULL, summary = TRUE, level = 0.95,
                                draws = NULL, seed = NULL, ...) {
  if (length(list(...))) stop("Unused arguments in predict.fitAFTBayes().")
  type <- match.arg(type)
  if (!is.logical(summary) || length(summary) != 1L || is.na(summary)) stop("'summary' must be TRUE or FALSE.")
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) stop("'level' must lie in (0,1).")
  X <- .aft_newdesign(object, newdata); p <- object$model$p; n <- nrow(X)
  id <- seq_len(nrow(object$draws))
  if (!is.null(draws)) {
    draws <- .fdb_scalar_count(draws, "draws")
    id <- .fdb_with_seed(seed, sample.int(length(id), draws, replace = draws > length(id)))
  }
  dr <- as.matrix(object$draws[id, object$model$parameters, drop = FALSE]); S <- nrow(dr)
  eta <- dr[, seq_len(p), drop = FALSE] %*% t(X)
  if (any(!is.finite(eta))) stop("Linear prediction exceeds the representable range.")
  alpha <- (1-level)/2
  summarize <- function(a, index) {
    q <- stats::quantile(a, c(alpha, .5, 1-alpha), names = FALSE, type = 8)
    data.frame(index, estimate = if (type == "survival") mean(a) else q[2L], median = q[2L], lower = q[1L], upper = q[3L], row.names = NULL)
  }
  if (type == "linear") {
    if (!is.null(times)) stop("'times' is used only for survival prediction.")
    if (!summary) return(eta)
    return(do.call(rbind, lapply(seq_len(n), function(i) summarize(eta[, i], data.frame(row = i)))))
  }
  if (!is.numeric(times) || is.complex(times) || !is.null(dim(times)) || !length(times) ||
      anyNA(times) || any(times < 0)) stop("Supply nonnegative prediction times (including 0 or Inf if desired).")
  values <- array(NA_real_, c(S, n, length(times)), dimnames = list(NULL, rownames(X), as.character(times)))
  family <- object$model$aft_family
  for (k in seq_len(S)) for (j in seq_along(times)) {
    if (is.null(family)) {
      values[k, , j] <- exp(.aft_logsurvival((log(times[j])-eta[k, ])/dr[k, p+1L], dr[k, p+2L]))
    } else {
      m <- family
      if (m$shape_estimated) {m$shape <- dr[k,p+1L]; scale <- 1} else
        scale <- if (is.null(m$fixed_scale)) dr[k,p+1L] else m$fixed_scale
      values[k,,j] <- exp(.aft_error_logsurvival((log(times[j])-eta[k,])/scale,m))
    }
  }
  if (any(!is.finite(values)) || any(values < 0 | values > 1)) stop("Nonfinite or invalid survival prediction; no clipping was applied.")
  if (!summary) return(values)
  do.call(rbind, lapply(seq_len(n), function(i) do.call(rbind, lapply(seq_along(times),
    function(j) summarize(values[, i, j], data.frame(row = i, time = times[j]))))))
}
