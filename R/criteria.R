# Optional post-processing. No function in this file is called during sampling.

.fdb_criteria_methods <- function(methods) {
  if (identical(methods, TRUE)) return(c("waic", "looic", "dic"))
  if (!is.character(methods) || !length(methods) || anyNA(methods))
    stop("Expected TRUE or criterion names: 'waic', 'looic', 'dic'. Disable fitting-time criteria with criteria = FALSE.", call. = FALSE)
  methods <- tolower(methods)
  methods[methods == "loo"] <- "looic"
  if (any(!methods %in% c("waic", "looic", "dic")))
    stop("Unknown criterion; use 'waic', 'looic' (or 'loo'), or 'dic'.", call. = FALSE)
  unique(methods)
}

.fdb_criteria_preflight <- function(methods, control) {
  methods <- .fdb_criteria_methods(methods)
  if (!isTRUE(control$store_callables))
    stop("Criteria require control$store_callables = TRUE.", call. = FALSE)
  if ("looic" %in% methods && !requireNamespace("loo", quietly = TRUE))
    stop("PSIS-LOO requires the optional 'loo' package. Install it with install.packages('loo'), or request 'waic'/'dic' only.", call. = FALSE)
  invisible(methods)
}

# Linear-time sufficient propriety checks for all leave-one-out training sets.
# These mirror the complete-event guards in .fdb_build_builtin(), but do not
# initialize samplers or rebuild models for each observation.
.fdb_loo_guard <- function(object) {
  x <- object$data
  event <- if (inherits(object, "fitcensBayes")) object$status == 1L else rep(TRUE, length(x))
  model <- object$model$name
  prior <- object$prior$key
  if (identical(model, "dhillon") && !identical(prior, "user-defined")) {
    ok <- .fdb_dhillon_loo(x, event)
    return(list(certified = ok, reason = if (all(ok))
      "Every training posterior satisfies the exact Dhillon observed-data propriety condition." else
      "At least one Dhillon training posterior is improper; PSIS-LOO is withheld."))
  }
  if (identical(prior, "user-defined")) return(list(
    certified = rep(FALSE, length(x)),
    reason = "Leave-one-out posterior propriety is not certified automatically for models or priors defined by the user."))
  if (model == "exponential") {
    # Nonnegative times: counting positive terms avoids cancellation in sum(x)-x.
    ok <- length(x) > 1L & (sum(x > 0) - (x > 0) > 0) &
      (sum(event) - event + if (prior == "mdi") 2 else 0) > 0
  } else {
    z <- x[event]
    if (model == "lognormal") z <- log(z)
    groups <- match(z, unique(z))
    freq <- tabulate(groups)
    ni <- sum(event) - event
    distinct <- rep(length(freq), length(x))
    maxcount <- rep(if (length(freq)) max(freq) else 0L, length(x))
    if (length(z)) {
      distinct[event] <- length(freq) - (freq[groups] == 1L)
      largest <- max(freq)
      maxcount[event] <- largest - (sum(freq == largest) == 1L & freq[groups] == largest)
    }
    ok <- switch(model,
      "chi-squared" = ni >= 1L,
      "geometric" = ni >= 1L,
      "negative binomial" = ni >= 1L,
      "Poisson" = ni >= 1L,
      "lomax" = ni >= 1L,
      "exponential-logarithmic" = ni >= if (prior == "mdi") 1L else 3L,
      "rician" = ni > 2L & distinct > 1L,
      "cauchy" = 2 * maxcount < ni + if (prior == "jeffreys") 0 else 1,
      "t" = {
        if (is.null(object$model$fixed$df)) ni > 1L & distinct == ni else {
          a <- if (prior == "jeffreys") 2 else 1
          ni + a > 2 & object$model$fixed$df * (ni - maxcount) - maxcount - a + 2 > 0
        }
      },
      "beta" =, "frechet" =, "gumbel" =, "gamma" =, "lognormal" =,
      "logistic" =, "nakagami" =, "normal" =, "weibull" =,
      "weighted lindley" = ni >= 2L & distinct > 1L,
      rep(FALSE, length(x)))
    ok <- ok & ni >= 1L
  }
  list(certified = ok, reason = if (all(ok))
    "Every training posterior satisfies the registered sufficient condition (exact-event subset under censoring)." else
    "At least one training posterior fails the sufficient propriety check. PSIS-LOO is withheld; this is not by itself a proof of impropriety under censoring.")
}

.fdb_criteria_identity <- function(object) {
  discrete <- c("geometric", "negative binomial", "Poisson")
  measure <- if (identical(object$prior$key, "user-defined"))
    paste0("custom:", object$model$name) else
    if (object$model$name %in% discrete) "counting" else "Lebesgue"
  list(data = as.numeric(object$data),
       status = if (inherits(object, "fitcensBayes")) as.integer(object$status) else rep(1L, length(object$data)),
       omitted = as.integer(object$omitted), measure = measure)
}

.fdb_ic_sd <- function(x) {
  if (length(x) < 2L || any(!is.finite(x))) return(NA_real_)
  # stats::var() can overflow while summing a large constant column on R
  # builds without long double. Center first (also preserving small deviations
  # from a large offset), then bound the inputs to var() by one in magnitude.
  centered <- x - x[1L]
  if (all(is.finite(centered))) {
    scale <- max(abs(centered))
    if (scale == 0) return(0)
    standardized <- centered / scale
  } else {
    # Opposite-sign finite extremes can overflow during subtraction itself.
    scale <- max(abs(x))
    standardized <- x / scale
  }
  sqrt(stats::var(standardized)) * scale
}

.fdb_ic_variance <- function(x) .fdb_ic_sd(x)^2

.fdb_ic_se <- function(pointwise) {
  sqrt(length(pointwise)) * .fdb_ic_sd(pointwise)
}

# Censoring at the lower endpoint of a positive continuous distribution is
# the sure event T > 0. Its predictive contribution and Monte Carlo error are
# exactly zero; fitting a Pareto tail to its constant importance weights is
# neither necessary nor meaningful.
.fdb_zero_information <- function(object) {
  if (!inherits(object, "fitcensBayes") || !object$model$name %in%
      c("beta", "chi-squared", "exponential", "gamma", "lognormal", "weibull",
        "frechet", "lomax", "nakagami", "exponential-logarithmic", "rician", "weighted lindley", "dhillon"))
    return(integer())
  which(object$status == 0L & object$data == 0)
}

.fdb_loo_restore_zero <- function(result, zero, S, n) {
  if (!length(zero)) return(result)
  informative <- setdiff(seq_len(n), zero)
  pw <- matrix(0, n, ncol(result$pointwise), dimnames = list(NULL, colnames(result$pointwise)))
  pw[informative, ] <- result$pointwise
  # -Inf is a sentinel for an exact degenerate tail, not an estimated k.
  pw[zero, "influence_pareto_k"] <- -Inf
  result$pointwise <- pw
  for (nm in c("pareto_k", "n_eff", "r_eff")) {
    v <- rep(switch(nm, pareto_k = -Inf, n_eff = S, r_eff = 1), n)
    v[informative] <- result$diagnostics[[nm]]
    result$diagnostics[[nm]] <- v
  }
  for (nm in c("elpd_loo", "p_loo", "looic")) {
    result$estimates[nm, ] <- c(sum(pw[, nm]), .fdb_ic_se(pw[, nm]))
    result[[nm]] <- result$estimates[nm, "Estimate"]
    result[[paste0("se_", nm)]] <- result$estimates[nm, "SE"]
  }
  attr(result, "dims") <- c(S, n)
  attr(result, "zero_information_observations") <- zero
  result
}

.fdb_criteria_after_fit <- function(object, methods) {
  methods <- .fdb_criteria_methods(methods)
  tryCatch(.fdb_compute_criteria(object, methods), error = function(e) {
    reason <- paste0("Posterior draws were retained, but optional criteria could not be computed: ", conditionMessage(e))
    warning(reason, call. = FALSE)
    ans <- list(estimates = data.frame(criterion = methods, estimate = NA_real_,
      se = NA_real_, p_eff = NA_real_, elpd = NA_real_, available = FALSE,
      reliable = FALSE, reason = reason), details = list(),
      diagnostics = list(sampling = object$diagnostics, loo_propriety = NULL,
        warnings = reason, computation_error = conditionMessage(e)),
      model = object$model$name, prior = object$prior$key,
      nobs = length(object$data), ndraws = nrow(object$draws),
      observation = .fdb_criteria_identity(object))
    class(ans) <- "fitdistrBayes_criteria"
    ans
  })
}

.fdb_compute_criteria <- function(object, methods, cores = 1L) {
  if (!inherits(object, c("fitdistrBayes", "fitcensBayes")))
    stop("'object' must be a fitdistrBayes or fitcensBayes fit.", call. = FALSE)
  methods <- .fdb_criteria_methods(methods)
  .fdb_criteria_preflight(methods, object$control)
  cores <- .fdb_scalar_count(cores, "cores", lower = 1L)
  if (cores > 2L) stop("Use at most two cores; the default is one.", call. = FALSE)
  if (!is.function(object$.loglik)) stop("Pointwise log-likelihood callback was not stored.", call. = FALSE)
  notes <- character()
  issue <- function(msg) {
    notes <<- unique(c(notes, msg))
    warning(msg, call. = FALSE)
  }
  if (!isTRUE(object$diagnostics$converged))
    issue("Sampling diagnostic targets were not met. Information criteria may be unreliable; inspect the chains before comparing models.")
  if (identical(object$prior$key, "user-defined"))
    issue("Criteria for models or priors defined by the user rely on normalized pointwise densities and valid posterior/moment declarations supplied by the user.")
  params <- object$model$parameters
  means_ok <- object$moment_status$mean_exists[match(params, object$moment_status$parameter)]
  dic_ok <- length(means_ok) == length(params) && all(means_ok %in% TRUE)
  guard <- if ("looic" %in% methods) .fdb_loo_guard(object) else NULL
  need_ll <- "waic" %in% methods || ("looic" %in% methods && all(guard$certified)) ||
    ("dic" %in% methods && dic_ok)
  ll <- NULL
  if (need_ll) {
    ll <- log_lik(object)
    if (!is.matrix(ll) || !is.numeric(ll) || nrow(ll) < 2L ||
        nrow(ll) != nrow(object$draws) || ncol(ll) != length(object$data) || any(!is.finite(ll)))
      stop("Criteria require a finite draws-by-observations log-likelihood matrix (at least two draws). No clipping was applied.", call. = FALSE)
  }
  S <- nrow(object$draws)
  n <- length(object$data)
  details <- list()
  rows <- list()
  row <- function(name, value = NA_real_, se = NA_real_, p = NA_real_, elpd = NA_real_,
                  reliable = FALSE, reason = "") {
    if (any(is.infinite(c(value, p, elpd))) || any(is.nan(c(value, p, elpd)))) {
      reason <- paste0(toupper(name), " aggregate exceeded the numerical range; no numeric score was returned.")
      issue(reason)
      value <- se <- p <- elpd <- NA_real_
    }
    if (is.finite(value) && name != "dic" && n > 1L && !is.finite(se)) {
      msg <- paste0(toupper(name), " standard error could not be represented reliably.")
      issue(msg); reason <- paste(reason, msg); se <- NA_real_; reliable <- FALSE
    }
    if (identical(object$prior$key, "user-defined") &&
        !identical(object$model$extension$propriety_source, "user-supplied")) {
      reliable <- FALSE
      reason <- paste(reason, "Posterior propriety was not declared for this model.")
    }
    data.frame(criterion = name, estimate = value, se = se, p_eff = p, elpd = elpd,
               available = is.finite(value), reliable = isTRUE(reliable) && is.finite(value), reason = trimws(reason),
               stringsAsFactors = FALSE)
  }
  if ("waic" %in% methods) {
    maxima <- apply(ll, 2L, max)
    lppd <- maxima + log(colMeans(exp(sweep(ll, 2L, maxima, "-"))))
    penalty <- apply(ll, 2L, .fdb_ic_variance)
    pw <- cbind(elpd_waic = lppd - penalty, p_waic = penalty, waic = -2 * (lppd - penalty))
    bad <- which(penalty > 0.4)
    outside <- which(rowSums(!is.finite(pw)) > 0L)
    if (length(outside)) {
      reason <- "WAIC pointwise arithmetic exceeded the numerical range; no numeric score was returned."
      issue(reason)
      pw[!is.finite(pw)] <- NA_real_
      details$waic <- list(pointwise = pw, high_p_waic = bad,
                          unrepresentable_observations = outside)
      rows$waic <- row("waic", reason = reason)
    } else {
      if (length(bad)) issue("Some pointwise p_waic values exceed 0.4; WAIC may be unreliable. Inspect PSIS-LOO diagnostics.")
      details$waic <- list(pointwise = pw, high_p_waic = bad)
      rows$waic <- row("waic", sum(pw[, "waic"]), .fdb_ic_se(pw[, "waic"]), sum(penalty),
        sum(pw[, "elpd_waic"]), isTRUE(object$diagnostics$converged) && !length(bad),
        if (length(bad)) "Pointwise p_waic > 0.4." else "")
    }
  }
  if ("looic" %in% methods) {
    if (!all(guard$certified)) {
      issue(guard$reason)
      rows$looic <- row("looic", reason = guard$reason)
      details$looic <- NULL
    } else {
      zero <- .fdb_zero_information(object)
      informative <- setdiff(seq_len(n), zero)
      ll_loo <- ll[, informative, drop = FALSE]
      if (isTRUE(object$diagnostics$exact_or_independent)) {
        r_eff <- rep(1, length(informative))
      } else {
        likelihood <- exp(sweep(ll_loo, 2L, apply(ll_loo, 2L, max), "-"))
        r_eff <- loo::relative_eff(likelihood, chain_id = object$draws$.chain, cores = cores)
        constant <- apply(ll_loo, 2L, function(z) all(z == z[1L]))
        r_eff[constant] <- 1
        if (length(r_eff) != length(informative) || any(!is.finite(r_eff) | r_eff <= 0))
          stop("Relative effective sample size could not be estimated. Inspect chains and likelihood values.", call. = FALSE)
      }
      result <- withCallingHandlers(loo::loo(ll_loo, r_eff = r_eff, cores = cores),
        warning = function(w) { notes <<- unique(c(notes, conditionMessage(w))) })
      result$diagnostics$r_eff <- r_eff
      result <- .fdb_loo_restore_zero(result, zero, S, n)
      k <- loo::pareto_k_values(result)
      threshold <- min(0.7, 1 - 1 / log10(S))
      bad <- setdiff(which(!is.finite(k) | k >= threshold), zero)
      est <- result$estimates
      details$looic <- result
      details$loo_r_eff <- result$diagnostics$r_eff
      details$loo_zero_information_observations <- zero
      details$loo_pareto_threshold <- threshold
      details$loo_problematic_observations <- bad
      rows$looic <- row("looic", est["looic", "Estimate"], est["looic", "SE"],
        est["p_loo", "Estimate"], est["elpd_loo", "Estimate"],
        isTRUE(object$diagnostics$converged) && !length(bad),
        if (length(bad)) "Pareto-k diagnostic threshold not satisfied; no automatic refitting." else "")
    }
  }
  if ("dic" %in% methods) {
    if (!dic_ok) {
      reason <- paste0("DIC requires certified finite posterior means for every parameter; unavailable for: ",
        paste(params[!means_ok %in% TRUE], collapse = ", "), ".")
      issue(reason)
      rows$dic <- row("dic", reason = reason)
    } else {
      dic <- tryCatch({
        theta <- colMeans(as.matrix(object$draws[, params, drop = FALSE]))
        at_mean <- object$.loglik(theta)
        if (!is.numeric(at_mean) || length(at_mean) != n || any(!is.finite(at_mean)))
          stop("DIC requires finite pointwise likelihood values at the posterior mean; no alternative point estimate was substituted.")
        deviance <- -2 * rowSums(ll)
        dbar <- mean(deviance)
        dhat <- -2 * sum(at_mean)
        pd <- dbar - dhat
        value <- dbar + pd
        if (any(!is.finite(c(deviance, dbar, dhat, pd, value))))
          stop("DIC arithmetic exceeded the numerical range.")
        list(Dbar = dbar, D_at_mean = dhat, p_D = pd, posterior_mean = theta,
          definition = "DIC = 2 * E[D(theta)] - D(E[theta]); D = -2 * sum(observed log likelihood)")
      }, error = function(e) e)
      if (inherits(dic, "error")) {
        reason <- paste0("DIC unavailable: ", conditionMessage(dic))
        issue(reason); rows$dic <- row("dic", reason = reason)
      } else {
        if (dic$p_D < 0) issue("DIC has negative p_D; inspect posterior geometry and parameterization before interpreting it.")
        details$dic <- dic
        rows$dic <- row("dic", dic$Dbar + dic$p_D, p = dic$p_D,
          reliable = isTRUE(object$diagnostics$converged) && dic$p_D >= 0,
          reason = "Parameterization-dependent; no standard predictive SE. Finite E[D(theta)] is additionally assumed.")
      }
    }
  }
  ans <- list(estimates = do.call(rbind, rows[methods]), details = details,
    diagnostics = list(sampling = object$diagnostics, loo_propriety = guard, warnings = notes),
    model = object$model$name, prior = object$prior$key, nobs = n, ndraws = S,
    observation = .fdb_criteria_identity(object))
  rownames(ans$estimates) <- NULL
  class(ans) <- "fitdistrBayes_criteria"
  ans
}

criteria <- function(object, methods = c("waic", "looic", "dic"), cores = 1L) {
  .fdb_compute_criteria(object, methods, cores)
}

WAIC <- function(object) criteria(object, "waic")
LOOIC <- function(object, cores = 1L) criteria(object, "looic", cores = cores)
DIC <- function(object) criteria(object, "dic")

print.fitdistrBayes_criteria <- function(x, digits = 4L, ...) {
  cat("\nBayesian information criteria (observed-data likelihood)\n")
  cat("Model:", x$model, "| Prior:", x$prior, "| Observations:", x$nobs,
      "| Posterior draws:", x$ndraws, "\n")
  print(x$estimates[, setdiff(names(x$estimates), "reason"), drop = FALSE],
        row.names = FALSE, digits = digits, ...)
  reasons <- unique(x$estimates$reason[nzchar(x$estimates$reason)])
  if (length(reasons)) cat(paste0("\n", reasons, collapse = "\n"), "\n")
  if (length(x$diagnostics$warnings)) cat("\nInspect $diagnostics$warnings and $details before model comparison.\n")
  cat("Lower IC is better; this is not a test or a posterior model probability.\n")
  invisible(x)
}

as.data.frame.fitdistrBayes_criteria <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(x$estimates, row.names = row.names, optional = optional, ...)
}

compare_models <- function(..., criterion = "looic", cores = 1L) {
  methods <- .fdb_criteria_methods(criterion)
  if (length(methods) != 1L) stop("Select exactly one comparison criterion.", call. = FALSE)
  objects <- list(...)
  if (length(objects) < 2L) stop("Supply at least two fits or criteria objects.", call. = FALSE)
  labels <- names(objects)
  if (is.null(labels)) labels <- rep("", length(objects))
  labels[!nzchar(labels)] <- paste0("model", which(!nzchar(labels)))
  if (anyDuplicated(labels)) stop("Supply distinct model labels.", call. = FALSE)
  identity <- lapply(objects, function(z) if (inherits(z, "fitdistrBayes_criteria")) z$observation else {
    if (!inherits(z, c("fitdistrBayes", "fitcensBayes"))) stop("Invalid comparison object.", call. = FALSE)
    .fdb_criteria_identity(z)
  })
  if (!all(vapply(identity, identical, logical(1), identity[[1L]])))
    stop("Comparison requires identical retained data, order, censoring indicators, omitted rows, and likelihood measure. Use the same data scale.", call. = FALSE)
  scores <- lapply(objects, function(z) {
    out <- if (inherits(z, "fitdistrBayes_criteria")) z else criteria(z, methods, cores)
    if (!methods %in% out$estimates$criterion) stop("A criteria object does not contain the requested criterion.", call. = FALSE)
    out
  })
  rows <- lapply(scores, function(z) z$estimates[z$estimates$criterion == methods, , drop = FALSE])
  tab <- do.call(rbind, rows)
  if (any(!tab$available)) stop("The requested criterion is unavailable for at least one model; inspect its criteria result.", call. = FALSE)
  if (any(!tab$reliable)) warning("At least one comparison has diagnostic concerns. The ordering is descriptive, not a reliable model-selection conclusion.", call. = FALSE)
  best <- which.min(tab$estimate)
  tab$model <- labels
  tab$delta <- tab$estimate - tab$estimate[best]
  tab$se_delta <- NA_real_
  if (methods != "dic") {
    pointwise <- lapply(scores, function(z) z$details[[methods]]$pointwise[, methods])
    tab$se_delta <- vapply(pointwise, function(z) .fdb_ic_se(z - pointwise[[best]]), numeric(1))
  }
  tab <- tab[order(tab$estimate), c("model", "criterion", "estimate", "delta", "se_delta", "se", "p_eff", "reliable")]
  rownames(tab) <- NULL
  ans <- list(table = tab, results = stats::setNames(scores, labels),
              criterion = methods, reference = labels[best])
  class(ans) <- "fitdistrBayes_comparison"
  ans
}

print.fitdistrBayes_comparison <- function(x, digits = 4L, ...) {
  cat("\nObserved-data model comparison (lower is better)\n")
  print(x$table, row.names = FALSE, digits = digits, ...)
  cat("\nDelta and its paired SE are relative to:", x$reference, "\n")
  if (any(!x$table$reliable)) cat("Diagnostic concerns: inspect $results before interpreting the ordering.\n")
  invisible(x)
}
