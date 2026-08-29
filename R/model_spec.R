# Extensible model specifications ---------------------------------------------

#' Define a user-supplied Bayesian distribution model
#'
#' `fitdistrBayes_model()` creates a self-contained specification that can be
#' passed as the `distr` argument of [fitdistrBayes()].  It is intended for
#' distributions or priors that are not in the built-in objective-prior
#' catalogue.  The specification may use the package's adaptive Metropolis or
#' univariate slice sampler, or a user-supplied posterior sampler.
#'
#' @param density A density function whose first argument is the observation
#'   vector.  Remaining named arguments are model parameters.  A logical
#'   `log` argument is recommended.
#' @param prior A prior-density function.  It may receive the parameters as
#'   named scalar arguments or as a named vector in its first argument.
#' @param start A finite, uniquely named numeric vector of starting values.
#' @param name A short model name used in printed output.
#' @param lower,upper Optional scalar or parameter-specific bounds.  Missing
#'   bounds default to negative or positive infinity, respectively.
#' @param fixed Optional named list of fixed model quantities.
#' @param engine Posterior engine.  `"adaptive_metropolis"` uses the generic
#'   component-wise adaptive sampler, `"slice"` uses the generic univariate
#'   slice sampler, and `"custom"` calls `sampler`.
#' @param sampler For `engine = "custom"`, a function returning posterior
#'   chains on the natural parameter scale.  It receives any matching formal
#'   arguments among `x`, `log_posterior`, `log_posterior_unconstrained`,
#'   `start`, `fixed`, `lower`, `upper`, `to_unconstrained`,
#'   `from_unconstrained`, `iter`, `warmup`, `thin`, `chains`, `n_save`, and
#'   `control`.  It must return a list of chain matrices, or a list containing
#'   a `chains` component.  Each chain must contain `n_save` rows and one column
#'   per parameter.
#' @param independent Logical indicator that draws returned by a custom sampler
#'   are independent.  If true, MCMC convergence diagnostics are marked as not
#'   applicable.
#' @param engine_label Optional label for a custom sampler.
#' @param propriety Optional declaration or check for posterior propriety.  It
#'   may be `NULL`, one logical value, a character explanation, or a function
#'   of `(x, fixed)` returning a logical value or a list with components
#'   `proper` and `message`.  This is user-supplied information and is not a
#'   mathematical certification by the package authors.
#' @param moments Optional posterior-moment information, supplied as a data
#'   frame or a function of `(x, fixed)` returning a data frame.  Required
#'   columns are `parameter`, `mean_exists`, and `variance_exists`; an optional
#'   `note` column explains the conditions.
#' @param rng Optional posterior-predictive generator whose first argument is
#'   the requested sample size and whose remaining arguments are named model
#'   parameters.
#' @param rng_validator Optional function returning whether values produced by
#'   `rng` lie on the required support.
#' @param validate Optional data-support function of `(x, fixed)`.  It should
#'   return `TRUE`; `FALSE` or a character string rejects the data.
#' @param density_is_log,prior_is_log Optional logical declarations used when
#'   the corresponding function does not expose a `log` argument.
#' @param prior_style Whether the prior accepts named scalar arguments, a named
#'   parameter vector, or should be detected automatically.
#' @param prior_label,prior_kernel Descriptive metadata stored in the fitted
#'   object.
#' @param reference Optional bibliographic or methodological note stored with
#'   the model specification.
#'
#' @return An object of class `"fitdistrBayes_model"`.  It stores the density,
#'   prior, parameter support, posterior engine, optional mathematical checks,
#'   predictive generator, and descriptive metadata.  The object contains no
#'   fitted values until it is supplied to [fitdistrBayes()].
#' @export
fitdistrBayes_model <- function(
    density, prior, start, name = "user-defined",
    lower = NULL, upper = NULL, fixed = NULL,
    engine = c("adaptive_metropolis", "slice", "custom"),
    sampler = NULL, independent = FALSE, engine_label = NULL,
    propriety = NULL, moments = NULL,
    rng = NULL, rng_validator = NULL, validate = NULL,
    density_is_log = NULL, prior_is_log = NULL,
    prior_style = c("auto", "scalar", "vector"),
    prior_label = "user-defined", prior_kernel = "user-supplied function",
    reference = NULL) {
  if (!is.function(density)) {
    stop("'density' must be a function.", call. = FALSE)
  }
  if (!is.function(prior)) {
    stop("'prior' must be a function.", call. = FALSE)
  }
  if (!is.numeric(start) || !length(start) || is.null(names(start)) ||
      any(!nzchar(names(start))) || anyDuplicated(names(start)) ||
      any(!is.finite(start))) {
    stop(
      "'start' must be a finite, uniquely named numeric vector.",
      call. = FALSE
    )
  }
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !nzchar(trimws(name))) {
    stop("'name' must be one nonempty character string.", call. = FALSE)
  }
  if (is.null(fixed)) fixed <- list()
  if (!is.list(fixed) || (length(fixed) &&
      (is.null(names(fixed)) || any(!nzchar(names(fixed))) ||
       anyDuplicated(names(fixed))))) {
    stop("'fixed' must be NULL or a uniquely named list.", call. = FALSE)
  }
  overlap <- intersect(names(start), names(fixed))
  if (length(overlap)) {
    stop(
      "Parameters cannot be both unknown and fixed: ",
      paste(overlap, collapse = ", "), ".", call. = FALSE
    )
  }

  engine <- match.arg(engine)
  prior_style <- match.arg(prior_style)
  if (engine == "slice" && length(start) != 1L) {
    stop("The generic slice engine requires exactly one unknown parameter.",
         call. = FALSE)
  }
  if (engine == "custom" && !is.function(sampler)) {
    stop("'sampler' must be supplied when engine = \"custom\".",
         call. = FALSE)
  }
  if (engine != "custom" && !is.null(sampler)) {
    stop("'sampler' is used only when engine = \"custom\".",
         call. = FALSE)
  }
  if (!is.logical(independent) || length(independent) != 1L ||
      is.na(independent)) {
    stop("'independent' must be TRUE or FALSE.", call. = FALSE)
  }
  for (item in c("density_is_log", "prior_is_log")) {
    value <- get(item)
    if (!is.null(value) &&
        (!is.logical(value) || length(value) != 1L || is.na(value))) {
      stop("'", item, "' must be NULL, TRUE, or FALSE.", call. = FALSE)
    }
  }
  for (item in c("rng", "rng_validator", "validate")) {
    value <- get(item)
    if (!is.null(value) && !is.function(value)) {
      stop("'", item, "' must be NULL or a function.", call. = FALSE)
    }
  }
  if (!(is.null(propriety) || is.function(propriety) ||
        (is.logical(propriety) && length(propriety) == 1L &&
         !is.na(propriety)) ||
        (is.character(propriety) && length(propriety) == 1L &&
         !is.na(propriety) && nzchar(propriety)))) {
    stop(
      "'propriety' must be NULL, one logical value, one explanation, or a function.",
      call. = FALSE
    )
  }
  if (!(is.null(moments) || is.function(moments) || is.data.frame(moments))) {
    stop("'moments' must be NULL, a data frame, or a function.",
         call. = FALSE)
  }
  for (item in c("engine_label", "reference")) {
    value <- get(item)
    if (!is.null(value) &&
        (!is.character(value) || length(value) != 1L || is.na(value))) {
      stop("'", item, "' must be NULL or one character string.", call. = FALSE)
    }
  }
  for (item in c("prior_label", "prior_kernel")) {
    value <- get(item)
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value)) {
      stop("'", item, "' must be one nonempty character string.",
           call. = FALSE)
    }
  }

  answer <- list(
    name = trimws(name), density = density, prior = prior,
    start = start, lower = lower, upper = upper, fixed = fixed,
    engine = engine, sampler = sampler, independent = independent,
    engine_label = engine_label, propriety = propriety, moments = moments,
    rng = rng, rng_validator = rng_validator, validate = validate,
    density_is_log = density_is_log, prior_is_log = prior_is_log,
    prior_style = prior_style, prior_label = prior_label,
    prior_kernel = prior_kernel, reference = reference
  )
  class(answer) <- "fitdistrBayes_model"
  answer
}

#' @export
print.fitdistrBayes_model <- function(x, ...) {
  cat("\nfitdistrBayes model specification\n")
  cat("----------------------------------\n")
  cat("Model:       ", x$name, "\n", sep = "")
  cat("Parameters:  ", paste(names(x$start), collapse = ", "), "\n",
      sep = "")
  cat("Prior:       ", x$prior_label, "\n", sep = "")
  cat("Engine:      ", x$engine, "\n", sep = "")
  cat("Propriety:   ",
      if (is.null(x$propriety)) "not supplied" else "user-supplied",
      "\n", sep = "")
  invisible(x)
}

.fdb_call_matching <- function(fun, arguments, label) {
  fml <- names(formals(fun))
  selected <- if ("..." %in% fml) arguments else
    arguments[intersect(names(arguments), fml)]
  tryCatch(
    do.call(fun, selected),
    error = function(e) {
      .fdb_stop("The user-supplied %s failed: %s", label,
                conditionMessage(e))
    }
  )
}

.fdb_extension_validate_data <- function(spec, x, fixed) {
  if (is.null(spec) || !is.function(spec$validate)) return(FALSE)
  result <- .fdb_call_matching(
    spec$validate, list(x = x, fixed = fixed), "data validator"
  )
  if (is.character(result) && length(result) == 1L && !is.na(result)) {
    .fdb_stop("The user-supplied data validator rejected the data: %s",
              result)
  }
  if (!is.logical(result) || length(result) != 1L || is.na(result)) {
    .fdb_stop("The user-supplied data validator must return TRUE, FALSE, or one explanatory character string.")
  }
  if (!result) {
    .fdb_stop("The user-supplied data validator rejected the data.")
  }
  TRUE
}

.fdb_extension_propriety <- function(spec, x, fixed) {
  declaration <- if (is.null(spec)) NULL else spec$propriety
  if (is.null(declaration)) {
    return(list(
      proper = NA,
      declared = FALSE,
      message = "not certified: user responsibility"
    ))
  }
  result <- if (is.function(declaration)) {
    .fdb_call_matching(
      declaration, list(x = x, fixed = fixed), "propriety check"
    )
  } else {
    declaration
  }
  message <- NULL
  if (is.list(result)) {
    if (is.null(result$proper)) {
      .fdb_stop("A user-supplied propriety check list must contain 'proper'.")
    }
    message <- result$message
    result <- result$proper
  }
  if (is.character(result) && length(result) == 1L && !is.na(result) &&
      nzchar(result)) {
    return(list(
      proper = TRUE,
      declared = TRUE,
      message = paste0("user-supplied propriety statement: ", result)
    ))
  }
  if (!is.logical(result) || length(result) != 1L || is.na(result)) {
    .fdb_stop("A user-supplied propriety declaration must resolve to TRUE or FALSE.")
  }
  if (!result) {
    reason <- if (is.character(message) && length(message) == 1L &&
                  !is.na(message) && nzchar(message)) message else
      "the user-supplied propriety condition was not satisfied"
    .fdb_stop("Posterior computation was stopped because %s.", reason)
  }
  explanation <- if (is.character(message) && length(message) == 1L &&
                     !is.na(message) && nzchar(message)) message else
    "condition verified by a user-supplied propriety check"
  list(
    proper = TRUE,
    declared = TRUE,
    message = paste0("user-supplied: ", explanation)
  )
}

.fdb_extension_moments <- function(spec, x, fixed, parameters) {
  value <- if (is.null(spec)) NULL else spec$moments
  if (is.null(value)) return(NULL)
  if (is.function(value)) {
    value <- .fdb_call_matching(
      value, list(x = x, fixed = fixed), "posterior-moment check"
    )
  }
  if (!is.data.frame(value)) {
    .fdb_stop("The user-supplied posterior-moment information must be a data frame.")
  }
  required <- c("parameter", "mean_exists", "variance_exists")
  missing <- setdiff(required, names(value))
  if (length(missing)) {
    .fdb_stop("The posterior-moment data frame is missing: %s.",
              paste(missing, collapse = ", "))
  }
  if (!is.character(value$parameter) || anyNA(value$parameter) ||
      any(!nzchar(value$parameter)) || anyDuplicated(value$parameter)) {
    .fdb_stop("The posterior-moment 'parameter' column must contain unique nonempty names.")
  }
  if (!is.logical(value$mean_exists) ||
      !is.logical(value$variance_exists)) {
    .fdb_stop("The posterior-moment existence columns must be logical and may contain NA.")
  }
  if (!setequal(value$parameter, parameters)) {
    .fdb_stop(
      "Posterior-moment information must contain exactly these parameters: %s.",
      paste(parameters, collapse = ", ")
    )
  }
  if (!"note" %in% names(value)) {
    value$note <- "user-supplied posterior-moment condition"
  }
  if (!is.character(value$note) || anyNA(value$note)) {
    .fdb_stop("The posterior-moment 'note' column must be character without NA.")
  }
  value <- value[match(parameters, value$parameter),
                 c("parameter", "mean_exists", "variance_exists", "note"),
                 drop = FALSE]
  row.names(value) <- NULL
  value
}

.fdb_extension_sampler <- function(
    fun, arguments, parameters, chains, n_save,
    independent, engine_label, initialization) {
  result <- .fdb_call_matching(fun, arguments, "posterior sampler")
  if (is.matrix(result) || is.data.frame(result)) {
    if (chains != 1L) {
      .fdb_stop("A matrix returned by a custom sampler is valid only when chains = 1; return a list of chain matrices otherwise.")
    }
    result <- list(chains = list(result))
  } else if (is.list(result) && is.null(result$chains) &&
             length(result) == chains &&
             all(vapply(result, function(z) is.matrix(z) || is.data.frame(z),
                        logical(1L)))) {
    result <- list(chains = result)
  }
  if (!is.list(result) || !is.list(result$chains) ||
      length(result$chains) != chains) {
    .fdb_stop("A custom sampler must return one posterior matrix per requested chain.")
  }
  result$chains <- lapply(seq_along(result$chains), function(i) {
    z <- as.matrix(result$chains[[i]])
    storage.mode(z) <- "double"
    if (nrow(z) != n_save || ncol(z) != length(parameters)) {
      .fdb_stop(
        "Custom sampler chain %d must have %d rows and %d columns.",
        i, n_save, length(parameters)
      )
    }
    if (any(!is.finite(z))) {
      .fdb_stop("Custom sampler chain %d contains non-finite values.", i)
    }
    if (is.null(colnames(z))) colnames(z) <- parameters
    if (!setequal(colnames(z), parameters) || anyDuplicated(colnames(z))) {
      .fdb_stop("Custom sampler chain %d has invalid parameter names.", i)
    }
    z[, parameters, drop = FALSE]
  })
  result$independent <- if (is.null(result$independent)) independent else
    isTRUE(result$independent)
  result$engine <- if (is.null(result$engine)) engine_label else
    as.character(result$engine)[1L]
  result$acceptance <- result$acceptance %||% NULL
  result$acceptance_warmup <- result$acceptance_warmup %||% NULL
  result$acceptance_all <- result$acceptance_all %||% NULL
  result$initialization <- result$initialization %||% initialization
  result
}

`%||%` <- function(x, y) if (is.null(x)) y else x
