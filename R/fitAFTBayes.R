.aft_control <- function(control) {
  a <- list(rhat_threshold = 1.01, ess_threshold = 400, warn_convergence = TRUE,
    store_callables = TRUE, slice_width = 1, slice_steps = 40L, max_shrink = 1000L,
    optimize_maxit = 1000L)
  if (!is.list(control) || (length(control) && (is.null(names(control)) || anyNA(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control))))) stop("'control' must be a uniquely named list.")
  if (length(setdiff(names(control), names(a)))) stop("Unknown AFT control option.")
  a[names(control)] <- control
  for (nm in c("rhat_threshold", "ess_threshold", "slice_width"))
    if (!is.numeric(a[[nm]]) || length(a[[nm]]) != 1L || !is.finite(a[[nm]]) || a[[nm]] <= 0)
      stop("Invalid AFT control$", nm, ".")
  if (a$rhat_threshold < 1) stop("rhat_threshold must be at least 1.")
  a$optimize_maxit <- .fdb_scalar_count(a$optimize_maxit, "optimize_maxit")
  a$slice_steps <- .fdb_scalar_count(a$slice_steps, "slice_steps")
  a$max_shrink <- .fdb_scalar_count(a$max_shrink, "max_shrink")
  for (nm in c("warn_convergence", "store_callables"))
    if (!is.logical(a[[nm]]) || length(a[[nm]]) != 1L || is.na(a[[nm]])) stop("control$", nm, " must be TRUE or FALSE.")
  a
}

.aft_mode <- function(target, initial, control, center_prob) {
  d <- length(initial); fn <- function(v) {a <- target(v); if (is.finite(a)) -a else 1e100}
  fits <- lapply(c(-1, 0, 1), function(q) {
    v <- initial; v[d] <- asinh(q); v[d-1L] <- initial[d-1L] + .aft_logcosh(v[d])
    v[1L] <- v[1L]+exp(v[d-1L])*.aft_quantile_offset(q, center_prob)
    tryCatch(stats::optim(v, fn, method = "BFGS", control = list(maxit = control$optimize_maxit, reltol = 1e-10)), error = function(e) NULL)
  })
  fits <- Filter(function(f) !is.null(f) && is.finite(f$value) && f$value < 1e90, fits)
  if (!length(fits)) stop("AFT posterior initialization failed at all starting values.")
  f <- fits[[which.min(vapply(fits, function(z) z$value, 0.0))]]
  H <- tryCatch(stats::optimHess(f$par, fn), error = function(e) NULL)
  if (is.null(H) || any(!is.finite(H))) stop("AFT proposal Hessian is nonfinite; inspect data scaling and prior choices.")
  H <- (H+t(H))/2; ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0 || max(ev)/min(ev) > 1e12) stop("AFT proposal Hessian is not positive definite or is ill-conditioned; no sampling was performed.")
  list(center = f$par, H = H, V = chol2inv(chol(H)), convergence = f$convergence,
       objective = f$value, starts = lapply(fits, function(z) z[c("par", "value", "convergence")]))
}

fitAFTBayes <- function(formula, data, distr = "gg", prior = "reference",
                        iter = 4000L, warmup = floor(iter/2), thin = 1L, chains = 4L,
                        seed = NULL, na.action = c("fail", "omit"), control = list(), ...,
                        criteria = FALSE, df = NULL, shape = NULL) {
  call <- match.call()
  if (length(list(...))) stop("Unused arguments in fitAFTBayes().")
  family <- .aft_distribution(distr, df, shape)
  if (family$key != "gg") return(.aft_fit_standard(formula, data, family, prior,
    iter, warmup, thin, chains, seed, na.action, control, criteria, call))
  if (!inherits(formula, "formula") || length(formula) != 3L) stop("Supply a two-sided regression formula.")
  na.action <- match.arg(na.action); control <- .aft_control(control)
  iter <- .fdb_scalar_count(iter, "iter", 20L); warmup <- .fdb_scalar_count(warmup, "warmup", 0L)
  thin <- .fdb_scalar_count(thin, "thin"); chains <- .fdb_scalar_count(chains, "chains", 2L)
  if (warmup >= iter || floor((iter-warmup)/thin) < 10L) stop("Keep at least 10 draws per chain after warmup/thinning.")
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed))) stop("Invalid seed.")
  # model.frame retains prediction metadata for transformations such as poly().
  mf <- stats::model.frame(formula, data = data, na.action = if (na.action == "fail") stats::na.fail else stats::na.omit)
  tt <- attr(mf, "terms")
  if (attr(tt, "intercept") != 1L) stop("The GG AFT interface requires an intercept for its propriety checks.")
  if (length(attr(tt, "offset"))) stop("Offsets are not supported by this AFT interface.")
  response <- stats::model.response(mf)
  if (inherits(response, "Surv")) {
    if (!identical(attr(response, "type"), "right")) stop("Only ordinary right-censored Surv responses are supported; left truncation, interval censoring and multi-state responses are not supported.")
    x <- as.numeric(response[, 1L]); status <- as.integer(response[, 2L])
  } else {
    if (!is.numeric(response) || is.complex(response) || !is.null(dim(response))) stop("Response must be numeric complete times or survival::Surv(time, status).")
    x <- as.numeric(response); status <- rep(1L, length(x))
  }
  if (!length(x) || any(!is.finite(x)) || any(x < 0) || any(x[status == 1L] <= 0) ||
      anyNA(status) || any(!status %in% 0:1)) stop("AFT times must be finite and nonnegative; failure times must be strictly positive, with status 0/1.")
  X <- stats::model.matrix(tt, mf)
  if (any(!is.finite(X))) stop("The AFT design contains nonfinite values.")
  original_X <- X
  design_scales <- apply(abs(X), 2L, max)
  design_scales[design_scales == 0] <- 1
  X <- sweep(X, 2L, design_scales, "/")
  p <- ncol(X); pr <- .aft_prior(prior)
  y <- log(x); check <- .aft_propriety(y, status, X)
  if (!check$ok) stop("AFT posterior propriety is not certified: ", check$reason,
    " These are sufficient conditions; failure does not prove impropriety for every censored sample.", call. = FALSE)
  if (!identical(criteria, FALSE)) .fdb_criteria_preflight(criteria, control)
  positive <- x > 0
  initial_fit <- if (p == 1L) survival::survreg(survival::Surv(x[positive], status[positive]) ~ 1, dist = "lognormal") else
    survival::survreg(survival::Surv(x[positive], status[positive]) ~ X[positive, -1L, drop = FALSE], dist = "lognormal")
  beta <- as.numeric(stats::coef(initial_fit))
  initial <- c(beta, log(initial_fit$scale), 0)
  # Fixed numerical centering probability; not a prior parameter. Once the
  # data are given, the intercept translation has unit Jacobian.
  center_prob <- max(.05, min(.5, mean(status)/2))
  target <- function(v) .aft_logtarget(.aft_uncenter(v, center_prob), y, status, X, pr)
  mode <- .aft_mode(target, initial, control, center_prob)
  parameters <- c(paste0("beta:", colnames(X)), "sigma", "Q")
  sampled <- .fdb_with_seed(seed, .aft_sample(target, mode, p, iter, warmup, thin, chains, control, parameters, center_prob))
  sampled$chains <- lapply(sampled$chains, function(ch) {
    ch[, seq_len(p)] <- sweep(ch[, seq_len(p), drop = FALSE], 2L, design_scales, "/")
    if (any(!is.finite(ch))) stop("AFT coefficients exceed the representable range on the original design scale.")
    ch
  })
  posterior <- .fdb_summarize(sampled$chains)
  precision <- .aft_quantile_precision(sampled$chains)
  moments <- data.frame(parameter = parameters, mean_exists = NA, variance_exists = NA,
    note = "Moment existence is not certified; use posterior medians and quantiles.", stringsAsFactors = FALSE)
  moments$mean_exists[p+2L] <- moments$variance_exists[p+2L] <- FALSE
  moments$note[p+2L] <- "E[abs(Q)] diverges under this prior; finite Monte Carlo means do not establish a posterior mean."
  posterior$mean_exists <- moments$mean_exists; posterior$variance_exists <- moments$variance_exists
  posterior$moment_note <- moments$note
  posterior$mean[!moments$mean_exists %in% TRUE] <- NA_real_
  posterior$sd[!moments$variance_exists %in% TRUE] <- NA_real_
  posterior$mcse_mean[!(moments$mean_exists %in% TRUE & moments$variance_exists %in% TRUE)] <- NA_real_
  good <- is.finite(posterior$rhat) & posterior$rhat <= control$rhat_threshold &
    is.finite(posterior$ess_bulk) & posterior$ess_bulk >= control$ess_threshold &
    is.finite(posterior$ess_tail) & posterior$ess_tail >= control$ess_threshold
  message <- if (all(good)) "Rank-normalized R-hat and bulk/tail ESS targets met." else
    "AFT sampling diagnostic targets were not met. Inspect trace plots and increase iterations; a returned fit is not proof of convergence."
  ans <- list(call = call, formula = formula, terms = tt, design = original_X,
    contrasts = attr(original_X, "contrasts"), xlevels = lapply(mf[vapply(mf, function(z) is.factor(z) || is.character(z), TRUE)],
      function(z) if (is.factor(z)) levels(z) else levels(factor(z))),
    model = list(name = "generalized gamma AFT", parameters = parameters, n = length(x), p = p),
    prior = c(pr, list(posterior_propriety = check$reason)), propriety = check,
    data = x, status = status, omitted = as.integer(attr(mf, "na.action")), na.action = attr(mf, "na.action"),
    chains = sampled$chains, draws = .fdb_long_draws(sampled$chains), summary = posterior,
    estimates = stats::setNames(posterior$median, parameters), moment_status = moments,
    quantile_precision = precision,
    initialization = c(mode, list(chain_starts = sampled$starts, design_scales = design_scales)),
    engine = list(algorithm = "quantile-centered mixed-basis slice sampling", iter = iter,
      warmup = warmup, thin = thin, chains = chains, center_probability = center_prob,
      coordinates = "quantile-centered intercept, slopes, log(sigma*cosh(zeta)), zeta=asinh(Q)",
      updates_per_iteration = p+2L),
    diagnostics = list(converged = all(good), exact_or_independent = FALSE, messages = message,
      movement_fraction = sampled$movement, target_evaluations = sampled$evaluations,
      max_quantile_mcse_ratio = if (all(is.finite(precision$max_mcse_interval_ratio)))
        max(precision$max_mcse_interval_ratio) else NA_real_,
      max_rhat = max(posterior$rhat), min_ess_bulk = min(posterior$ess_bulk),
      min_ess_tail = min(posterior$ess_tail)), control = control,
    censoring = list(events = sum(status), censored = sum(status == 0L), fraction = mean(status == 0L)))
  if (control$store_callables) ans$.loglik <- .aft_loglik_factory(original_X, y, status)
  class(ans) <- "fitAFTBayes"
  if (!identical(criteria, FALSE)) ans$criteria <- .fdb_criteria_after_fit(ans, criteria)
  if (!all(good) && control$warn_convergence) warning(message, call. = FALSE)
  ans
}
