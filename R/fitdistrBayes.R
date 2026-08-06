# fitdistrBayes ---------------------------------------------------------------
#
# Dependency-free implementation of objective Bayesian distribution fitting.
# The public entry point is fitdistrBayes(). All internal helpers use the
# .fdb_ prefix and are deliberately kept outside the exported namespace.

.fdb_stop <- function(..., call. = FALSE) {
  stop(sprintf(...), call. = call.)
}

.fdb_warn <- function(...) {
  warning(sprintf(...), call. = FALSE, immediate. = TRUE)
}

.fdb_merge_control <- function(control) {
  defaults <- list(
    rhat_threshold = 1.01,
    ess_threshold = 400,
    target_accept = 0.44,
    adapt_interval = 50L,
    proposal_scale = 0.35,
    init_jitter = 0.35,
    slice_width = 1,
    slice_steps = 100L,
    entropy_tol = 1e-10,
    entropy_exact_limit = 2000L,
    max_init_tries = 100L,
    warn_convergence = TRUE,
    store_callables = TRUE,
    rng = NULL,
    rng_validator = NULL,
    density_is_log = NULL,
    prior_is_log = NULL,
    prior_style = "auto",
    lower = NULL,
    upper = NULL
  )
  if (is.null(control) || (is.list(control) && length(control) == 0L)) {
    return(defaults)
  }
  if (!is.list(control) || is.null(names(control)) ||
      any(!nzchar(names(control)))) {
    .fdb_stop("'control' must be a named list.")
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    .fdb_stop("Unknown control option%s: %s.",
              if (length(unknown) > 1L) "s" else "",
              paste(unknown, collapse = ", "))
  }
  defaults[names(control)] <- control
  positive_scalars <- c("rhat_threshold", "ess_threshold",
                        "proposal_scale", "init_jitter", "slice_width",
                        "entropy_tol")
  for (nm in positive_scalars) {
    value <- defaults[[nm]]
    if (length(value) != 1L || !is.numeric(value) ||
        !is.finite(value) || value <= 0) {
      .fdb_stop("control$%s must be one positive finite number.", nm)
    }
  }
  if (defaults$entropy_tol >= 0.1) {
    .fdb_stop("control$entropy_tol must be smaller than 0.1.")
  }
  if (defaults$rhat_threshold < 1) {
    .fdb_stop("control$rhat_threshold must be at least 1.")
  }
  if (length(defaults$target_accept) != 1L ||
      !is.numeric(defaults$target_accept) ||
      !is.finite(defaults$target_accept) ||
      defaults$target_accept <= 0 || defaults$target_accept >= 1) {
    .fdb_stop("control$target_accept must lie strictly between 0 and 1.")
  }
  integer_controls <- c("adapt_interval", "slice_steps",
                        "entropy_exact_limit", "max_init_tries")
  for (nm in integer_controls) {
    value <- defaults[[nm]]
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
        value < 1 || value != as.integer(value)) {
      .fdb_stop("control$%s must be one positive integer.", nm)
    }
    defaults[[nm]] <- as.integer(value)
  }
  logical_controls <- c("warn_convergence", "store_callables")
  for (nm in logical_controls) {
    value <- defaults[[nm]]
    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      .fdb_stop("control$%s must be TRUE or FALSE.", nm)
    }
  }
  if (!is.null(defaults$rng) && !is.function(defaults$rng)) {
    .fdb_stop("control$rng must be NULL or a simulation function.")
  }
  if (!is.null(defaults$rng_validator) &&
      !is.function(defaults$rng_validator)) {
    .fdb_stop("control$rng_validator must be NULL or a validation function.")
  }
  for (nm in c("density_is_log", "prior_is_log")) {
    value <- defaults[[nm]]
    if (!is.null(value) &&
        (length(value) != 1L || !is.logical(value) || is.na(value))) {
      .fdb_stop("control$%s must be NULL, TRUE, or FALSE.", nm)
    }
  }
  if (length(defaults$prior_style) != 1L ||
      !is.character(defaults$prior_style) ||
      is.na(defaults$prior_style) ||
      !defaults$prior_style %in% c("auto", "scalar", "vector")) {
    .fdb_stop(
      "control$prior_style must be one of 'auto', 'scalar', or 'vector'."
    )
  }
  defaults
}

.fdb_scalar_count <- function(value, name, lower = 1L,
                              strict_lower = FALSE) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
    .fdb_stop("'%s' must be one finite number.", name)
  }
  if (value > .Machine$integer.max) {
    .fdb_stop("'%s' must not exceed %d.", name, .Machine$integer.max)
  }
  ok <- if (strict_lower) value > lower else value >= lower
  if (!ok || value != as.integer(value)) {
    op <- if (strict_lower) "greater than" else "at least"
    .fdb_stop("'%s' must be an integer %s %s.", name, op, lower)
  }
  as.integer(value)
}

.fdb_validate_x <- function(x, na.action) {
  if (!is.numeric(x) || is.complex(x) || is.matrix(x) ||
      is.data.frame(x)) {
    .fdb_stop("'x' must be a numeric vector.")
  }
  if (!length(x)) {
    .fdb_stop("'x' must contain at least one observation.")
  }
  bad_na <- is.na(x) | is.nan(x)
  omitted <- which(bad_na)
  if (length(omitted)) {
    if (na.action == "fail") {
      .fdb_stop("'x' contains %d missing value%s; use na.action = \"omit\" to remove them.",
                length(omitted), if (length(omitted) == 1L) "" else "s")
    }
    x <- x[!bad_na]
  }
  if (!length(x)) {
    .fdb_stop("No observations remain after removing missing values.")
  }
  if (any(!is.finite(x))) {
    .fdb_stop("'x' must not contain Inf or -Inf.")
  }
  list(x = as.numeric(x), omitted = omitted)
}

.fdb_nonconstant <- function(x) {
  if (length(x) < 2L) return(FALSE)
  # The question is whether at least two supplied floating-point values are
  # distinct, not whether their range is large relative to their location.
  # A relative-to-|x| tolerance falsely labels translated samples such as
  # 1e12 + c(0, 1e-3) as constant, although the values are representable.
  any(x != x[1L])
}

.fdb_min_positive <- .Machine$double.xmin * .Machine$double.eps

.fdb_finite_magnitude <- function(x, what = "A generated value") {
  if (any(!is.finite(x))) {
    .fdb_stop(
      "%s exceeded the finite floating-point range; no clipping was applied.",
      what
    )
  }
  x
}

.fdb_positive_finite <- function(x, what = "A generated positive value") {
  if (any(!is.finite(x)) || any(x <= 0)) {
    .fdb_stop(
      paste0(
        "%s could not be represented as a strictly positive finite number; ",
        "no clipping was applied."
      ),
      what
    )
  }
  x
}

.fdb_open_probability <- function(x) {
  pmax(.fdb_min_positive,
       pmin(1 - .Machine$double.eps / 2, x))
}

.fdb_rgamma_positive <- function(n, shape, rate = 1) {
  .fdb_positive_finite(
    stats::rgamma(n, shape = shape, rate = rate),
    "A Gamma random variate"
  )
}

.fdb_exp_positive <- function(log_x) {
  lower <- log(.fdb_min_positive)
  upper <- log(.Machine$double.xmax)
  if (any(!is.finite(log_x)) || any(log_x < lower | log_x > upper)) {
    .fdb_stop(
      paste0(
        "A positive natural-scale draw exceeded the floating-point range; ",
        "no clipping was applied."
      )
    )
  }
  .fdb_positive_finite(exp(log_x), "A transformed positive draw")
}

.fdb_discrete <- function(x, lower = 0) {
  all(x >= lower & x == trunc(x))
}

.fdb_safe_scale <- function(x) {
  z <- stats::mad(x, constant = 1.4826)
  if (!is.finite(z) || z <= 0) z <- stats::IQR(x) / 1.349
  if (!is.finite(z) || z <= 0) z <- stats::sd(x)
  if (!is.finite(z) || z <= 0) z <- diff(range(x)) / 2
  if (!is.finite(z) || z <= 0) z <- max(1, abs(stats::median(x)))
  z
}

.fdb_second_central_moment <- function(x, center = mean(x)) {
  mean((x - center)^2)
}

.fdb_softplus <- function(x) {
  ans <- numeric(length(x))
  large <- x > 0
  ans[large] <- x[large] + log1p(exp(-x[large]))
  ans[!large] <- log1p(exp(x[!large]))
  ans
}

.fdb_log_bessel_i0 <- function(x) {
  if (any(x < 0) || any(!is.finite(x))) {
    return(rep(NA_real_, length(x)))
  }
  log(besselI(x, nu = 0, expon.scaled = TRUE)) + x
}

.fdb_polylog2_unit <- function(z) {
  if (length(z) != 1L || !is.finite(z) || z < 0 || z > 1) {
    .fdb_stop("Internal dilogarithm evaluation requires one value in [0,1].")
  }
  if (z == 0) return(0)
  if (z == 1) return(pi^2 / 6)
  series <- function(y) {
    total <- 0
    power <- y
    k <- 1L
    repeat {
      add <- power / k^2
      total <- total + add
      if (abs(add) <= 2e-16 * max(1, abs(total)) || k >= 100000L) break
      k <- k + 1L
      power <- power * y
    }
    total
  }
  if (z <= 0.5) return(series(z))
  pi^2 / 6 - log(z) * log1p(-z) - series(1 - z)
}

.fdb_polylog3_unit <- function(z) {
  if (length(z) != 1L || !is.finite(z) || z < 0 || z > 1) {
    .fdb_stop("Internal trilogarithm evaluation requires one value in [0,1].")
  }
  if (z == 0) return(0)
  if (z == 1) return(1.2020569031595942854)
  total <- 0
  power <- z
  k <- 1L
  repeat {
    add <- power / k^3
    total <- total + add
    if (abs(add) <= 2e-15 * max(1, abs(total)) || k >= 100000L) break
    k <- k + 1L
    power <- power * z
  }
  total
}

.fdb_sample_lmoments <- function(x) {
  n <- length(x)
  if (n < 2L) return(c(l1 = mean(x), l2 = NA_real_, tau = NA_real_))
  ordered <- sort(x)
  b0 <- mean(ordered)
  b1 <- sum(((seq_len(n) - 1) / (n - 1)) * ordered) / n
  l2 <- 2 * b1 - b0
  c(l1 = b0, l2 = l2, tau = l2 / b0)
}

.fdb_frechet_lmoment_start <- function(x) {
  lm <- .fdb_sample_lmoments(x)
  tau <- unname(lm["tau"])
  shape <- log(2) / log1p(tau)
  scale_parameter <- unname((lm["l1"] / gamma(1 - 1 / shape))^shape)
  method <- "Frechet L-moments (closed form)"
  if (!is.finite(shape) || shape <= 1 ||
      !is.finite(scale_parameter) || scale_parameter <= 0) {
    q <- stats::quantile(x, c(0.25, 0.75), names = FALSE, type = 8)
    shape <- log(log(4) / log(4 / 3)) / log(q[2L] / q[1L])
    median_x <- stats::median(x)
    scale_parameter <- log(2) * median_x^shape
    method <- paste0(method, "; quantile-matching fallback")
  }
  if (!is.finite(shape) || shape <= 0 ||
      !is.finite(scale_parameter) || scale_parameter <= 0) {
    .fdb_stop("Could not compute a finite Frechet initialization.")
  }
  list(values = c(shape = shape, scale = scale_parameter), method = method,
       sample_lmoments = lm)
}

.fdb_lomax_lmoment_start <- function(x) {
  if (length(x) < 2L) {
    return(list(
      values = c(shape = 2, scale = max(x, .fdb_min_positive)),
      method = "Lomax single-observation scale anchor (L-moments require n >= 2)"
    ))
  }
  lm <- .fdb_sample_lmoments(x)
  tau <- unname(lm["tau"])
  method <- "Lomax L-moments (closed form)"
  if (!is.finite(tau) || tau <= 0.5 || tau >= 1) {
    tau <- min(1 - 1e-6, max(0.5005, tau))
    method <- paste0(method, "; interior boundary correction")
  }
  shape <- tau / (2 * tau - 1)
  scale <- unname(lm["l1"] * (shape - 1))
  if (!is.finite(shape) || shape <= 1 ||
      !is.finite(scale) || scale <= 0) {
    shape <- 2
    scale <- max(mean(x), .fdb_min_positive)
    method <- paste0(method, "; finite moment fallback")
  }
  list(values = c(shape = shape, scale = scale), method = method,
       sample_lmoments = lm)
}

.fdb_nakagami_moment_start <- function(x) {
  y <- x^2
  omega <- mean(y)
  denominator <- mean(y^2) - omega^2
  shape <- omega^2 / denominator
  method <- "Nakagami method of moments (closed form)"
  if (!is.finite(shape) || shape <= 0) {
    shape <- 1e4
    method <- paste0(method, "; finite degenerate-sample limit")
  }
  list(values = c(shape = min(shape, 1e4), spread = omega), method = method)
}

.fdb_el_moment_start <- function(x) {
  mean_x <- mean(x)
  cv2 <- .fdb_second_central_moment(x, mean_x) / mean_x^2
  cv2_el <- function(theta) {
    log_theta <- log(theta)
    d2 <- .fdb_polylog2_unit(1 - theta)
    d3 <- .fdb_polylog3_unit(1 - theta)
    -2 * d3 * log_theta / d2^2 - 1
  }
  lower <- 1e-6
  upper <- 1 - 1e-8
  if (!is.finite(cv2) || cv2 <= 1) {
    theta <- upper
    method <- paste0(
      "Exponential-Logarithmic moments; exponential-boundary correction"
    )
  } else if (cv2 >= cv2_el(lower)) {
    theta <- lower
    method <- paste0(
      "Exponential-Logarithmic moments; heavy-tail boundary correction"
    )
  } else {
    theta <- stats::uniroot(
      function(z) cv2_el(z) - cv2,
      interval = c(lower, upper), tol = 1e-9
    )$root
    method <- "Exponential-Logarithmic method of moments (one-dimensional root)"
  }
  rate <- -.fdb_polylog2_unit(1 - theta) / (mean_x * log(theta))
  if (!is.finite(rate) || rate <= 0) rate <- 1 / mean_x
  list(values = c(theta = theta, rate = rate), method = method)
}

.fdb_rician_moment_start <- function(x) {
  m2 <- mean(x^2)
  m4 <- mean(x^4)
  eta4 <- 2 * m2^2 - m4
  method <- "Rician second/fourth method of moments (closed form)"
  if (is.finite(eta4) && eta4 > 0) {
    eta <- eta4^(1 / 4)
  } else {
    eta <- max(0.05 * sqrt(m2), .fdb_min_positive)
    method <- paste0(method, "; positive boundary correction")
  }
  alpha2 <- (m2 - eta^2) / 2
  if (!is.finite(alpha2) || alpha2 <= 0) {
    alpha2 <- max(m2 / 4, .fdb_min_positive)
    method <- paste0(method, "; scale fallback")
  }
  list(values = c(noncentrality = eta, scale = sqrt(alpha2)), method = method)
}

.fdb_gauss_laguerre <- function(order = 96L) {
  j <- seq_len(order)
  jacobi <- diag(2 * j - 1)
  off <- seq_len(order - 1L)
  jacobi[cbind(j[-1L], j[-order])] <- off
  jacobi[cbind(j[-order], j[-1L])] <- off
  eig <- eigen(jacobi, symmetric = TRUE)
  index <- order(eig$values)
  list(nodes = eig$values[index], weights = eig$vectors[1L, index]^2)
}

.fdb_rician_logq_factory <- function() {
  quadrature <- .fdb_gauss_laguerre(96L)
  log_grid <- seq(log(1e-4), log(25), length.out = 321L)
  q_one <- function(rho) {
    y <- sqrt(2 * rho * quadrature$nodes)
    i0e <- besselI(y, nu = 0, expon.scaled = TRUE)
    i1e <- besselI(y, nu = 1, expon.scaled = TRUE)
    integral <- sum(
      quadrature$weights * 2 * quadrature$nodes *
        exp(-rho / 2 + y) * i1e^2 / i0e
    )
    psi <- integral - rho
    pmax((rho + 1) * psi - rho, .fdb_min_positive)
  }
  log_values <- log(vapply(exp(log_grid), q_one, numeric(1)))
  interpolator <- stats::splinefun(log_grid, log_values, method = "monoH.FC")
  function(log_rho) {
    if (!is.finite(log_rho)) {
      return(if (log_rho > 0) log(0.5) else -Inf)
    }
    if (log_rho < log_grid[1L]) {
      rho <- exp(log_rho)
      correction <- pmax(1 - 3 * rho, .Machine$double.eps)
      return(3 * log_rho - log(4) + log(correction))
    }
    if (log_rho > log_grid[length(log_grid)]) {
      if (log_rho > log(.Machine$double.xmax)) return(log(0.5))
      inverse <- exp(-log_rho)
      q <- 0.5 - 0.75 * inverse - 0.75 * inverse^2 -
        (19 / 8) * inverse^3
      return(log(pmax(q, .fdb_min_positive)))
    }
    unname(interpolator(log_rho))
  }
}

.fdb_el_prior_terms <- function(log_theta) {
  log_theta <- unname(log_theta)
  theta <- exp(log_theta)
  one_minus <- -expm1(log_theta)
  if (!is.finite(theta) || theta <= 0 || one_minus <= 0) {
    return(c(log_dilog = NA_real_, log_zeta = NA_real_,
             log_information_theta = NA_real_))
  }
  dilog <- .fdb_polylog2_unit(one_minus)
  if (theta < 1e-4) {
    log_i <- -log(2) - 2 * log_theta - log(-log_theta)
    log_zeta <- log(pi^2 / 12 - 1 / 4) -
      2 * log_theta - 2 * log(-log_theta)
  } else if (one_minus < 1e-4) {
    log_i <- -log(12)
    log_zeta <- -log(48)
  } else {
    numerator <- 2 * one_minus^2 + (1 - theta^2) * log_theta
    denominator <- 2 * theta^2 * one_minus^2 * log_theta^2
    information_theta <- -numerator / denominator
    cross_term <- (one_minus + theta * log_theta)^2 /
      (4 * theta^2 * one_minus^2 * log_theta^2)
    zeta <- information_theta * (-dilog / log_theta) - cross_term
    if (!is.finite(information_theta) || information_theta <= 0 ||
        !is.finite(zeta) || zeta <= 0) {
      return(c(log_dilog = NA_real_, log_zeta = NA_real_,
               log_information_theta = NA_real_))
    }
    log_i <- log(information_theta)
    log_zeta <- log(zeta)
  }
  c(log_dilog = log(dilog), log_zeta = log_zeta,
    log_information_theta = log_i)
}

.fdb_cauchy_quantile_start <- function(x) {
  q <- stats::quantile(x, c(0.25, 0.5, 0.75), names = FALSE,
                       type = 8)
  scale <- (q[3L] - q[1L]) / 2
  method <- "Cauchy quantile matching (median and half-IQR)"
  if (!is.finite(scale) || scale <= 0) {
    scale <- .fdb_safe_scale(x)
    method <- paste0(method, "; robust-scale fallback")
  }
  list(values = c(location = q[2L], scale = scale), method = method)
}

.fdb_student_start <- function(x, df = NULL) {
  if (!is.null(df) && df > 2) {
    location <- mean(x)
    m2 <- .fdb_second_central_moment(x, location)
    scale <- sqrt(m2 * (df - 2) / df)
    return(list(
      values = c(location = location, scale = scale),
      method = "Student-t method of moments (fixed df)"
    ))
  }

  if (!is.null(df)) {
    q <- stats::quantile(x, c(0.25, 0.5, 0.75), names = FALSE,
                         type = 8)
    theoretical_iqr <- stats::qt(0.75, df = df) -
      stats::qt(0.25, df = df)
    scale <- (q[3L] - q[1L]) / theoretical_iqr
    if (!is.finite(scale) || scale <= 0) scale <- .fdb_safe_scale(x)
    return(list(
      values = c(location = q[2L], scale = scale),
      method = paste0(
        "Student-t quantile matching (fixed df <= 2; moments unavailable)"
      )
    ))
  }

  location <- mean(x)
  centered <- x - location
  m2 <- mean(centered^2)
  m4 <- mean(centered^4)
  excess <- m4 / m2^2 - 3
  if (is.finite(excess) && excess > 0) {
    df_start <- 4 + 6 / excess
    method <- paste0(
      "Student-t method of moments (variance and excess kurtosis)"
    )
  } else {
    df_start <- 100
    method <- paste0(
      "Student-t method of moments; finite normal-limit fallback for df"
    )
  }
  # Very small positive empirical excess kurtosis implies an arbitrarily
  # large moment estimate.  A finite cap represents the normal limit while
  # keeping the transformed MCMC initialization numerically useful.
  df_start <- min(max(df_start, 4 + sqrt(.Machine$double.eps)), 100)
  scale <- sqrt(m2 * (df_start - 2) / df_start)
  if (!is.finite(scale) || scale <= 0) scale <- .fdb_safe_scale(x)
  list(
    values = c(location = location, scale = scale, df = df_start),
    method = method
  )
}

.fdb_weibull_lmoment_start <- function(x) {
  n <- length(x)
  ordered <- sort(x)
  b0 <- mean(ordered)
  b1 <- sum(((seq_len(n) - 1) / (n - 1)) * ordered) / n
  l1 <- b0
  l2 <- 2 * b1 - b0
  tau <- l2 / l1
  shape <- -log(2) / log1p(-tau)
  scale <- l1 / gamma(1 + 1 / shape)
  method <- "Weibull L-moments (closed form)"

  if (!is.finite(shape) || shape <= 0 ||
      !is.finite(scale) || scale <= 0) {
    q <- stats::quantile(x, c(0.25, 0.75), names = FALSE, type = 8)
    shape <- log(log(4) / log(4 / 3)) / log(q[2L] / q[1L])
    scale <- l1 / gamma(1 + 1 / shape)
    method <- paste0(method, "; quantile-matching fallback")
  }
  if (!is.finite(shape) || shape <= 0 ||
      !is.finite(scale) || scale <= 0) {
    .fdb_stop("Could not compute a finite Weibull L-moment initialization.")
  }
  list(values = c(shape = shape, scale = scale), method = method,
       sample_lmoments = c(l1 = l1, l2 = l2, tau = tau))
}

.fdb_normalize_name <- function(x) {
  y <- tolower(trimws(x))
  y <- gsub("[_.]+", "-", y)
  y <- gsub("[[:space:]]+", "-", y)
  y <- gsub("-+", "-", y)
  y
}

.fdb_model_name <- function(distr) {
  key <- .fdb_normalize_name(distr)
  aliases <- c(
    "beta" = "beta",
    "cauchy" = "cauchy",
    "chi-squared" = "chi-squared",
    "chi-square" = "chi-squared",
    "chisquared" = "chi-squared",
    "chisq" = "chi-squared",
    "exponential" = "exponential",
    "exp" = "exponential",
    "exponential-logarithmic" = "exponential-logarithmic",
    "exponential-logarithmic-distribution" = "exponential-logarithmic",
    "el" = "exponential-logarithmic",
    "frechet" = "frechet",
    "frechet-distribution" = "frechet",
    "gamma" = "gamma",
    "geometric" = "geometric",
    "geom" = "geometric",
    "gumbel" = "gumbel",
    "gumbel-maximum" = "gumbel",
    "log-normal" = "lognormal",
    "lognormal" = "lognormal",
    "lnorm" = "lognormal",
    "logistic" = "logistic",
    "lomax" = "lomax",
    "pareto-ii" = "lomax",
    "pareto-type-ii" = "lomax",
    "nakagami" = "nakagami",
    "nakagami-m" = "nakagami",
    "negative-binomial" = "negative binomial",
    "negativebinomial" = "negative binomial",
    "nbinom" = "negative binomial",
    "normal" = "normal",
    "gaussian" = "normal",
    "norm" = "normal",
    "poisson" = "Poisson",
    "rice" = "rician",
    "rician" = "rician",
    "t" = "t",
    "student-t" = "t",
    "student" = "t",
    "weibull" = "weibull"
  )
  if (!key %in% names(aliases)) {
    .fdb_stop(
      "Unknown distribution '%s'. Recognized models are beta, cauchy, chi-squared, exponential, exponential-logarithmic, Frechet, gamma, geometric, Gumbel, log-normal/lognormal, logistic, Lomax, Nakagami-m, negative binomial, normal, Poisson, Rician, t, and weibull.",
      distr
    )
  }
  unname(aliases[[key]])
}

.fdb_prior_name <- function(prior) {
  key <- .fdb_normalize_name(prior)
  aliases <- c(
    "j" = "jeffreys",
    "jeffreys" = "jeffreys",
    "jeffreys-rule" = "jeffreys",
    "joint-jeffreys" = "jeffreys",
    "reference" = "reference",
    "ref" = "reference",
    "r" = "reference",
    "mdi" = "mdi",
    "maximal-data-information" = "mdi",
    "maximal-information" = "mdi",
    "independence-jeffreys" = "independence-jeffreys",
    "jeffreys-independence" = "independence-jeffreys",
    "first-rule" = "first-rule",
    "jeffreys-first-rule" = "first-rule",
    "reference-shape" = "reference-shape",
    "reference-rate" = "reference-rate",
    "reference-theta" = "reference-theta",
    "reference-lambda" = "reference-rate"
  )
  if (!key %in% names(aliases)) {
    .fdb_stop(
      "Unknown objective prior '%s'. Use Jeffreys, Jeffreys first-rule, reference, MDI, independence-Jeffreys, reference-shape, reference-theta, or reference-rate as appropriate.",
      prior
    )
  }
  unname(aliases[[key]])
}

.fdb_route_catalog <- function() {
  rows <- function(model, priors, parameters, fixed, engine, condition) {
    if (length(engine) == 1L) engine <- rep(engine, length(priors))
    if (length(condition) == 1L) condition <- rep(condition, length(priors))
    if (length(engine) != length(priors) ||
        length(condition) != length(priors)) {
      .fdb_stop("Internal route catalogue error for model '%s'.", model)
    }
    data.frame(
      model = rep(model, length(priors)),
      prior = priors,
      parameters = rep(parameters, length(priors)),
      required_fixed = rep(fixed, length(priors)),
      engine = engine,
      posterior_condition = condition,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  do.call(rbind, list(
    rows("beta", c("jeffreys", "reference"), "shape1, shape2", "none",
         "adaptive Metropolis", "n >= 2; 0 < x < 1; nonconstant sample"),
    rows("cauchy", c("jeffreys", "reference", "mdi"), "location, scale", "none",
         "adaptive Metropolis",
         c("largest multiplicity m satisfies 2m < n",
           "largest multiplicity m satisfies 2m < n + 1",
           "largest multiplicity m satisfies 2m < n + 1")),
    rows("chi-squared", c("jeffreys", "reference"), "df", "none",
         "univariate slice sampling", "n >= 1; x > 0"),
    rows("exponential", c("jeffreys", "reference", "mdi"), "rate", "none",
         "exact Gamma posterior", "n >= 1; x >= 0; sum(x) > 0"),
    rows("exponential-logarithmic",
         c("jeffreys", "mdi", "reference-theta", "reference-rate"),
         "theta, rate", "none", "adaptive Metropolis",
         "x > 0; n > 2 for Jeffreys/reference, n >= 1 for MDI"),
    rows("frechet", c("jeffreys", "reference"), "shape, scale", "none",
         "marginal slice + exact conditional Gamma",
         "n >= 2; x > 0; nonconstant sample"),
    rows("gamma", c("jeffreys", "first-rule", "reference-shape", "reference-rate"),
         "shape, rate", "none", "marginal slice + exact conditional Gamma",
         "n >= 2; x > 0; nonconstant sample"),
    rows("geometric", c("jeffreys", "reference", "mdi"), "prob", "none",
         c("exact Beta posterior", "exact Beta posterior", "univariate slice sampling"),
         "n >= 1; nonnegative integer x"),
    rows("gumbel", c("jeffreys", "reference", "mdi"), "location, scale", "none",
         "adaptive Metropolis", "n >= 2; nonconstant sample"),
    rows("lognormal", c("jeffreys", "reference"), "meanlog, sdlog", "none",
         "exact Normal-inverse-Gamma posterior",
         "n >= 2; x > 0; nonconstant log(x)"),
    rows("logistic", c("jeffreys", "reference", "mdi"), "location, scale", "none",
         "adaptive Metropolis", "n >= 2; nonconstant sample"),
    rows("lomax", "jeffreys", "shape, scale", "none", "adaptive Metropolis",
         "n >= 1; x > 0"),
    rows("nakagami-m", c("jeffreys", "reference"), "shape, spread", "none",
         "marginal slice + exact conditional inverse-Gamma",
         "n >= 2; x > 0; nonconstant sample"),
    rows("negative binomial", c("jeffreys", "reference", "mdi"), "mu", "size",
         c("exact transformed-Beta posterior", "exact transformed-Beta posterior",
           "univariate slice sampling with entropy evaluation"),
         "n >= 1; nonnegative integer x; fixed size > 0"),
    rows("normal", c("jeffreys", "reference", "mdi"), "mean, sd", "none",
         "exact Normal-inverse-Gamma posterior", "n >= 2; nonconstant sample"),
    rows("Poisson", c("jeffreys", "reference", "mdi"), "lambda", "none",
         c("exact Gamma posterior", "exact Gamma posterior",
           "univariate slice sampling with entropy evaluation"),
         "n >= 1; nonnegative integer x"),
    rows("rician", "jeffreys", "noncentrality, scale", "none",
         "adaptive Metropolis with Gauss-Laguerre Jeffreys factor",
         "n > 2; x > 0; nonconstant sample"),
    rows("t", c("jeffreys", "reference", "mdi"), "location, scale", "df",
         "adaptive Metropolis",
         c("fixed df > 0; n+2 > 2 and df*(n-m)-m > 0",
           "fixed df > 0; n+1 > 2 and df*(n-m)-m+1 > 0",
           "fixed df > 0; n+1 > 2 and df*(n-m)-m+1 > 0")),
    rows("t", "independence-jeffreys", "location, scale, df", "none",
         "adaptive Metropolis", "n >= 2; nonconstant sample"),
    rows("weibull", c("jeffreys", "reference"), "shape, scale", "none",
         "marginal slice + exact conditional Gamma transform",
         "n >= 2; x > 0; nonconstant sample")
  ))
}

#' List the built-in objective Bayesian fitting routes
#'
#' Returns the machine-readable catalogue used to document and test the
#' available model--prior combinations. The condition column is a concise
#' summary; the fitting function remains the authoritative executable check.
#'
#' @param model Optional recognized distribution name. If supplied, only its
#'   routes are returned. Use `model = "t"` to display both fixed- and
#'   unknown-degrees-of-freedom Student-t routes.
#' @return A data frame with one row for every enabled model--prior route.
#' @export
fitdistrBayes_routes <- function(model = NULL) {
  catalog <- .fdb_route_catalog()
  if (is.null(model)) return(catalog)
  if (!is.character(model) || length(model) != 1L || is.na(model) ||
      !nzchar(trimws(model))) {
    .fdb_stop("'model' must be NULL or one recognized distribution name.")
  }
  key <- .fdb_model_name(model)
  if (identical(key, "nakagami")) key <- "nakagami-m"
  catalog[catalog$model == key, , drop = FALSE]
}

.fdb_moment_status <- function(model, prior, parameters, n, fixed, x = NULL) {
  answer <- data.frame(
    parameter = parameters,
    mean_exists = rep(NA, length(parameters)),
    variance_exists = rep(NA, length(parameters)),
    note = rep("not certified by the current moment registry", length(parameters)),
    stringsAsFactors = FALSE
  )
  certify <- function(result, par, mean, variance, note) {
    index <- match(par, result$parameter)
    result$mean_exists[index] <- mean
    result$variance_exists[index] <- variance
    result$note[index] <- note
    result
  }

  if (model %in% c("beta", "chi-squared", "exponential", "gamma",
                   "geometric", "Poisson")) {
    answer <- certify(
      answer, parameters, TRUE, TRUE,
      "all positive integer posterior moments are finite"
    )
  } else if (model %in% c("normal", "lognormal")) {
    a <- if (prior == "jeffreys") 2 else 1
    df_marginal <- n + a - 2
    answer <- certify(
      answer, parameters,
      df_marginal > 1,
      df_marginal > 2,
      sprintf("absolute posterior moments of order q exist iff q < %g",
              df_marginal)
    )
  } else if (model == "negative binomial") {
    nr <- n * fixed$size
    answer <- certify(
      answer, parameters,
      nr > 1,
      nr > 2,
      sprintf("posterior moments of mu of order q exist iff q < n*size = %g",
              nr)
    )
  } else if (model == "weibull") {
    answer <- certify(
      answer, "shape", TRUE, TRUE,
      "all positive integer posterior moments are finite"
    )
    answer <- certify(
      answer, "scale", FALSE, FALSE,
      "no positive posterior moment of scale exists under this improper objective prior"
    )
  } else if (model == "t" && is.null(fixed$df)) {
    answer <- certify(
      answer, "df", FALSE, FALSE,
      "positive integer posterior moments of df do not exist"
    )
  } else if (model == "frechet") {
    answer <- certify(
      answer, "shape", TRUE, TRUE,
      "all positive posterior moments of shape are finite for a positive nonconstant sample"
    )
    minimum <- min(x)
    log_ratio_product <- sum(log(x / minimum))
    mean_ok <- log(minimum) < log_ratio_product
    variance_ok <- 2 * log(minimum) < log_ratio_product
    answer <- certify(
      answer, "scale", mean_ok, variance_ok,
      paste0(
        "the q-th scale moment exists iff q*log(min(x)) < ",
        "sum(log(x/min(x))); this is data dependent"
      )
    )
  } else if (model == "nakagami") {
    answer <- certify(
      answer, "shape", TRUE, TRUE,
      paste0(
        "all positive posterior moments of shape are finite for ",
        "n >= 2 and a positive nonconstant sample"
      )
    )
    answer <- certify(
      answer, "spread", FALSE, FALSE,
      "no positive spread moment exists on the full shape support (0,infinity)"
    )
  } else if (model == "exponential-logarithmic") {
    answer <- certify(
      answer, parameters, TRUE, TRUE,
      "all positive integer posterior moments are finite"
    )
  }
  answer
}

.fdb_check_fixed <- function(fixed, allowed) {
  if (is.null(fixed)) return(list())
  if (!is.list(fixed) || is.null(names(fixed)) ||
      any(!nzchar(names(fixed))) || anyDuplicated(names(fixed))) {
    .fdb_stop("'fixed' must be NULL or a uniquely named list.")
  }
  bad <- setdiff(names(fixed), allowed)
  if (length(bad)) {
    .fdb_stop("Unknown fixed parameter%s: %s.",
              if (length(bad) > 1L) "s" else "",
              paste(bad, collapse = ", "))
  }
  fixed
}

.fdb_logsumexp <- function(x) {
  z <- max(x)
  if (!is.finite(z)) return(z)
  z + log(sum(exp(x - z)))
}

.fdb_gamma_joint_term <- function(shape) {
  ans <- shape * trigamma(shape) - 1
  large <- shape > 1e5
  if (any(large)) {
    a <- shape[large]
    ans[large] <- 1 / (2 * a) + 1 / (6 * a^2) - 1 / (30 * a^4)
  }
  pmax(ans, .Machine$double.xmin)
}

.fdb_t_B <- function(df) {
  ans <- trigamma(df / 2) - trigamma((df + 1) / 2) -
    2 * (df + 3) / (df * (df + 1)^2)
  # Direct subtraction loses precision in the far tail.  The leading
  # positive expansion below is used only where that cancellation matters.
  large <- df > 1e4 | !is.finite(ans) | ans <= 0
  if (any(large)) {
    v <- df[large]
    ans[large] <- 6 / v^4 - 12 / v^5 + 14 / v^6 -
      12 / v^7 + 22 / v^8 - 60 / v^9
  }
  pmax(ans, .Machine$double.xmin)
}

.fdb_poisson_entropy <- function(lambda, control) {
  if (!is.finite(lambda) || lambda <= 0) {
    return(if (identical(lambda, 0)) 0 else Inf)
  }
  if (lambda > control$entropy_exact_limit) {
    return(0.5 * log(2 * pi * exp(1) * lambda) -
             1 / (12 * lambda) - 1 / (24 * lambda^2) -
             19 / (360 * lambda^3))
  }
  lo <- max(0, stats::qpois(control$entropy_tol / 2, lambda))
  hi <- stats::qpois(1 - control$entropy_tol / 2, lambda)
  k <- seq.int(lo, hi)
  lp <- stats::dpois(k, lambda, log = TRUE)
  p <- exp(lp)
  -sum(p * lp)
}

.fdb_nbinom_entropy <- function(size, mu, control) {
  if (!is.finite(mu) || mu <= 0) {
    return(if (identical(mu, 0)) 0 else Inf)
  }

  # Exact integral representation of the negative-binomial entropy; see
  # Cheraghchi (2019), IEEE Transactions on Information Theory 65, 3999--4009.
  # Numerical integration has essentially constant cost as mu grows, unlike
  # direct truncation of the infinite probability sum.  The u = -log(1-z)
  # form below also removes the endpoint singularities of the z-integral.
  log1p_ratio <- function(a, b) {
    ratio <- a / b
    answer <- log1p(ratio)
    overflow <- is.infinite(ratio)
    if (any(overflow)) {
      aa <- rep_len(a, length(ratio))
      bb <- rep_len(b, length(ratio))
      answer[overflow] <- log(aa[overflow]) - log(bb[overflow])
    }
    answer
  }

  base <- mu * log1p_ratio(size, mu) +
    size * log1p_ratio(mu, size)
  if (size == 1) return(base)  # geometric entropy, available in closed form

  integrand <- function(u) {
    one_minus_exp <- -expm1(-u)
    denominator <- one_minus_exp * u
    numerator <- if (size < 1) {
      exp(-size * u) * expm1(-(1 - size) * u)
    } else {
      exp(-u) * (-expm1(-(size - 1) * u))
    }
    scaled_mean <- one_minus_exp * mu
    bracket <- expm1(-size * log1p_ratio(scaled_mean, size))
    numerator / denominator * bracket
  }

  integral <- stats::integrate(
    integrand, lower = 0, upper = Inf,
    rel.tol = control$entropy_tol, abs.tol = control$entropy_tol,
    subdivisions = 200L, stop.on.error = FALSE
  )
  value <- base + integral$value
  if (!identical(integral$message, "OK") || !is.finite(value) || value < 0) {
    return(Inf)
  }
  value
}

.fdb_with_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    .fdb_stop("'seed' must be NULL or one integer between 0 and %d.",
              .Machine$integer.max)
  }
  set.seed(as.integer(seed))
  force(code)
}

.fdb_slice_one <- function(current, log_density, width, max_steps) {
  log_y <- log_density(current) - stats::rexp(1)
  if (!is.finite(log_y)) {
    .fdb_stop("The slice sampler encountered a non-finite density at its current state.")
  }
  u <- stats::runif(1)
  left <- current - width * u
  right <- left + width
  j <- floor(stats::runif(1, 0, max_steps))
  k <- (max_steps - 1L) - j
  while (j > 0L && log_density(left) > log_y) {
    left <- left - width
    j <- j - 1L
  }
  while (k > 0L && log_density(right) > log_y) {
    right <- right + width
    k <- k - 1L
  }
  repeat {
    proposal <- stats::runif(1, left, right)
    lp <- log_density(proposal)
    if (is.finite(lp) && lp >= log_y) return(proposal)
    if (proposal < current) left <- proposal else right <- proposal
    if (!is.finite(left + right) || right - left <
        .Machine$double.eps * max(1, abs(current))) {
      .fdb_stop("The slice sampler collapsed to a numerically empty interval.")
    }
  }
}

.fdb_slice_chains <- function(log_density, init, iter, warmup, thin, chains,
                              control) {
  n_save <- floor((iter - warmup) / thin)
  out <- vector("list", chains)
  for (ch in seq_len(chains)) {
    current <- init + stats::rnorm(1, 0, control$init_jitter)
    tries <- 0L
    while (!is.finite(log_density(current)) &&
           tries < control$max_init_tries) {
      current <- init + stats::rnorm(1, 0, control$init_jitter)
      tries <- tries + 1L
    }
    if (!is.finite(log_density(current))) {
      .fdb_stop("Could not find a finite initial state for slice chain %d.", ch)
    }
    saved <- numeric(n_save)
    pos <- 0L
    for (i in seq_len(iter)) {
      current <- .fdb_slice_one(
        current, log_density,
        width = control$slice_width,
        max_steps = control$slice_steps
      )
      if (i > warmup && ((i - warmup) %% thin == 0L)) {
        pos <- pos + 1L
        saved[pos] <- current
      }
    }
    out[[ch]] <- saved
  }
  out
}

.fdb_amwg <- function(log_density, init, iter, warmup, thin, chains, control) {
  p <- length(init)
  n_save <- floor((iter - warmup) / thin)
  out <- vector("list", chains)
  acceptance_post <- matrix(0, nrow = chains, ncol = p)
  acceptance_warmup <- matrix(0, nrow = chains, ncol = p)
  acceptance_all <- matrix(0, nrow = chains, ncol = p)
  colnames(acceptance_post) <- names(init)
  colnames(acceptance_warmup) <- names(init)
  colnames(acceptance_all) <- names(init)

  # Work in model-relevant units. For built-in location-scale families the
  # unconstrained location coordinate has the units of the observations,
  # whereas log(scale) is dimensionless. A common absolute proposal standard
  # deviation makes small-scale or large-scale data mix extremely poorly.
  coordinate_scale <- rep(1, p)
  location_index <- match("location", names(init))
  log_scale_index <- match("log_scale", names(init))
  if (!is.na(location_index) && !is.na(log_scale_index)) {
    initial_scale <- exp(init[log_scale_index])
    if (is.finite(initial_scale) && initial_scale > 0) {
      coordinate_scale[location_index] <- initial_scale
    }
  }
  jitter_sd <- pmax(
    .fdb_min_positive, control$init_jitter * coordinate_scale
  )
  initial_proposal_sd <- pmax(
    .fdb_min_positive, control$proposal_scale * coordinate_scale
  )

  for (ch in seq_len(chains)) {
    current <- init + stats::rnorm(p, 0, jitter_sd)
    names(current) <- names(init)
    lp <- log_density(current)
    tries <- 0L
    while (!is.finite(lp) && tries < control$max_init_tries) {
      current <- init + stats::rnorm(p, 0, jitter_sd)
      names(current) <- names(init)
      lp <- log_density(current)
      tries <- tries + 1L
    }
    if (!is.finite(lp)) {
      current <- init
      lp <- log_density(current)
    }
    if (!is.finite(lp)) {
      .fdb_stop("Could not find a finite initial state for Metropolis chain %d.", ch)
    }

    log_sd <- log(initial_proposal_sd)
    block_accept <- numeric(p)
    warmup_accept <- numeric(p)
    post_accept <- numeric(p)
    saved <- matrix(NA_real_, nrow = n_save, ncol = p,
                    dimnames = list(NULL, names(init)))
    pos <- 0L

    for (i in seq_len(iter)) {
      for (j in seq_len(p)) {
        proposal <- current
        proposal[j] <- stats::rnorm(1, current[j], exp(log_sd[j]))
        lp_new <- log_density(proposal)
        accepted <- is.finite(lp_new) &&
          log(stats::runif(1)) < (lp_new - lp)
        if (accepted) {
          current <- proposal
          lp <- lp_new
          block_accept[j] <- block_accept[j] + 1
          if (i <= warmup) {
            warmup_accept[j] <- warmup_accept[j] + 1
          } else {
            post_accept[j] <- post_accept[j] + 1
          }
        }
      }

      if (i <= warmup && i %% control$adapt_interval == 0L) {
        batch <- i / control$adapt_interval
        gain <- min(0.1, 1 / sqrt(batch))
        rate <- block_accept / control$adapt_interval
        log_sd <- log_sd + gain * (rate - control$target_accept)
        log_sd <- pmin(
          log(.Machine$double.xmax) / 2,
          pmax(log(.Machine$double.xmin) / 2, log_sd)
        )
        block_accept[] <- 0
      }

      if (i > warmup && ((i - warmup) %% thin == 0L)) {
        pos <- pos + 1L
        saved[pos, ] <- current
      }
    }
    out[[ch]] <- saved
    acceptance_warmup[ch, ] <- if (warmup > 0L) {
      warmup_accept / warmup
    } else {
      NA_real_
    }
    acceptance_post[ch, ] <- post_accept / (iter - warmup)
    acceptance_all[ch, ] <- (warmup_accept + post_accept) / iter
  }
  list(
    chains = out,
    acceptance = acceptance_post,
    acceptance_warmup = acceptance_warmup,
    acceptance_all = acceptance_all
  )
}

.fdb_split_matrix <- function(chains, parameter) {
  n <- min(vapply(chains, nrow, integer(1)))
  half <- floor(n / 2)
  if (half < 2L) return(NULL)
  ans <- matrix(NA_real_, nrow = half, ncol = 2L * length(chains))
  col <- 0L
  for (ch in chains) {
    z <- ch[seq_len(n), parameter]
    col <- col + 1L
    ans[, col] <- z[seq_len(half)]
    col <- col + 1L
    ans[, col] <- z[n - half + seq_len(half)]
  }
  ans
}

.fdb_rank_normalize <- function(x) {
  n <- length(x)
  stats::qnorm((rank(x, ties.method = "average") - 3 / 8) / (n + 1 / 4))
}

.fdb_basic_rhat <- function(x) {
  n <- nrow(x)
  chain_var <- apply(x, 2, stats::var)
  W <- mean(chain_var)
  between <- stats::var(colMeans(x))
  scale <- max(1, mean(x^2))
  tolerance <- 100 * .Machine$double.eps * scale
  if (!is.finite(W) || !is.finite(between)) return(NA_real_)
  if (W <= tolerance) {
    return(if (between <= tolerance) 1 else Inf)
  }
  B <- n * between
  var_plus <- (n - 1) / n * W + B / n
  sqrt(max(var_plus / W, 1))
}

.fdb_rhat <- function(x) {
  if (is.null(x)) return(NA_real_)
  z <- matrix(.fdb_rank_normalize(as.vector(x)),
              nrow = nrow(x), ncol = ncol(x))
  folded_raw <- abs(x - stats::median(x))
  folded <- matrix(.fdb_rank_normalize(as.vector(folded_raw)),
                   nrow = nrow(x), ncol = ncol(x))
  max(.fdb_basic_rhat(z), .fdb_basic_rhat(folded))
}

.fdb_ess_matrix <- function(x) {
  if (is.null(x)) return(NA_real_)
  n <- nrow(x)
  m <- ncol(x)
  W <- mean(apply(x, 2, stats::var))
  between <- stats::var(colMeans(x))
  scale <- max(1, mean(x^2))
  tolerance <- 100 * .Machine$double.eps * scale
  if (!is.finite(W) || !is.finite(between)) return(NA_real_)
  if (W <= tolerance) {
    return(if (between <= tolerance) n * m else 0)
  }
  B <- n * between
  var_plus <- (n - 1) / n * W + B / n
  if (!is.finite(var_plus) || var_plus <= 0) return(NA_real_)

  max_lag <- min(n - 1L, 1000L)
  acov <- vapply(seq_len(m), function(j) {
    as.numeric(stats::acf(x[, j], lag.max = max_lag,
                          type = "covariance", plot = FALSE,
                          demean = TRUE)$acf)
  }, numeric(max_lag + 1L))
  mean_acov <- rowMeans(acov)
  rho <- 1 - (W - mean_acov) / var_plus
  rho[1L] <- 1

  pair_sums <- numeric()
  k <- 1L
  while ((2L * k) <= length(rho)) {
    pair <- rho[2L * k - 1L] + rho[2L * k]
    if (!is.finite(pair) || pair < 0) break
    pair_sums <- c(pair_sums, pair)
    k <- k + 1L
  }
  if (length(pair_sums) > 1L) {
    for (j in 2:length(pair_sums)) {
      pair_sums[j] <- min(pair_sums[j], pair_sums[j - 1L])
    }
  }
  tau <- if (length(pair_sums)) -1 + 2 * sum(pair_sums) else 1
  tau <- max(tau, 1 / log10(n * m))
  min(n * m, n * m / tau)
}

.fdb_ess <- function(x, type = c("bulk", "tail")) {
  type <- match.arg(type)
  if (is.null(x)) return(NA_real_)
  if (type == "bulk") {
    z <- matrix(.fdb_rank_normalize(as.vector(x)),
                nrow = nrow(x), ncol = ncol(x))
    return(.fdb_ess_matrix(z))
  }
  q <- stats::quantile(x, c(0.05, 0.95), names = FALSE, type = 8)
  low <- (x <= q[1]) * 1
  high <- (x >= q[2]) * 1
  min(.fdb_ess_matrix(low), .fdb_ess_matrix(high))
}

.fdb_summarize <- function(chains, independent = FALSE) {
  parameters <- colnames(chains[[1L]])
  total <- do.call(rbind, chains)
  ans <- lapply(seq_along(parameters), function(j) {
    values <- total[, j]
    q <- stats::quantile(values, c(0.025, 0.25, 0.5, 0.75, 0.975),
                         names = FALSE, type = 8)
    if (independent) {
      rhat <- 1
      ess_mean <- length(values)
      ess_bulk <- length(values)
      ess_tail <- length(values)
    } else {
      split <- .fdb_split_matrix(chains, j)
      rhat <- .fdb_rhat(split)
      ess_mean <- .fdb_ess_matrix(split)
      ess_bulk <- .fdb_ess(split, "bulk")
      ess_tail <- .fdb_ess(split, "tail")
    }
    sdev <- stats::sd(values)
    data.frame(
      parameter = parameters[j],
      mean = mean(values),
      sd = sdev,
      median = q[3],
      mad = stats::mad(values),
      q2.5 = q[1],
      q25 = q[2],
      q75 = q[4],
      q97.5 = q[5],
      mcse_mean = sdev / sqrt(ess_mean),
      rhat = rhat,
      ess_mean = ess_mean,
      ess_bulk = ess_bulk,
      ess_tail = ess_tail,
      row.names = NULL,
      check.names = FALSE
    )
  })
  do.call(rbind, ans)
}

.fdb_long_draws <- function(chains) {
  out <- lapply(seq_along(chains), function(ch) {
    z <- as.data.frame(chains[[ch]], check.names = FALSE)
    data.frame(
      .chain = ch,
      .iteration = seq_len(nrow(z)),
      .draw = (ch - 1L) * nrow(z) + seq_len(nrow(z)),
      z,
      check.names = FALSE
    )
  })
  do.call(rbind, out)
}

.fdb_transform_matrix <- function(chains_u, transform, parameter_names) {
  lapply(chains_u, function(ch) {
    values <- apply(ch, 1L, transform)
    ans <- if (length(parameter_names) == 1L) {
      matrix(values, ncol = 1L)
    } else {
      t(values)
    }
    colnames(ans) <- parameter_names
    ans
  })
}

.fdb_direct_chains <- function(draw_one_chain, chains, n_save) {
  lapply(seq_len(chains), function(i) draw_one_chain(n_save))
}

.fdb_custom_bounds <- function(start, lower, upper) {
  p <- length(start)
  nm <- names(start)
  expand <- function(z, default, label) {
    if (is.null(z)) return(stats::setNames(rep(default, p), nm))
    if (!is.numeric(z) || is.complex(z) || anyNA(z)) {
      .fdb_stop("control$%s must be real numeric and must not contain NA.",
                label)
    }
    if (!is.null(names(z))) {
      if (any(!nzchar(names(z))) || anyDuplicated(names(z))) {
        .fdb_stop("Names in control$%s must be unique and nonempty.", label)
      }
      bad <- setdiff(names(z), nm)
      if (length(bad)) .fdb_stop("Unknown name in control$%s: %s.",
                                 label, paste(bad, collapse = ", "))
      ans <- stats::setNames(rep(default, p), nm)
      ans[names(z)] <- z
      return(ans)
    }
    if (length(z) == 1L) z <- rep(z, p)
    if (length(z) != p) {
      .fdb_stop("control$%s must have length 1 or length(start).", label)
    }
    stats::setNames(z, nm)
  }
  lo <- expand(lower, -Inf, "lower")
  hi <- expand(upper, Inf, "upper")
  if (any(lo >= hi)) .fdb_stop("Every custom lower bound must be below its upper bound.")
  if (any(start <= lo | start >= hi)) {
    .fdb_stop("Every custom starting value must lie strictly inside its bounds.")
  }
  list(lower = lo, upper = hi)
}

.fdb_to_unconstrained <- function(theta, lower, upper) {
  ans <- numeric(length(theta))
  for (j in seq_along(theta)) {
    if (is.finite(lower[j]) && is.finite(upper[j])) {
      ans[j] <- stats::qlogis((theta[j] - lower[j]) /
                               (upper[j] - lower[j]))
    } else if (is.finite(lower[j])) {
      ans[j] <- log(theta[j] - lower[j])
    } else if (is.finite(upper[j])) {
      ans[j] <- log(upper[j] - theta[j])
    } else {
      ans[j] <- theta[j]
    }
  }
  names(ans) <- names(theta)
  ans
}

.fdb_from_unconstrained <- function(u, lower, upper, jacobian = FALSE) {
  theta <- numeric(length(u))
  log_jac <- 0
  for (j in seq_along(u)) {
    if (is.finite(lower[j]) && is.finite(upper[j])) {
      p <- stats::plogis(u[j])
      theta[j] <- lower[j] + (upper[j] - lower[j]) * p
      log_jac <- log_jac + log(upper[j] - lower[j]) +
        stats::plogis(u[j], log.p = TRUE) +
        stats::plogis(u[j], lower.tail = FALSE, log.p = TRUE)
    } else if (is.finite(lower[j])) {
      theta[j] <- lower[j] + exp(u[j])
      log_jac <- log_jac + u[j]
    } else if (is.finite(upper[j])) {
      theta[j] <- upper[j] - exp(u[j])
      log_jac <- log_jac + u[j]
    } else {
      theta[j] <- u[j]
    }
  }
  names(theta) <- names(u)
  if (jacobian) list(theta = theta, log_jacobian = log_jac) else theta
}

.fdb_resolve_log_mode <- function(fun, requested, label) {
  fml <- names(formals(fun))
  explicit <- "log" %in% fml
  if (is.null(requested)) return(explicit)
  if (isTRUE(requested) && !explicit && !"..." %in% fml) {
    .fdb_stop(
      paste0(
        "control$%s is TRUE, but the custom %s accepts neither an explicit ",
        "'log' argument nor '...'."
      ),
      if (label == "density") "density_is_log" else "prior_is_log",
      label
    )
  }
  isTRUE(requested)
}

.fdb_valid_log_values <- function(ans, expected_length, label) {
  if (!is.numeric(ans) || length(ans) != expected_length ||
      anyNA(ans) || any(is.nan(ans)) || any(ans == Inf)) {
    .fdb_stop(
      paste0(
        "The custom %s must return %s numeric log-density value%s; ",
        "-Inf is allowed outside the support, but +Inf, NA, and NaN are not."
      ),
      label,
      if (expected_length == 1L) "one" else expected_length,
      if (expected_length == 1L) "" else "s"
    )
  }
  ans
}

.fdb_call_density_pointwise <- function(fun, x, theta, fixed, dots,
                                        log_mode) {
  fml <- names(formals(fun))
  args <- c(list(x), as.list(theta), fixed, dots)
  if (log_mode) {
    ans <- do.call(fun, c(args, list(log = TRUE)))
    return(.fdb_valid_log_values(ans, length(x), "density"))
  }
  ordinary_args <- if ("log" %in% fml) {
    c(args, list(log = FALSE))
  } else {
    args
  }
  dens <- do.call(fun, ordinary_args)
  if (!is.numeric(dens) || length(dens) != length(x) ||
      anyNA(dens) || any(is.nan(dens)) || any(!is.finite(dens)) ||
      any(dens < 0)) {
    .fdb_stop("The custom density must return one nonnegative finite density per observation.")
  }
  log(dens)
}

.fdb_call_density <- function(fun, x, theta, fixed, dots, log_mode) {
  sum(.fdb_call_density_pointwise(
    fun, x, theta, fixed, dots, log_mode = log_mode
  ))
}

.fdb_call_prior <- function(fun, theta, log_mode, vector_style) {
  fml <- names(formals(fun))
  if (vector_style) {
    args <- list(theta)
  } else {
    args <- as.list(theta)
  }
  if (log_mode) {
    ans <- do.call(fun, c(args, list(log = TRUE)))
  } else {
    ordinary_args <- if ("log" %in% fml) {
      c(args, list(log = FALSE))
    } else {
      args
    }
    ans <- do.call(fun, ordinary_args)
    if (length(ans) == 1L && is.numeric(ans) && !is.na(ans) &&
        !is.nan(ans) && is.finite(ans) && ans > 0) {
      ans <- log(ans)
    } else if (length(ans) == 1L && is.numeric(ans) &&
               !is.na(ans) && !is.nan(ans) && is.finite(ans) && ans == 0) {
      ans <- -Inf
    } else {
      .fdb_stop("A custom prior without a 'log' argument must return one nonnegative density.")
    }
  }
  .fdb_valid_log_values(ans, 1L, "prior")
}

.fdb_log_contract_equal <- function(log_mode, ordinary_mode,
                                    tolerance = 1e-8) {
  same_infinity <- is.infinite(log_mode) & is.infinite(ordinary_mode) &
    sign(log_mode) == sign(ordinary_mode)
  finite <- is.finite(log_mode) & is.finite(ordinary_mode)
  close <- rep(FALSE, length(log_mode))
  close[finite] <- abs(log_mode[finite] - ordinary_mode[finite]) <=
    tolerance * (1 + abs(ordinary_mode[finite]))
  all(same_infinity | close)
}

.fdb_build_builtin <- function(x, model, prior, fixed, start = NULL) {
  n <- length(x)
  S <- sum(x)
  sum_log_x <- if (all(x > 0)) sum(log(x)) else NA_real_
  m <- max(tabulate(match(x, unique(x))))

  reject_mdi <- function(model_label) {
    .fdb_stop(
      "The MDI prior for the %s model is disabled because its posterior is improper for every finite sample.",
      model_label
    )
  }
  require_prior <- function(allowed, model_label) {
    if (!prior %in% allowed) {
      .fdb_stop("Prior '%s' is not available for %s. Available choice%s: %s.",
                prior, model_label, if (length(allowed) == 1L) " is" else "s are",
                paste(allowed, collapse = ", "))
    }
  }
  initialization_state <- new.env(parent = emptyenv())
  initialization_state$value <- list(
    source = "not applicable",
    method = "independent exact posterior simulation",
    automatic = NULL,
    center = NULL,
    supplied = NULL,
    sampled_parameters = character()
  )
  resolve_start <- function(defaults, method, positive = character(),
                            probability = character(),
                            sampled_parameters = names(defaults)) {
    automatic <- defaults
    supplied <- NULL
    if (!is.null(start)) {
      supplied <- if (is.list(start)) unlist(start, use.names = TRUE) else start
    }
    if (is.null(supplied)) {
      initialization_state$value <- list(
        source = "automatic",
        method = method,
        automatic = automatic,
        center = defaults,
        supplied = NULL,
        sampled_parameters = sampled_parameters
      )
      return(defaults)
    }
    if (!is.numeric(supplied) || any(!is.finite(supplied))) {
      .fdb_stop("'start' must contain finite numeric values.")
    }
    if (is.null(names(supplied))) {
      if (length(supplied) != length(defaults)) {
        .fdb_stop("Unnamed 'start' must have length %d for this model.",
                  length(defaults))
      }
      names(supplied) <- names(defaults)
    }
    if (any(!nzchar(names(supplied))) || anyDuplicated(names(supplied))) {
      .fdb_stop("Named starting values must have unique, nonempty names.")
    }
    bad <- setdiff(names(supplied), names(defaults))
    if (length(bad)) {
      .fdb_stop("Unknown starting parameter%s: %s.",
                if (length(bad) > 1L) "s" else "",
                paste(bad, collapse = ", "))
    }
    defaults[names(supplied)] <- supplied
    if (length(positive) && any(defaults[positive] <= 0)) {
      .fdb_stop("Positive starting value%s required for: %s.",
                if (length(positive) > 1L) "s are" else " is",
                paste(positive, collapse = ", "))
    }
    if (length(probability) &&
        any(defaults[probability] <= 0 | defaults[probability] >= 1)) {
      .fdb_stop("Starting value%s must lie in (0,1): %s.",
                if (length(probability) > 1L) "s" else "",
                 paste(probability, collapse = ", "))
    }
    initialization_state$value <- list(
      source = if (setequal(names(supplied), names(defaults)))
        "user-supplied" else "automatic with user overrides",
      method = method,
      automatic = automatic,
      center = defaults,
      supplied = supplied,
      sampled_parameters = sampled_parameters
    )
    defaults
  }
  make_mh <- function(logpost, init, transform, params, label) {
    function(iter, warmup, thin, chains, control) {
      ans <- .fdb_amwg(logpost, init, iter, warmup, thin, chains, control)
      list(
        chains = .fdb_transform_matrix(ans$chains, transform, params),
        acceptance = ans$acceptance,
        acceptance_warmup = ans$acceptance_warmup,
        acceptance_all = ans$acceptance_all,
        independent = FALSE,
        engine = label,
        initialization = initialization_state$value
      )
    }
  }
  make_slice <- function(logpost, init, transform, params, label) {
    function(iter, warmup, thin, chains, control) {
      raw <- .fdb_slice_chains(logpost, init, iter, warmup, thin,
                               chains, control)
      out <- lapply(raw, function(z) {
        values <- vapply(z, transform, numeric(length(params)))
        ans <- if (length(params) == 1L) {
          matrix(values, ncol = 1L)
        } else {
          t(values)
        }
        colnames(ans) <- params
        ans
      })
      list(chains = out, acceptance = NULL, independent = FALSE,
           acceptance_warmup = NULL, acceptance_all = NULL,
           engine = label, initialization = initialization_state$value)
    }
  }
  make_direct <- function(draw_chain, label) {
    function(iter, warmup, thin, chains, control) {
      n_save <- floor((iter - warmup) / thin)
      list(
        chains = .fdb_direct_chains(draw_chain, chains, n_save),
        acceptance = NULL,
        acceptance_warmup = NULL,
        acceptance_all = NULL,
        independent = TRUE,
        engine = label,
        initialization = initialization_state$value
      )
    }
  }

  if (model == "beta") {
    fixed <- .fdb_check_fixed(fixed, character())
    if (prior == "mdi") reject_mdi("Beta")
    require_prior(c("jeffreys", "reference"), "Beta")
    if (n < 2L || any(x <= 0 | x >= 1) || !.fdb_nonconstant(x)) {
      .fdb_stop(
        "Beta with Jeffreys/reference prior requires n >= 2, 0 < x[i] < 1, and a nonconstant sample."
      )
    }
    mu <- mean(x)
    v <- .fdb_second_central_moment(x, mu)
    common <- mu * (1 - mu) / v - 1
    sv <- resolve_start(
      c(shape1 = max(0.1, mu * common),
        shape2 = max(0.1, (1 - mu) * common)),
      method = "Beta method of moments",
      positive = c("shape1", "shape2")
    )
    beta_mean <- sv["shape1"] / sum(sv)
    concentration <- sum(sv)
    init <- c(logit_mean = stats::qlogis(beta_mean),
              log_concentration = log(concentration))
    logpost <- function(u) {
      beta_mean <- stats::plogis(u[1])
      concentration <- exp(u[2])
      a <- beta_mean * concentration
      b <- (1 - beta_mean) * concentration
      if (!is.finite(a + b) || a <= 0 || b <= 0) return(-Inf)
      ia <- trigamma(a) - trigamma(a + b)
      ib <- trigamma(b) - trigamma(a + b)
      iab <- -trigamma(a + b)
      determinant <- ia * ib - iab^2
      if (!is.finite(determinant) || determinant <= 0) return(-Inf)
      sum(stats::dbeta(x, a, b, log = TRUE)) +
        0.5 * log(determinant) + 2 * u[2] +
        log(beta_mean) + log1p(-beta_mean)
    }
    transform <- function(u) {
      beta_mean <- stats::plogis(u[1])
      concentration <- exp(u[2])
      c(shape1 = beta_mean * concentration,
        shape2 = (1 - beta_mean) * concentration)
    }
    sampler <- make_mh(logpost, init, transform,
                       c("shape1", "shape2"), "adaptive Metropolis")
    loglik <- function(theta) stats::dbeta(x, theta["shape1"],
                                           theta["shape2"], log = TRUE)
    rng <- function(theta, size) stats::rbeta(size, theta["shape1"],
                                              theta["shape2"])
    return(list(
      sampler = sampler, parameters = c("shape1", "shape2"),
      prior_label = if (prior == "reference") "reference (regular block)" else "joint Jeffreys",
      prior_kernel = "sqrt(det(I(shape1, shape2)))",
      propriety = "proper: n >= 2, observations in (0,1), nonconstant sample",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "cauchy") {
    fixed <- .fdb_check_fixed(fixed, character())
    require_prior(c("jeffreys", "reference", "mdi"), "Cauchy")
    a <- if (prior == "jeffreys") 2 else 1
    proper <- if (a == 2) 2 * m < n else 2 * m < n + 1
    if (!proper) {
      .fdb_stop(
        "The Cauchy posterior is improper for this sample/prior. The largest multiplicity is m = %d; %s requires %s.",
        m, prior, if (a == 2) "2m < n" else "2m < n + 1"
      )
    }
    classical <- .fdb_cauchy_quantile_start(x)
    sv <- resolve_start(classical$values,
                        method = classical$method,
                        positive = "scale")
    init <- c(location = unname(sv["location"]),
              log_scale = log(unname(sv["scale"])))
    logpost <- function(u) {
      scale <- exp(u[2])
      if (!is.finite(scale) || scale <= 0) return(-Inf)
      sum(stats::dcauchy(x, u[1], scale, log = TRUE)) -
        a * log(scale) + u[2]
    }
    transform <- function(u) c(location = u[1], scale = exp(u[2]))
    sampler <- make_mh(logpost, init, transform,
                       c("location", "scale"), "adaptive Metropolis")
    loglik <- function(theta) stats::dcauchy(x, theta["location"],
                                             theta["scale"], log = TRUE)
    rng <- function(theta, size) stats::rcauchy(size, theta["location"],
                                                theta["scale"])
    return(list(
      sampler = sampler, parameters = c("location", "scale"),
      prior_label = if (a == 2) "joint Jeffreys" else if (prior == "mdi") "MDI" else "reference",
      prior_kernel = sprintf("scale^(-%d)", a),
      propriety = sprintf("proper for this sample: m = %d and %s", m,
                          if (a == 2) "2m < n" else "2m < n + 1"),
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "chi-squared") {
    fixed <- .fdb_check_fixed(fixed, character())
    if (prior == "mdi") reject_mdi("chi-squared")
    require_prior(c("jeffreys", "reference"), "chi-squared")
    if (any(x <= 0)) {
      .fdb_stop("Chi-squared data must be strictly positive.")
    }
    sv <- resolve_start(c(df = mean(x)),
                        method = "chi-squared method of moments",
                        positive = "df")
    init <- log(sv["df"])
    logpost <- function(u) {
      eta <- u
      df <- exp(eta)
      if (!is.finite(df) || df <= 0) return(-Inf)
      sum(stats::dchisq(x, df, log = TRUE)) +
        0.5 * log(trigamma(df / 2)) + eta
    }
    sampler <- make_slice(logpost, init,
                          function(z) c(df = exp(z)), "df",
                          "univariate slice sampling")
    loglik <- function(theta) stats::dchisq(x, theta["df"], log = TRUE)
    rng <- function(theta, size) stats::rchisq(size, theta["df"])
    return(list(
      sampler = sampler, parameters = "df",
      prior_label = "Jeffreys/reference",
      prior_kernel = "sqrt(trigamma(df/2))",
      propriety = "proper: n >= 1 and all observations are positive",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "exponential") {
    fixed <- .fdb_check_fixed(fixed, character())
    require_prior(c("jeffreys", "reference", "mdi"), "exponential")
    if (any(x < 0) || S <= 0) {
      .fdb_stop("Exponential inference requires x[i] >= 0 and sum(x) > 0.")
    }
    shape_post <- if (prior == "mdi") n + 2 else n
    draw_chain <- function(n_save) {
      ans <- matrix(.fdb_rgamma_positive(
        n_save, shape = shape_post, rate = S
      ),
                    ncol = 1L)
      colnames(ans) <- "rate"
      ans
    }
    sampler <- make_direct(draw_chain, "exact Gamma posterior")
    loglik <- function(theta) stats::dexp(x, theta["rate"], log = TRUE)
    rng <- function(theta, size) stats::rexp(size, theta["rate"])
    return(list(
      sampler = sampler, parameters = "rate",
      prior_label = if (prior == "mdi") "MDI" else "Jeffreys/reference",
      prior_kernel = if (prior == "mdi") "rate" else "1/rate",
      propriety = "proper: n >= 1 and sum(x) > 0",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "exponential-logarithmic") {
    fixed <- .fdb_check_fixed(fixed, character())
    require_prior(
      c("jeffreys", "reference", "reference-theta", "reference-rate", "mdi"),
      "Exponential-Logarithmic"
    )
    if (any(x <= 0)) {
      .fdb_stop("Exponential-Logarithmic data must be strictly positive.")
    }
    if (prior != "mdi" && n <= 2L) {
      .fdb_stop(
        "Exponential-Logarithmic Jeffreys/reference posteriors require n > 2. The MDI route is proper for n >= 1."
      )
    }
    classical <- .fdb_el_moment_start(x)
    sv <- resolve_start(
      classical$values,
      method = classical$method,
      positive = "rate",
      probability = "theta"
    )
    init <- c(log_neglog_theta = log(-log(sv["theta"])),
              log_rate = log(sv["rate"]))
    logpost <- function(u) {
      neg_log_theta <- exp(u[1L])
      rate <- exp(u[2L])
      if (!is.finite(neg_log_theta + rate) ||
          neg_log_theta <= 0 || rate <= 0) return(-Inf)
      log_theta <- -neg_log_theta
      theta <- exp(log_theta)
      one_minus <- -expm1(log_theta)
      if (!is.finite(theta) || theta <= 0 ||
          !is.finite(one_minus) || one_minus <= 0) return(-Inf)
      attenuation <- exp(-rate * x)
      product <- one_minus * attenuation
      if (any(!is.finite(product)) || any(product >= 1)) return(-Inf)
      ll <- n * u[2L] + n * log(one_minus) - n * u[1L] -
        rate * S - sum(log1p(-product))
      if (!is.finite(ll)) return(-Inf)
      terms <- .fdb_el_prior_terms(log_theta)
      if (any(!is.finite(terms))) return(-Inf)
      log_prior <- switch(
        prior,
        "jeffreys" = -u[2L] + 0.5 * terms["log_zeta"],
        "mdi" = u[2L] +
          0.5 * (log(one_minus) - log(theta) - u[1L]) +
          exp(terms["log_dilog"]) / log_theta,
        "reference" = -u[2L] + 0.5 * (
          terms["log_zeta"] + u[1L] - terms["log_dilog"]
        ),
        "reference-theta" = -u[2L] + 0.5 * (
          terms["log_zeta"] + u[1L] - terms["log_dilog"]
        ),
        "reference-rate" = -u[2L] +
          0.5 * terms["log_information_theta"]
      )
      # theta = exp{-exp(u1)} and rate = exp(u2).
      ll + log_prior + log_theta + u[1L] + u[2L]
    }
    transform <- function(u) {
      c(theta = exp(-exp(u[1L])), rate = exp(u[2L]))
    }
    sampler <- make_mh(
      logpost, init, transform, c("theta", "rate"), "adaptive Metropolis"
    )
    loglik <- function(theta) {
      probability <- theta["theta"]
      rate <- theta["rate"]
      log_probability <- log(probability)
      one_minus <- -expm1(log_probability)
      log(rate) + log(one_minus) - log(-log_probability) - rate * x -
        log1p(-one_minus * exp(-rate * x))
    }
    rng <- function(theta, size) {
      probability <- theta["theta"]
      rate <- theta["rate"]
      log_probability <- log(probability)
      u <- .fdb_open_probability(stats::runif(size))
      numerator <- -expm1((1 - u) * log_probability)
      denominator <- -expm1(log_probability)
      .fdb_positive_finite(
        -(log(numerator) - log(denominator)) / rate
      )
    }
    labels <- c(
      "jeffreys" = "joint Jeffreys",
      "mdi" = "MDI",
      "reference" = "reference (theta interest)",
      "reference-theta" = "reference (theta interest)",
      "reference-rate" = "reference (rate interest)"
    )
    kernels <- c(
      "jeffreys" = "sqrt(zeta(theta))/rate",
      "mdi" = "rate*sqrt((1-theta)/(-theta*log(theta)))*exp(dilog(theta)/log(theta))",
      "reference" = "sqrt(-zeta(theta)*log(theta)/dilog(theta))/rate",
      "reference-theta" = "sqrt(-zeta(theta)*log(theta)/dilog(theta))/rate",
      "reference-rate" = "sqrt(I_theta_theta(theta))/rate"
    )
    return(list(
      sampler = sampler, parameters = c("theta", "rate"),
      prior_label = unname(labels[prior]),
      prior_kernel = unname(kernels[prior]),
      propriety = if (prior == "mdi")
        "proper: n >= 1 and all observations are positive; all positive moments are finite" else
        "proper: n > 2 and all observations are positive; all positive moments are finite",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "frechet") {
    fixed <- .fdb_check_fixed(fixed, character())
    if (prior == "mdi") {
      .fdb_stop("An MDI route is not available for Frechet in the certified source registry.")
    }
    require_prior(c("jeffreys", "reference"), "Frechet")
    if (n < 2L || any(x <= 0) || !.fdb_nonconstant(x)) {
      .fdb_stop(
        "Frechet Jeffreys/reference inference requires n >= 2, positive observations, and a nonconstant sample."
      )
    }
    lx <- log(x)
    classical <- .fdb_frechet_lmoment_start(x)
    sv <- resolve_start(
      classical$values,
      method = classical$method,
      positive = c("shape", "scale"),
      sampled_parameters = "shape"
    )
    init_shape <- sv["shape"]
    logpost <- function(u) {
      shape <- exp(u)
      if (!is.finite(shape) || shape <= 0) return(-Inf)
      log_sum_inverse <- .fdb_logsumexp(-shape * lx)
      # Integrating scale from
      # L(shape, scale) / (shape * scale) leaves shape^(n - 1).
      (n - 1) * log(shape) - shape * sum(lx) -
        n * log_sum_inverse + u
    }
    sampler <- function(iter, warmup, thin, chains, control) {
      raw <- .fdb_slice_chains(
        logpost, log(init_shape), iter, warmup, thin, chains, control
      )
      out <- lapply(raw, function(eta) {
        shape <- exp(eta)
        log_sum_inverse <- vapply(
          shape,
          function(a) .fdb_logsumexp(-a * lx),
          numeric(1)
        )
        log_scale <- log(.fdb_rgamma_positive(length(shape), shape = n)) -
          log_sum_inverse
        cbind(shape = shape, scale = .fdb_exp_positive(log_scale))
      })
      list(
        chains = out, acceptance = NULL, independent = FALSE,
        engine = "marginal slice + exact conditional Gamma",
        initialization = initialization_state$value
      )
    }
    loglik <- function(theta) {
      shape <- theta["shape"]
      log_scale <- log(theta["scale"])
      exponent <- log_scale - shape * lx
      base <- log_scale + log(shape) - (shape + 1) * lx
      answer <- rep(-Inf, length(x))
      representable <- is.finite(exponent) & is.finite(base) &
        exponent <= log(.Machine$double.xmax)
      answer[representable] <- base[representable] -
        exp(exponent[representable])
      answer[is.nan(answer) | answer == Inf] <- -Inf
      answer
    }
    rng <- function(theta, size) {
      log_x <- (log(theta["scale"]) -
                  log(.fdb_positive_finite(stats::rexp(size)))) /
        theta["shape"]
      .fdb_exp_positive(log_x)
    }
    return(list(
      sampler = sampler, parameters = c("shape", "scale"),
      prior_label = "Jeffreys/reference",
      prior_kernel = "1/(shape*scale)",
      propriety = "proper: n >= 2, positive nonconstant sample; scale moments are data dependent",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "gumbel") {
    fixed <- .fdb_check_fixed(fixed, character())
    require_prior(c("jeffreys", "reference", "mdi"), "Gumbel")
    if (n < 2L || !.fdb_nonconstant(x)) {
      .fdb_stop(
        "Gumbel objective inference requires n >= 2 and a nonconstant sample."
      )
    }
    location <- mean(x)
    scale <- sqrt(6 * .fdb_second_central_moment(x, location)) / pi
    location <- location - 0.5772156649015329 * scale
    sv <- resolve_start(
      c(location = location, scale = scale),
      method = "Gumbel method of moments (closed form)",
      positive = "scale"
    )
    init <- c(location = unname(sv["location"]),
              log_scale = log(unname(sv["scale"])))
    prior_power <- if (prior == "jeffreys") 2 else 1
    logpost <- function(u) {
      scale <- exp(u[2L])
      if (!is.finite(scale) || scale <= 0) return(-Inf)
      z <- (x - u[1L]) / scale
      if (any(!is.finite(z)) ||
          any(-z > log(.Machine$double.xmax))) return(-Inf)
      ll <- -n * u[2L] - sum(z) - sum(exp(-z))
      ll - prior_power * u[2L] + u[2L]
    }
    transform <- function(u) c(location = u[1L], scale = exp(u[2L]))
    sampler <- make_mh(
      logpost, init, transform, c("location", "scale"),
      "adaptive Metropolis"
    )
    loglik <- function(theta) {
      z <- (x - theta["location"]) / theta["scale"]
      -log(theta["scale"]) - z - exp(-z)
    }
    rng <- function(theta, size) {
      theta["location"] - theta["scale"] *
        log(-log(.fdb_open_probability(stats::runif(size))))
    }
    return(list(
      sampler = sampler, parameters = c("location", "scale"),
      prior_label = if (prior == "jeffreys") "joint Jeffreys" else
        if (prior == "mdi") "MDI" else "reference",
      prior_kernel = if (prior == "jeffreys") "scale^(-2)" else "scale^(-1)",
      propriety = "proper: n >= 2 and nonconstant sample",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "gamma") {
    fixed <- .fdb_check_fixed(fixed, character())
    if (prior == "mdi") reject_mdi("Gamma")
    if (prior == "reference") {
      .fdb_stop(
        "Gamma reference priors depend on the parameter of interest; use 'reference-shape' or 'reference-rate'."
      )
    }
    require_prior(c("jeffreys", "first-rule",
                    "reference-shape", "reference-rate"), "Gamma")
    if (n < 2L || any(x <= 0) || !.fdb_nonconstant(x)) {
      .fdb_stop(
        "Gamma objective posteriors require n >= 2, positive observations, and a nonconstant sample."
      )
    }
    mean_x <- mean(x)
    var_x <- .fdb_second_central_moment(x, mean_x)
    init_shape <- mean_x^2 / var_x
    sv <- resolve_start(c(shape = init_shape, rate = init_shape / mean_x),
                        method = "Gamma method of moments",
                        positive = c("shape", "rate"),
                        sampled_parameters = "shape")
    init_shape <- sv["shape"]
    log_h <- switch(
      prior,
      "jeffreys" = function(a) 0.5 * log(.fdb_gamma_joint_term(a)),
      "first-rule" = function(a) -log(a),
      "reference-shape" = function(a) {
        0.5 * (log(.fdb_gamma_joint_term(a)) - log(a))
      },
      "reference-rate" = function(a) 0.5 * log(trigamma(a))
    )
    logpost_eta <- function(eta) {
      shape <- exp(eta)
      if (!is.finite(shape) || shape <= 0) return(-Inf)
      log_h(shape) + lgamma(n * shape) - n * lgamma(shape) +
        (shape - 1) * sum_log_x - n * shape * log(S) + eta
    }
    sampler <- function(iter, warmup, thin, chains, control) {
      raw <- .fdb_slice_chains(logpost_eta, log(init_shape), iter, warmup,
                               thin, chains, control)
      out <- lapply(raw, function(eta) {
        shape <- exp(eta)
        rate <- .fdb_rgamma_positive(
          length(shape), shape = n * shape, rate = S
        )
        cbind(shape = shape, rate = rate)
      })
      list(chains = out, acceptance = NULL, independent = FALSE,
           engine = "marginal slice + exact conditional Gamma",
           initialization = initialization_state$value)
    }
    loglik <- function(theta) stats::dgamma(x, shape = theta["shape"],
                                            rate = theta["rate"], log = TRUE)
    rng <- function(theta, size) .fdb_rgamma_positive(
      size, shape = theta["shape"], rate = theta["rate"]
    )
    labels <- c(
      "jeffreys" = "joint Jeffreys",
      "first-rule" = "Jeffreys' first-rule prior",
      "reference-shape" = "reference (shape interest)",
      "reference-rate" = "reference (rate interest)"
    )
    kernels <- c(
      "jeffreys" = "sqrt(shape*trigamma(shape)-1)/rate",
      "first-rule" = "1/(shape*rate)",
      "reference-shape" = "sqrt((shape*trigamma(shape)-1)/shape)/rate",
      "reference-rate" = "sqrt(trigamma(shape))/rate"
    )
    return(list(
      sampler = sampler, parameters = c("shape", "rate"),
      prior_label = unname(labels[prior]),
      prior_kernel = unname(kernels[prior]),
      propriety = "proper: n >= 2, positive nonconstant sample",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "geometric") {
    fixed <- .fdb_check_fixed(fixed, character())
    require_prior(c("jeffreys", "reference", "mdi"), "geometric")
    if (!.fdb_discrete(x, 0)) {
      .fdb_stop("Geometric data must be nonnegative integers (R convention: 0,1,...).")
    }
    x <- round(x)
    S <- sum(x)
    if (prior != "mdi") {
      draw_chain <- function(n_save) {
        probability <- .fdb_open_probability(
          stats::rbeta(n_save, n, S + 0.5)
        )
        ans <- matrix(probability, ncol = 1L)
        colnames(ans) <- "prob"
        ans
      }
      sampler <- make_direct(draw_chain, "exact Beta posterior")
    } else {
      prob_mom <- 1 / (1 + mean(x))
      geometric_method <- "geometric method of moments"
      if (prob_mom >= 1) {
        prob_mom <- (n + 0.5) / (n + 1)
        geometric_method <- paste0(
          geometric_method, "; interior continuity correction"
        )
      }
      sv <- resolve_start(c(prob = prob_mom),
                          method = geometric_method,
                          probability = "prob")
      init <- stats::qlogis(sv["prob"])
      logpost <- function(u) {
        eta <- u
        p <- stats::plogis(eta)
        if (!is.finite(p) || p <= 0 || p >= 1) return(-Inf)
        log_p <- log(p)
        log_q <- log1p(-p)
        n * log_p + S * log_q +
          log_p + ((1 - p) / p) * log_q +
          log_p + log_q
      }
      sampler <- make_slice(logpost, init,
                            function(z) c(prob = stats::plogis(z)),
                            "prob", "univariate slice sampling")
    }
    loglik <- function(theta) stats::dgeom(x, theta["prob"], log = TRUE)
    rng <- function(theta, size) stats::rgeom(size, theta["prob"])
    return(list(
      sampler = sampler, parameters = "prob",
      prior_label = if (prior == "mdi") "MDI" else "Jeffreys/reference",
      prior_kernel = if (prior == "mdi")
        "prob*(1-prob)^((1-prob)/prob)" else
        "1/(prob*sqrt(1-prob))",
      propriety = "proper for every n >= 1",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "lognormal") {
    fixed <- .fdb_check_fixed(fixed, character())
    if (prior == "mdi") reject_mdi("lognormal")
    require_prior(c("jeffreys", "reference"), "lognormal")
    if (any(x <= 0)) .fdb_stop("Lognormal data must be strictly positive.")
    y <- log(x)
    if (n < 2L || !.fdb_nonconstant(y)) {
      .fdb_stop("Lognormal Jeffreys/reference inference requires n >= 2 and nonconstant log-data.")
    }
    a <- if (prior == "jeffreys") 2 else 1
    q <- sum((y - mean(y))^2)
    shape_s2 <- (n + a - 2) / 2
    draw_chain <- function(n_save) {
      log_precision <- log(.fdb_rgamma_positive(
        n_save, shape = shape_s2, rate = q / 2
      ))
      s2 <- .fdb_exp_positive(-log_precision)
      meanlog <- .fdb_finite_magnitude(
        stats::rnorm(n_save, mean(y), sqrt(s2 / n))
      )
      cbind(meanlog = meanlog, sdlog = sqrt(s2))
    }
    sampler <- make_direct(draw_chain, "exact Normal-inverse-Gamma posterior")
    loglik <- function(theta) stats::dlnorm(x, theta["meanlog"],
                                            theta["sdlog"], log = TRUE)
    rng <- function(theta, size) stats::rlnorm(size, theta["meanlog"],
                                               theta["sdlog"])
    return(list(
      sampler = sampler, parameters = c("meanlog", "sdlog"),
      prior_label = if (a == 2) "joint Jeffreys" else "reference",
      prior_kernel = sprintf("sdlog^(-%d)", a),
      propriety = "proper: n >= 2 and nonconstant log-data",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "logistic") {
    fixed <- .fdb_check_fixed(fixed, character())
    require_prior(c("jeffreys", "reference", "mdi"), "logistic")
    if (n < 2L || !.fdb_nonconstant(x)) {
      .fdb_stop("Logistic objective inference requires n >= 2 and a nonconstant sample.")
    }
    a <- if (prior == "jeffreys") 2 else 1
    logistic_location <- mean(x)
    logistic_scale <- sqrt(
      3 * .fdb_second_central_moment(x, logistic_location)
    ) / pi
    sv <- resolve_start(
      c(location = logistic_location, scale = logistic_scale),
      method = "logistic method of moments",
      positive = "scale"
    )
    init <- c(location = unname(sv["location"]),
              log_scale = log(unname(sv["scale"])))
    logpost <- function(u) {
      scale <- exp(u[2])
      if (!is.finite(scale) || scale <= 0) return(-Inf)
      sum(stats::dlogis(x, u[1], scale, log = TRUE)) -
        a * log(scale) + u[2]
    }
    transform <- function(u) c(location = u[1], scale = exp(u[2]))
    sampler <- make_mh(logpost, init, transform,
                       c("location", "scale"), "adaptive Metropolis")
    loglik <- function(theta) stats::dlogis(x, theta["location"],
                                            theta["scale"], log = TRUE)
    rng <- function(theta, size) stats::rlogis(size, theta["location"],
                                               theta["scale"])
    return(list(
      sampler = sampler, parameters = c("location", "scale"),
      prior_label = if (a == 2) "joint Jeffreys" else if (prior == "mdi") "MDI" else "reference",
      prior_kernel = sprintf("scale^(-%d)", a),
      propriety = "proper: n >= 2 and nonconstant sample",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "lomax") {
    fixed <- .fdb_check_fixed(fixed, character())
    if (prior == "reference") {
      .fdb_stop(
        "The independent Jeffreys/reference posterior for Lomax is improper for every finite sample. Use 'jeffreys'."
      )
    }
    if (prior == "mdi") {
      .fdb_stop("An MDI route is not available for Lomax in the certified source registry.")
    }
    require_prior("jeffreys", "Lomax")
    if (any(x <= 0)) {
      .fdb_stop(
        "Lomax Jeffreys inference requires strictly positive observations; an exact zero makes this improper-prior route nonregular."
      )
    }
    classical <- .fdb_lomax_lmoment_start(x)
    sv <- resolve_start(
      classical$values,
      method = classical$method,
      positive = c("shape", "scale")
    )
    # The likelihood is nearly ridge-shaped in (log shape, log scale).  The
    # ratio scale/shape is much better identified and substantially reduces
    # posterior correlation without changing the public parameterization.
    init <- c(log_shape = log(sv["shape"]),
              log_scale_per_shape = log(sv["scale"] / sv["shape"]))
    logpost <- function(u) {
      shape <- exp(u[1L])
      log_scale <- u[1L] + u[2L]
      scale <- exp(log_scale)
      if (!is.finite(shape + scale) || shape <= 0 || scale <= 0) {
        return(-Inf)
      }
      log_ratio <- log(x) - log_scale
      log_survival_base <- .fdb_softplus(log_ratio)
      ll <- n * u[1L] - n * log_scale -
        (shape + 1) * sum(log_survival_base)
      log_prior <- -log_scale - 0.5 * u[1L] - log1p(shape) -
        0.5 * log(shape + 2)
      # Jacobian for shape=exp(u1), scale=exp(u1+u2).
      ll + log_prior + u[1L] + log_scale
    }
    transform <- function(u) {
      c(shape = exp(u[1L]), scale = exp(u[1L] + u[2L]))
    }
    sampler <- make_mh(
      logpost, init, transform, c("shape", "scale"),
      "adaptive Metropolis"
    )
    loglik <- function(theta) {
      log(theta["shape"]) - log(theta["scale"]) -
        (theta["shape"] + 1) *
        .fdb_softplus(log(x) - log(theta["scale"]))
    }
    rng <- function(theta, size) {
      u <- .fdb_open_probability(stats::runif(size))
      .fdb_positive_finite(
        theta["scale"] * expm1(-log(u) / theta["shape"])
      )
    }
    return(list(
      sampler = sampler, parameters = c("shape", "scale"),
      prior_label = "joint Jeffreys (dependent)",
      prior_kernel = "1/(scale*sqrt(shape)*(shape+1)*sqrt(shape+2))",
      propriety = "proper: n >= 1 and all observations are strictly positive; reference is improper",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "nakagami") {
    fixed <- .fdb_check_fixed(fixed, character())
    if (prior == "mdi") reject_mdi("Nakagami-m")
    require_prior(c("jeffreys", "reference"), "Nakagami-m")
    if (n < 2L || any(x <= 0) || !.fdb_nonconstant(x)) {
      .fdb_stop(
        paste0(
          "Nakagami-m Jeffreys/reference inference requires n >= 2, ",
          "strictly positive observations, and a nonconstant sample."
        )
      )
    }
    lx <- log(x)
    sum_squares <- sum(x^2)
    classical <- .fdb_nakagami_moment_start(x)
    sv <- resolve_start(
      classical$values,
      method = classical$method,
      positive = c("shape", "spread"),
      sampled_parameters = "shape"
    )
    init_shape <- sv["shape"]
    log_h <- if (prior == "jeffreys") {
      function(shape) 0.5 * log(.fdb_gamma_joint_term(shape))
    } else {
      function(shape) 0.5 * (
        log(.fdb_gamma_joint_term(shape)) - log(shape)
      )
    }
    logpost <- function(u) {
      shape <- exp(u)
      if (!is.finite(shape) || shape <= 0) return(-Inf)
      log_h(shape) + lgamma(n * shape) - n * lgamma(shape) +
        (2 * shape - 1) * sum(lx) - n * shape * log(sum_squares) + u
    }
    sampler <- function(iter, warmup, thin, chains, control) {
      raw <- .fdb_slice_chains(
        logpost, log(init_shape), iter, warmup, thin, chains, control
      )
      out <- lapply(raw, function(eta) {
        shape <- exp(eta)
        gamma_draw <- .fdb_rgamma_positive(length(shape), shape = n * shape)
        log_spread <- log(shape) + log(sum_squares) - log(gamma_draw)
        cbind(shape = shape, spread = .fdb_exp_positive(log_spread))
      })
      list(
        chains = out, acceptance = NULL, independent = FALSE,
        engine = "marginal slice + exact conditional inverse-Gamma",
        initialization = initialization_state$value
      )
    }
    loglik <- function(theta) {
      shape <- theta["shape"]
      spread <- theta["spread"]
      log(2) - lgamma(shape) + shape * (log(shape) - log(spread)) +
        (2 * shape - 1) * lx - shape * x^2 / spread
    }
    rng <- function(theta, size) {
      sqrt(.fdb_rgamma_positive(
        size, shape = theta["shape"],
        rate = theta["shape"] / theta["spread"]
      ))
    }
    return(list(
      sampler = sampler, parameters = c("shape", "spread"),
      prior_label = if (prior == "jeffreys") "joint Jeffreys" else
        "overall reference",
      prior_kernel = if (prior == "jeffreys")
        "sqrt(shape*trigamma(shape)-1)/spread" else
        "sqrt((shape*trigamma(shape)-1)/shape)/spread",
      propriety = paste0(
        "proper: n >= 2 and a positive nonconstant sample; ",
        "MDI is improper"
      ),
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "negative binomial") {
    fixed <- .fdb_check_fixed(fixed, "size")
    if (is.null(fixed$size)) {
      if (prior == "mdi") {
        .fdb_stop("The full negative-binomial MDI posterior is improper. Supply a fixed positive size.")
      }
      .fdb_stop(
        "The exact full-model %s prior is not in the certified registry. Supply fixed = list(size = ...) for the audited one-parameter model.",
        prior
      )
    }
    size <- fixed$size
    if (length(size) != 1L || !is.numeric(size) || !is.finite(size) || size <= 0) {
      .fdb_stop("fixed$size must be one positive finite number.")
    }
    require_prior(c("jeffreys", "reference", "mdi"),
                  "negative binomial with fixed size")
    if (!.fdb_discrete(x, 0)) {
      .fdb_stop("Negative-binomial data must be nonnegative integers.")
    }
    x <- round(x)
    S <- sum(x)
    if (prior != "mdi") {
      draw_chain <- function(n_save) {
        p <- stats::rbeta(n_save, n * size, S + 0.5)
        log_mu <- log(size) + log1p(-p) - log(p)
        mu <- .fdb_exp_positive(log_mu)
        ans <- matrix(mu, ncol = 1L)
        colnames(ans) <- "mu"
        ans
      }
      sampler <- make_direct(draw_chain, "exact transformed-Beta posterior")
    } else {
      mu_mom <- mean(x)
      nb_method <- "negative-binomial method of moments (fixed size)"
      if (mu_mom <= 0) {
        mu_mom <- 0.5 / n
        nb_method <- paste0(nb_method, "; interior continuity correction")
      }
      sv <- resolve_start(c(mu = mu_mom), method = nb_method,
                          positive = "mu")
      init <- log(sv["mu"])
      logpost_nbinom_mdi <- function(u, entropy_control) {
        eta <- u
        mu <- exp(eta)
        if (!is.finite(mu) || mu <= 0) return(-Inf)
        sum(stats::dnbinom(x, size = size, mu = mu, log = TRUE)) -
          .fdb_nbinom_entropy(size, mu, entropy_control) + eta
      }
      sampler <- function(iter, warmup, thin, chains, control) {
        current_logpost <- function(u) logpost_nbinom_mdi(u, control)
        raw <- .fdb_slice_chains(current_logpost, init, iter, warmup, thin,
                                 chains, control)
        out <- lapply(raw, function(z) {
          ans <- matrix(exp(z), ncol = 1L)
          colnames(ans) <- "mu"
          ans
        })
        list(chains = out, acceptance = NULL, independent = FALSE,
             engine = "univariate slice sampling with entropy evaluation",
             initialization = initialization_state$value)
      }
    }
    loglik <- function(theta) stats::dnbinom(x, size = size,
                                             mu = theta["mu"], log = TRUE)
    nb_size <- size
    rng <- function(theta, size) stats::rnbinom(size, size = nb_size,
                                                  mu = theta["mu"])
    return(list(
      sampler = sampler, parameters = "mu",
      prior_label = if (prior == "mdi") "MDI (fixed size)" else "Jeffreys/reference (fixed size)",
      prior_kernel = if (prior == "mdi") "exp(-H_NB(size, mu))" else
        "1/sqrt(mu*(size+mu))",
      propriety = "proper for every n >= 1 with fixed positive size",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "normal") {
    fixed <- .fdb_check_fixed(fixed, character())
    require_prior(c("jeffreys", "reference", "mdi"), "normal")
    if (n < 2L || !.fdb_nonconstant(x)) {
      .fdb_stop("Normal objective inference requires n >= 2 and a nonconstant sample.")
    }
    a <- if (prior == "jeffreys") 2 else 1
    q <- sum((x - mean(x))^2)
    shape_s2 <- (n + a - 2) / 2
    draw_chain <- function(n_save) {
      log_precision <- log(.fdb_rgamma_positive(
        n_save, shape = shape_s2, rate = q / 2
      ))
      s2 <- .fdb_exp_positive(-log_precision)
      mu <- .fdb_finite_magnitude(
        stats::rnorm(n_save, mean(x), sqrt(s2 / n))
      )
      cbind(mean = mu, sd = sqrt(s2))
    }
    sampler <- make_direct(draw_chain, "exact Normal-inverse-Gamma posterior")
    loglik <- function(theta) stats::dnorm(x, theta["mean"],
                                           theta["sd"], log = TRUE)
    rng <- function(theta, size) stats::rnorm(size, theta["mean"],
                                              theta["sd"])
    return(list(
      sampler = sampler, parameters = c("mean", "sd"),
      prior_label = if (a == 2) "joint Jeffreys" else if (prior == "mdi") "MDI" else "reference",
      prior_kernel = sprintf("sd^(-%d)", a),
      propriety = "proper: n >= 2 and nonconstant sample",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "Poisson") {
    fixed <- .fdb_check_fixed(fixed, character())
    require_prior(c("jeffreys", "reference", "mdi"), "Poisson")
    if (!.fdb_discrete(x, 0)) {
      .fdb_stop("Poisson data must be nonnegative integers.")
    }
    x <- round(x)
    S <- sum(x)
    if (prior != "mdi") {
      draw_chain <- function(n_save) {
        ans <- matrix(.fdb_rgamma_positive(
          n_save, shape = S + 0.5, rate = n
        ),
                      ncol = 1L)
        colnames(ans) <- "lambda"
        ans
      }
      sampler <- make_direct(draw_chain, "exact Gamma posterior")
    } else {
      lambda_mom <- mean(x)
      poisson_method <- "Poisson method of moments"
      if (lambda_mom <= 0) {
        lambda_mom <- 0.5 / n
        poisson_method <- paste0(
          poisson_method, "; interior continuity correction"
        )
      }
      sv <- resolve_start(c(lambda = lambda_mom),
                          method = poisson_method,
                          positive = "lambda")
      init <- log(sv["lambda"])
      logpost_poisson_mdi <- function(u, entropy_control) {
        eta <- u
        lambda <- exp(eta)
        if (!is.finite(lambda) || lambda <= 0) return(-Inf)
        S * eta - n * lambda -
          .fdb_poisson_entropy(lambda, entropy_control) + eta
      }
      sampler <- function(iter, warmup, thin, chains, control) {
        current_logpost <- function(u) logpost_poisson_mdi(u, control)
        raw <- .fdb_slice_chains(current_logpost, init, iter, warmup, thin,
                                 chains, control)
        out <- lapply(raw, function(z) {
          ans <- matrix(exp(z), ncol = 1L)
          colnames(ans) <- "lambda"
          ans
        })
        list(chains = out, acceptance = NULL, independent = FALSE,
             engine = "univariate slice sampling with entropy evaluation",
             initialization = initialization_state$value)
      }
    }
    loglik <- function(theta) stats::dpois(x, theta["lambda"], log = TRUE)
    rng <- function(theta, size) stats::rpois(size, theta["lambda"])
    return(list(
      sampler = sampler, parameters = "lambda",
      prior_label = if (prior == "mdi") "MDI" else "Jeffreys/reference",
      prior_kernel = if (prior == "mdi") "exp(-H_Poisson(lambda))" else
        "lambda^(-1/2)",
      propriety = "proper for every n >= 1, including an all-zero sample",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "rician") {
    fixed <- .fdb_check_fixed(fixed, character())
    if (prior == "reference") {
      .fdb_stop(
        "A reference prior for the Rician model is not available in the certified source registry."
      )
    }
    if (prior == "mdi") {
      .fdb_stop(
        "An MDI route for the Rician model is not available in the certified source registry."
      )
    }
    require_prior("jeffreys", "Rician")
    if (n <= 2L || any(x <= 0) || !.fdb_nonconstant(x)) {
      .fdb_stop(
        "Rician Jeffreys inference requires n > 2, positive observations, and a nonconstant sample."
      )
    }
    classical <- .fdb_rician_moment_start(x)
    sv <- resolve_start(
      classical$values,
      method = classical$method,
      positive = c("noncentrality", "scale")
    )
    init <- c(log_noncentrality = log(sv["noncentrality"]),
              log_scale = log(sv["scale"]))
    log_q <- .fdb_rician_logq_factory()
    pointwise_loglik <- function(log_noncentrality, log_scale) {
      noncentrality <- exp(log_noncentrality)
      scale <- exp(log_scale)
      if (!is.finite(noncentrality + scale) ||
          noncentrality <= 0 || scale <= 0) return(NULL)
      log_argument <- log_noncentrality + log(x) - 2 * log_scale
      if (any(!is.finite(log_argument)) ||
          any(log_argument > log(.Machine$double.xmax))) return(NULL)
      argument <- exp(log_argument)
      scaled_i0 <- besselI(argument, nu = 0, expon.scaled = TRUE)
      standardized_difference <- (x - noncentrality) / scale
      if (any(!is.finite(scaled_i0)) || any(scaled_i0 <= 0) ||
          any(!is.finite(standardized_difference))) return(NULL)
      log(x) - 2 * log_scale + log(scaled_i0) -
        0.5 * standardized_difference^2
    }
    logpost <- function(u) {
      ll <- pointwise_loglik(u[1L], u[2L])
      if (is.null(ll) || any(!is.finite(ll))) return(-Inf)
      log_rho <- 2 * (u[1L] - u[2L])
      log_q_value <- log_q(log_rho)
      if (!is.finite(log_q_value)) return(-Inf)
      sum(ll) + 0.5 * log_q_value + u[1L] - u[2L]
    }
    transform <- function(u) {
      c(noncentrality = exp(u[1L]), scale = exp(u[2L]))
    }
    sampler <- make_mh(
      logpost, init, transform, c("noncentrality", "scale"),
      "adaptive Metropolis with Gauss-Laguerre Jeffreys factor"
    )
    loglik <- function(theta) {
      answer <- pointwise_loglik(
        log(theta["noncentrality"]), log(theta["scale"])
      )
      if (is.null(answer)) rep(-Inf, length(x)) else answer
    }
    rng <- function(theta, size) {
      sqrt(
        stats::rnorm(size, theta["noncentrality"], theta["scale"])^2 +
          stats::rnorm(size, 0, theta["scale"])^2
      )
    }
    return(list(
      sampler = sampler, parameters = c("noncentrality", "scale"),
      prior_label = "joint Jeffreys",
      prior_kernel = "sqrt((rho+1)*Psi(rho)-rho)/scale^2, rho=noncentrality^2/scale^2",
      propriety = "proper: n > 2, positive nonconstant sample; reference and MDI routes are unavailable",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "t") {
    fixed <- .fdb_check_fixed(fixed, "df")
    if (!is.null(fixed$df)) {
      df <- fixed$df
      if (length(df) != 1L || !is.numeric(df) || !is.finite(df) || df <= 0) {
        .fdb_stop("fixed$df must be one positive finite number.")
      }
      require_prior(c("jeffreys", "reference", "mdi"),
                    "Student-t with fixed df")
      a <- if (prior == "jeffreys") 2 else 1
      if (!(n + a > 2 && df * (n - m) - m - a + 2 > 0)) {
        .fdb_stop(
          "The fixed-df Student-t posterior is improper: it requires n+a>2 and df*(n-m)-m-a+2>0 (here n=%d, m=%d, a=%d, df=%g).",
          n, m, a, df
        )
      }
      classical <- .fdb_student_start(x, df = df)
      sv <- resolve_start(classical$values,
                          method = classical$method,
                          positive = "scale")
      init <- c(location = unname(sv["location"]),
                log_scale = log(unname(sv["scale"])))
      logpost <- function(u) {
        scale <- exp(u[2])
        if (!is.finite(scale) || scale <= 0) return(-Inf)
        sum(stats::dt((x - u[1]) / scale, df, log = TRUE) -
              log(scale)) - a * log(scale) + u[2]
      }
      transform <- function(u) c(location = u[1], scale = exp(u[2]))
      sampler <- make_mh(logpost, init, transform,
                         c("location", "scale"), "adaptive Metropolis")
      loglik <- function(theta) stats::dt(
        (x - theta["location"]) / theta["scale"], df, log = TRUE
      ) - log(theta["scale"])
      rng <- function(theta, size) theta["location"] +
        theta["scale"] * stats::rt(size, df)
      return(list(
        sampler = sampler, parameters = c("location", "scale"),
        prior_label = if (a == 2) "joint Jeffreys (fixed df)" else if (prior == "mdi") "MDI (fixed df)" else "reference (fixed df)",
        prior_kernel = sprintf("scale^(-%d)", a),
        propriety = sprintf("proper for this sample: n=%d, m=%d, df=%g", n, m, df),
        loglik = loglik, rng = rng, fixed = fixed
      ))
    }

    if (prior == "jeffreys") {
      .fdb_stop("Joint Jeffreys for Student-t with unknown df has an improper posterior. Use 'independence-jeffreys' or fix df.")
    }
    if (prior == "mdi") {
      .fdb_stop("The MDI posterior for Student-t with unknown df is improper. Fix df to use MDI.")
    }
    if (prior %in% c("reference", "reference-shape", "reference-rate")) {
      .fdb_stop("Unknown-df Student-t reference priors are not yet in the certified full-joint registry. Fix df or use 'independence-jeffreys'.")
    }
    require_prior("independence-jeffreys", "Student-t with unknown df")
    if (n <= 1L || !.fdb_nonconstant(x)) {
      .fdb_stop("Unknown-df Student-t independence Jeffreys requires n > 1 and a nondegenerate sample.")
    }
    classical <- .fdb_student_start(x)
    sv <- resolve_start(classical$values,
                        method = classical$method,
                        positive = c("scale", "df"))
    init <- c(location = unname(sv["location"]),
              log_scale = log(unname(sv["scale"])),
              log_df = log(unname(sv["df"])))
    logpost <- function(u) {
      scale <- exp(u[2]); df <- exp(u[3])
      if (!is.finite(scale + df) || scale <= 0 || df <= 0) return(-Inf)
      B <- .fdb_t_B(df)
      log_prior <- -log(scale) +
        0.5 * (log(df) - log(df + 3) + log(B))
      sum(stats::dt((x - u[1]) / scale, df, log = TRUE) -
            log(scale)) + log_prior + u[2] + u[3]
    }
    transform <- function(u) c(location = u[1], scale = exp(u[2]),
                               df = exp(u[3]))
    sampler <- make_mh(logpost, init, transform,
                       c("location", "scale", "df"),
                       "adaptive Metropolis")
    loglik <- function(theta) stats::dt(
      (x - theta["location"]) / theta["scale"], theta["df"], log = TRUE
    ) - log(theta["scale"])
    rng <- function(theta, size) theta["location"] +
      theta["scale"] * stats::rt(size, theta["df"])
    return(list(
      sampler = sampler, parameters = c("location", "scale", "df"),
      prior_label = "independence Jeffreys",
      prior_kernel = "scale^(-1)*sqrt(df*B(df)/(df+3))",
      propriety = "proper for intercept-only model with n > 1 outside degenerate samples; positive df moments do not exist",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  if (model == "weibull") {
    fixed <- .fdb_check_fixed(fixed, character())
    if (prior == "mdi") reject_mdi("Weibull")
    require_prior(c("jeffreys", "reference"), "Weibull")
    if (n < 2L || any(x <= 0) || !.fdb_nonconstant(x)) {
      .fdb_stop("Weibull objective posteriors require n >= 2, positive observations, and a nonconstant sample.")
    }
    lx <- log(x)
    classical <- .fdb_weibull_lmoment_start(x)
    sv <- resolve_start(classical$values,
                        method = classical$method,
                        positive = c("shape", "scale"),
                        sampled_parameters = "shape")
    init_shape <- sv["shape"]
    power <- if (prior == "jeffreys") n - 1 else n - 2
    logpost <- function(u) {
      eta <- u
      shape <- exp(eta)
      if (!is.finite(shape) || shape <= 0) return(-Inf)
      log_A <- .fdb_logsumexp(shape * lx)
      power * log(shape) + (shape - 1) * sum(lx) -
        n * log_A + eta
    }
    sampler <- function(iter, warmup, thin, chains, control) {
      raw <- .fdb_slice_chains(logpost, log(init_shape), iter, warmup,
                               thin, chains, control)
      out <- lapply(raw, function(eta) {
        shape <- exp(eta)
        log_A <- vapply(shape, function(k) .fdb_logsumexp(k * lx),
                        numeric(1))
        log_g <- log(.fdb_rgamma_positive(
          length(shape), shape = n, rate = 1
        ))
        scale <- .fdb_exp_positive((log_A - log_g) / shape)
        cbind(shape = shape, scale = scale)
      })
      list(chains = out, acceptance = NULL, independent = FALSE,
           engine = "marginal slice + exact conditional Gamma transform",
           initialization = initialization_state$value)
    }
    loglik <- function(theta) {
      shape <- theta["shape"]
      log_scale <- log(theta["scale"])
      log_ratio <- lx - log_scale
      powered_log_ratio <- shape * log_ratio
      answer <- rep(-Inf, length(x))
      ordinary <- is.finite(powered_log_ratio) &
        is.finite(log_ratio) &
        powered_log_ratio <= log(.Machine$double.xmax)
      answer[ordinary] <- log(shape) - log_scale +
        powered_log_ratio[ordinary] - log_ratio[ordinary] -
        exp(powered_log_ratio[ordinary])
      answer[is.nan(answer) | answer == Inf] <- -Inf
      answer
    }
    rng <- function(theta, size) stats::rweibull(size, theta["shape"],
                                                 theta["scale"])
    return(list(
      sampler = sampler, parameters = c("shape", "scale"),
      prior_label = if (prior == "jeffreys") "joint Jeffreys" else "reference",
      prior_kernel = if (prior == "jeffreys") "1/scale" else "1/(shape*scale)",
      propriety = "proper: n >= 2, positive nonconstant sample",
      loglik = loglik, rng = rng, fixed = fixed
    ))
  }

  .fdb_stop("Internal error: model '%s' has no implementation.", model)
}

.fdb_build_custom <- function(x, density, prior, start, fixed, dots, control) {
  if (!is.function(prior)) {
    .fdb_stop("A user-defined density requires a user-defined prior function.")
  }
  if (is.null(start) || !is.numeric(start) || is.null(names(start)) ||
      any(!nzchar(names(start))) || anyDuplicated(names(start)) ||
      any(!is.finite(start))) {
    .fdb_stop(
      "A user-defined density requires 'start' as a finite, uniquely named numeric vector."
    )
  }
  if (is.null(fixed)) fixed <- list()
  if (!is.list(fixed) || (length(fixed) &&
      (is.null(names(fixed)) || any(!nzchar(names(fixed))) ||
       anyDuplicated(names(fixed))))) {
    .fdb_stop("'fixed' must be NULL or a uniquely named list.")
  }
  overlap <- intersect(names(start), names(fixed))
  if (length(overlap)) {
    .fdb_stop("Parameters cannot be both unknown and fixed: %s.",
              paste(overlap, collapse = ", "))
  }
  bounds <- .fdb_custom_bounds(start, control$lower, control$upper)
  density_log_mode <- .fdb_resolve_log_mode(
    density, control$density_is_log, "density"
  )
  prior_log_mode <- .fdb_resolve_log_mode(
    prior, control$prior_is_log, "prior"
  )
  prior_formals <- names(formals(prior))
  prior_vector_style <- switch(
    control$prior_style,
    "vector" = TRUE,
    "scalar" = FALSE,
    "auto" = length(prior_formals) > 0L &&
      prior_formals[1L] %in% c("theta", "par", "parameters")
  )

  if (density_log_mode) {
    logged <- .fdb_call_density_pointwise(
      density, x, start, fixed, dots, log_mode = TRUE
    )
    ordinary <- .fdb_call_density_pointwise(
      density, x, start, fixed, dots, log_mode = FALSE
    )
    if (!.fdb_log_contract_equal(logged, ordinary)) {
      .fdb_stop(
        paste0(
          "The custom density does not satisfy its log-density contract at ",
          "'start': values from log=TRUE differ from log(values from log=FALSE)."
        )
      )
    }
  }
  if (prior_log_mode) {
    logged <- .fdb_call_prior(
      prior, start, log_mode = TRUE, vector_style = prior_vector_style
    )
    ordinary <- .fdb_call_prior(
      prior, start, log_mode = FALSE, vector_style = prior_vector_style
    )
    if (!.fdb_log_contract_equal(logged, ordinary)) {
      .fdb_stop(
        paste0(
          "The custom prior does not satisfy its log-density contract at ",
          "'start': log=TRUE differs from log(value from log=FALSE)."
        )
      )
    }
  }
  init <- .fdb_to_unconstrained(start, bounds$lower, bounds$upper)
  initialization <- list(
    source = "user-supplied",
    method = "user-specified starting values for a custom density",
    automatic = NULL,
    center = start,
    supplied = start,
    sampled_parameters = names(start)
  )
  logpost <- function(u) {
    tr <- .fdb_from_unconstrained(u, bounds$lower, bounds$upper,
                                  jacobian = TRUE)
    ll <- .fdb_call_density(
      density, x, tr$theta, fixed, dots, log_mode = density_log_mode
    )
    if (!is.finite(ll)) return(-Inf)
    lp <- .fdb_call_prior(
      prior, tr$theta, log_mode = prior_log_mode,
      vector_style = prior_vector_style
    )
    if (!is.finite(lp)) return(-Inf)
    ll + lp + tr$log_jacobian
  }
  transform <- function(u) {
    .fdb_from_unconstrained(u, bounds$lower, bounds$upper,
                            jacobian = FALSE)
  }
  sampler <- function(iter, warmup, thin, chains, control) {
    ans <- .fdb_amwg(logpost, init, iter, warmup, thin, chains, control)
    list(
      chains = .fdb_transform_matrix(ans$chains, transform, names(start)),
      acceptance = ans$acceptance,
      acceptance_warmup = ans$acceptance_warmup,
      acceptance_all = ans$acceptance_all,
      independent = FALSE,
      engine = "generic adaptive Metropolis",
      initialization = initialization
    )
  }
  loglik <- function(theta) {
    .fdb_call_density_pointwise(
      density, x, theta, fixed, dots, log_mode = density_log_mode
    )
  }
  rng <- NULL
  if (is.function(control$rng)) {
    user_rng <- control$rng
    rng <- function(theta, size) {
      ans <- do.call(user_rng,
                     c(list(size), as.list(theta), fixed, dots))
      if (!is.numeric(ans) || length(ans) != size || anyNA(ans) ||
          any(is.nan(ans)) || any(!is.finite(ans))) {
        .fdb_stop(
          paste0(
            "control$rng must return exactly 'size' finite numeric values ",
            "without NA or NaN."
          )
        )
      }
      if (is.function(control$rng_validator)) {
        valid <- control$rng_validator(ans)
        if (!is.logical(valid) || !length(valid) || anyNA(valid) ||
            !(length(valid) %in% c(1L, size)) || !all(valid)) {
          .fdb_stop(
            paste0(
              "control$rng_validator must return TRUE, or one TRUE value ",
              "per simulated observation."
            )
          )
        }
      }
      ans
    }
  }
  list(
    sampler = sampler,
    parameters = names(start),
    prior_label = "user-defined",
    prior_kernel = "user-supplied function",
    propriety = "not certified: user responsibility",
    loglik = loglik,
    rng = rng,
    fixed = fixed,
    bounds = bounds,
    custom_contract = list(
      density_is_log = density_log_mode,
      prior_is_log = prior_log_mode,
      prior_style = if (prior_vector_style) "vector" else "scalar",
      rng_validated = is.function(control$rng_validator)
    )
  )
}

#' Objective Bayesian fitting of parametric distributions
#'
#' @param x Numeric vector of observations.
#' @param distr A recognized distribution name (case-insensitive) or a density
#'   function whose first argument is the observation vector. Version 0.2.0
#'   includes 19 built-in distributions, including Gumbel, Frechet, Lomax,
#'   Nakagami-m, Exponential-Logarithmic, and Rician.
#' @param prior An objective-prior name or, for a custom density, a prior
#'   function.  A custom prior may accept named scalar parameters or a named
#'   vector as its first argument, and may optionally accept log = TRUE.
#' @param start Optional named starting values. Required for custom densities.
#' @param fixed Optional named list of fixed parameters. Currently used for
#'   negative-binomial size and Student-t degrees of freedom.
#' @param iter Total iterations per chain, including warmup.
#' @param warmup Warmup iterations per chain.
#' @param thin Positive thinning interval.
#' @param chains Number of chains.
#' @param seed Optional reproducibility seed passed to `set.seed()` before
#'   posterior simulation.
#' @param na.action Either "fail" or "omit".
#' @param control Named list of computational controls.
#' @param ... Extra fixed arguments forwarded to a custom density/RNG.
#'
#' @return An object of class "fitdistrBayes".
#' @export
fitdistrBayes <- function(x, distr, prior, start = NULL, fixed = NULL,
                          iter = 4000L, warmup = floor(iter / 2),
                          thin = 1L, chains = 4L, seed = NULL,
                          na.action = c("fail", "omit"),
                          control = list(), ...) {
  call <- match.call()
  na.action <- match.arg(na.action)
  iter <- .fdb_scalar_count(iter, "iter", lower = 20L)
  warmup <- .fdb_scalar_count(warmup, "warmup", lower = 0L)
  thin <- .fdb_scalar_count(thin, "thin", lower = 1L)
  chains <- .fdb_scalar_count(chains, "chains", lower = 1L)
  if (warmup >= iter) .fdb_stop("'warmup' must be smaller than 'iter'.")
  n_save <- floor((iter - warmup) / thin)
  if (n_save < 20L) {
    .fdb_stop("The settings retain only %d draws per chain; at least 20 are required.",
              n_save)
  }

  control <- .fdb_merge_control(control)
  validated <- .fdb_validate_x(x, na.action)
  x_used <- validated$x
  dots <- list(...)

  custom <- is.function(distr)
  if (!custom && (!is.character(distr) || length(distr) != 1L ||
                  is.na(distr) || !nzchar(trimws(distr)))) {
    .fdb_stop("'distr' must be one recognized character string or a density function.")
  }

  if (custom) {
    model <- "user-defined"
    built <- .fdb_build_custom(x_used, distr, prior, start, fixed,
                               dots, control)
    prior_key <- "user-defined"
    .fdb_warn(
      "Posterior propriety for a user-defined density/prior cannot be certified automatically."
    )
  } else {
    if (!is.character(prior) || length(prior) != 1L ||
        is.na(prior) || !nzchar(trimws(prior))) {
      .fdb_stop("Built-in distributions require 'prior' as one character string.")
    }
    model <- .fdb_model_name(distr)
    prior_key <- .fdb_prior_name(prior)
    built <- .fdb_build_builtin(x_used, model, prior_key, fixed, start)
  }

  sampled <- .fdb_with_seed(
    seed,
    built$sampler(iter, warmup, thin, chains, control)
  )
  if (!is.null(start) &&
      identical(sampled$initialization$source, "not applicable")) {
    .fdb_warn(
      "'start' was ignored because this model-prior route uses independent exact posterior simulation."
    )
  }
  chain_list <- sampled$chains
  if (length(chain_list) != chains ||
      any(vapply(chain_list, nrow, integer(1)) != n_save) ||
      any(vapply(chain_list, function(z) any(!is.finite(z)), logical(1)))) {
    .fdb_stop("Internal sampling error: invalid posterior chains were produced.")
  }
  posterior_summary <- .fdb_summarize(
    chain_list, independent = sampled$independent
  )
  moment_status <- .fdb_moment_status(
    model, prior_key, built$parameters, length(x_used), built$fixed, x_used
  )
  moment_index <- match(posterior_summary$parameter,
                        moment_status$parameter)
  posterior_summary$mean_exists <- moment_status$mean_exists[moment_index]
  posterior_summary$variance_exists <-
    moment_status$variance_exists[moment_index]
  posterior_summary$moment_note <- moment_status$note[moment_index]
  certified_mean <- posterior_summary$mean_exists %in% TRUE
  certified_variance <- posterior_summary$variance_exists %in% TRUE
  posterior_summary$mean[!certified_mean] <- NA_real_
  posterior_summary$sd[!certified_variance] <- NA_real_
  posterior_summary$mcse_mean[!(certified_mean & certified_variance)] <-
    NA_real_
  draws <- .fdb_long_draws(chain_list)

  if (sampled$independent) {
    converged <- TRUE
    diagnostic_messages <- "Independent posterior simulation: MCMC convergence is not applicable."
  } else {
    good_rhat <- is.finite(posterior_summary$rhat) &
      posterior_summary$rhat <= control$rhat_threshold
    good_ess <- is.finite(posterior_summary$ess_bulk) &
      posterior_summary$ess_bulk >= control$ess_threshold
    good_tail_ess <- is.finite(posterior_summary$ess_tail) &
      posterior_summary$ess_tail >= control$ess_threshold
    converged <- all(good_rhat & good_ess & good_tail_ess)
    diagnostic_messages <- character()
    if (any(!good_rhat)) {
      diagnostic_messages <- c(
        diagnostic_messages,
        sprintf("R-hat above %.3f for: %s.",
                control$rhat_threshold,
                paste(posterior_summary$parameter[!good_rhat],
                      collapse = ", "))
      )
    }
    if (any(!good_ess)) {
      diagnostic_messages <- c(
        diagnostic_messages,
        sprintf("Bulk ESS below %d for: %s.",
                control$ess_threshold,
                paste(posterior_summary$parameter[!good_ess],
                      collapse = ", "))
      )
    }
    if (any(!good_tail_ess)) {
      diagnostic_messages <- c(
        diagnostic_messages,
        sprintf("Tail ESS below %d for: %s.",
                control$ess_threshold,
                paste(posterior_summary$parameter[!good_tail_ess],
                      collapse = ", "))
      )
    }
    if (!length(diagnostic_messages)) {
      diagnostic_messages <- "Rank-normalized split R-hat and ESS thresholds were satisfied."
    }
  }

  acceptance <- sampled$acceptance
  if (!is.null(acceptance)) {
    warmup_acceptance <- sampled$acceptance_warmup
    all_acceptance <- sampled$acceptance_all
    acceptance_summary <- data.frame(
      chain = seq_len(nrow(acceptance)),
      overall = rowMeans(acceptance),
      warmup_overall = if (is.null(warmup_acceptance)) NA_real_ else
        rowMeans(warmup_acceptance),
      all_iterations_overall = if (is.null(all_acceptance)) NA_real_ else
        rowMeans(all_acceptance),
      acceptance,
      check.names = FALSE,
      row.names = NULL
    )
  } else {
    acceptance_summary <- NULL
  }

  diagnostics <- list(
    converged = converged,
    exact_or_independent = sampled$independent,
    rhat_threshold = control$rhat_threshold,
    ess_threshold = control$ess_threshold,
    max_rhat = max(posterior_summary$rhat, na.rm = TRUE),
    min_ess_bulk = min(posterior_summary$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(posterior_summary$ess_tail, na.rm = TRUE),
    acceptance = acceptance_summary,
    messages = diagnostic_messages,
    definition = paste0(
      "rank-normalized split/folded R-hat; mean/bulk/tail ESS; MCSE(mean); ",
      "acceptance rates labeled 'overall' use post-warmup iterations"
    )
  )

  answer <- list(
    call = call,
    model = list(
      name = model,
      distribution_input = if (custom) "<function>" else distr,
      parameters = built$parameters,
      fixed = built$fixed,
      n = length(x_used),
      support_checked = !custom,
      custom_contract = if (custom) built$custom_contract else NULL
    ),
    prior = list(
      input = if (is.function(prior)) "<function>" else prior,
      key = prior_key,
      label = built$prior_label,
      kernel = built$prior_kernel,
      posterior_propriety = built$propriety
    ),
    initialization = sampled$initialization,
    engine = list(
      algorithm = sampled$engine,
      chains = chains,
      iterations = iter,
      warmup = warmup,
      thin = thin,
      saved_per_chain = n_save,
      seed = seed
    ),
    estimates = stats::setNames(posterior_summary$median,
                                posterior_summary$parameter),
    summary = posterior_summary,
    moment_status = moment_status,
    diagnostics = diagnostics,
    capabilities = list(
      callables_stored = isTRUE(control$store_callables),
      prediction = isTRUE(control$store_callables) && !is.null(built$rng),
      pointwise_log_likelihood = isTRUE(control$store_callables)
    ),
    chains = chain_list,
    draws = draws,
    data = x_used,
    omitted = validated$omitted,
    control = control
  )
  if (control$store_callables) {
    answer$.loglik <- built$loglik
    answer$.rng <- built$rng
  }
  class(answer) <- "fitdistrBayes"

  if (!converged && isTRUE(control$warn_convergence)) {
    .fdb_warn(
      "The convergence targets were not all met. Inspect fit$diagnostics and plot(fit, type = \"trace\"); increase 'iter' if needed."
    )
  }
  answer
}

# Methods ----------------------------------------------------------------------

#' @export
print.fitdistrBayes <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("\nObjective Bayesian distribution fit\n")
  cat("-----------------------------------\n")
  cat("Model:       ", x$model$name, "\n", sep = "")
  if (length(x$model$fixed)) {
    fixed_text <- paste(
      sprintf("%s=%s", names(x$model$fixed),
              vapply(x$model$fixed,
                     function(z) paste(as.character(z), collapse = "/"),
                     character(1))),
      collapse = ", "
    )
    cat("Fixed:       ", fixed_text, "\n", sep = "")
  }
  cat("Prior:       ", x$prior$label, "\n", sep = "")
  cat("Engine:      ", x$engine$algorithm, "\n", sep = "")
  if (identical(x$initialization$source, "not applicable")) {
    cat("Initialization: not required (independent exact draws)\n")
  } else {
    cat("Initialization: ", x$initialization$source, "; ",
        x$initialization$method, "\n", sep = "")
    center_text <- paste(
      sprintf("%s=%.*g", names(x$initialization$center), digits,
              x$initialization$center),
      collapse = ", "
    )
    cat("Initial center:", center_text, "\n")
  }
  cat("Observations:", x$model$n, "\n")
  propriety_parts <- strsplit(
    x$prior$posterior_propriety, ": ", fixed = TRUE
  )[[1L]]
  cat("Propriety:   ", propriety_parts[1L], "\n", sep = "")
  if (length(propriety_parts) > 1L) {
    cat("Condition:    ",
        paste(propriety_parts[-1L], collapse = ": "), "\n", sep = "")
  }
  cat("Diagnostic:  ",
      if (x$diagnostics$converged) "OK" else "attention required",
      "\n\n", sep = "")
  shown <- x$summary[, c("parameter", "mean", "sd", "q2.5", "median",
                         "q97.5", "rhat", "ess_bulk")]
  print(shown, row.names = FALSE, digits = digits)
  if (any(!(x$summary$mean_exists %in% TRUE) |
          !(x$summary$variance_exists %in% TRUE))) {
    cat("\nNA posterior moments are intentionally suppressed when they are infinite or not certified.\n")
    cat("Posterior medians are used by coef() and fit$estimates.\n")
  }
  invisible(x)
}

#' @export
summary.fitdistrBayes <- function(object, ...) {
  answer <- list(
    call = object$call,
    model = object$model,
    prior = object$prior,
    initialization = object$initialization,
    engine = object$engine,
    posterior = object$summary,
    moment_status = object$moment_status,
    diagnostics = object$diagnostics,
    capabilities = object$capabilities
  )
  class(answer) <- "summary.fitdistrBayes"
  answer
}

#' @export
print.summary.fitdistrBayes <- function(x, digits = max(3L, getOption("digits") - 3L),
                                        ...) {
  cat("\nSummary of objective Bayesian distribution fit\n")
  cat("Model:", x$model$name, " | Prior:", x$prior$label,
      " | Engine:", x$engine$algorithm, "\n")
  if (identical(x$initialization$source, "not applicable")) {
    cat("Initialization: not required (independent exact draws)\n")
  } else {
    cat("Initialization:", x$initialization$source, "|",
        x$initialization$method, "\n")
    print(x$initialization$center, digits = digits)
  }
  cat("Posterior propriety:", x$prior$posterior_propriety, "\n\n")
  print(x$posterior, row.names = FALSE, digits = digits)
  cat("\nDiagnostics:\n")
  cat(paste0("- ", x$diagnostics$messages, collapse = "\n"), "\n")
  if (!is.null(x$diagnostics$acceptance)) {
    cat("\nAcceptance rates:\n")
    print(x$diagnostics$acceptance, row.names = FALSE, digits = digits)
  }
  if (any(!(x$posterior$mean_exists %in% TRUE) |
          !(x$posterior$variance_exists %in% TRUE))) {
    cat("\nMoment audit:\n")
    moment_lines <- unique(sprintf(
      "- %s: %s",
      x$posterior$parameter,
      x$posterior$moment_note
    ))
    cat(paste(moment_lines, collapse = "\n"), "\n")
    cat("Posterior medians are the default point estimates.\n")
  }
  invisible(x)
}

#' @export
coef.fitdistrBayes <- function(object, ...) {
  object$estimates
}

#' @export
confint.fitdistrBayes <- function(object, parm = object$model$parameters,
                                  level = 0.95, ...) {
  if (length(level) != 1L || !is.numeric(level) || !is.finite(level) ||
      level <= 0 || level >= 1) {
    .fdb_stop("'level' must lie strictly between 0 and 1.")
  }
  bad <- setdiff(parm, object$model$parameters)
  if (length(bad)) .fdb_stop("Unknown parameter%s: %s.",
                             if (length(bad) > 1L) "s" else "",
                             paste(bad, collapse = ", "))
  alpha <- (1 - level) / 2
  values <- object$draws[, parm, drop = FALSE]
  ans <- t(vapply(values, stats::quantile, numeric(2),
                  probs = c(alpha, 1 - alpha), names = FALSE, type = 8))
  colnames(ans) <- paste0(format(100 * c(alpha, 1 - alpha),
                                 trim = TRUE), "%")
  ans
}

#' @export
as.data.frame.fitdistrBayes <- function(x, row.names = NULL, optional = FALSE,
                                        ...) {
  x$draws
}

#' @export
plot.fitdistrBayes <- function(x,
                               type = c("trace", "density", "acf", "pairs"),
                               pars = x$model$parameters, ...) {
  type <- match.arg(type)
  bad <- setdiff(pars, x$model$parameters)
  if (length(bad)) .fdb_stop("Unknown plotting parameter%s: %s.",
                             if (length(bad) > 1L) "s" else "",
                             paste(bad, collapse = ", "))
  if (type == "pairs") {
    if (length(pars) < 2L) .fdb_stop("'pairs' requires at least two parameters.")
    graphics::pairs(x$draws[, pars, drop = FALSE], ...)
    return(invisible(x))
  }
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  panel_count <- if (type == "acf") {
    length(pars) * length(x$chains)
  } else {
    length(pars)
  }
  graphics::par(mfrow = grDevices::n2mfrow(panel_count),
                mar = c(3.2, 3.2, 2.2, 1))
  colors <- seq_along(x$chains)
  for (p in pars) {
    if (type == "trace") {
      yr <- range(vapply(x$chains, function(z) range(z[, p]),
                         numeric(2)))
      graphics::plot(x$chains[[1L]][, p], type = "l", col = colors[1L],
                     xlab = "Saved iteration", ylab = p, ylim = yr,
                     main = paste("Trace:", p), ...)
      if (length(x$chains) > 1L) {
        for (ch in 2:length(x$chains)) {
          graphics::lines(x$chains[[ch]][, p], col = colors[ch])
        }
      }
    } else if (type == "density") {
      dens <- lapply(x$chains, function(z) stats::density(z[, p]))
      xr <- range(vapply(dens, function(z) range(z$x), numeric(2)))
      yr <- range(vapply(dens, function(z) range(z$y), numeric(2)))
      graphics::plot(dens[[1L]], col = colors[1L], xlim = xr, ylim = yr,
                     xlab = p, main = paste("Posterior:", p), ...)
      if (length(dens) > 1L) {
        for (ch in 2:length(dens)) {
          graphics::lines(dens[[ch]], col = colors[ch])
        }
      }
    } else {
      for (ch in seq_along(x$chains)) {
        stats::acf(
          x$chains[[ch]][, p],
          main = sprintf("ACF: %s, chain %d", p, ch),
          ...
        )
      }
    }
  }
  invisible(x)
}

#' Draw from the posterior predictive distribution
#' @export
predict.fitdistrBayes <- function(object, draws = 1000L, size = 1L,
                                  seed = NULL, ...) {
  if (is.null(object$.rng)) {
    .fdb_stop(
      "Posterior prediction is unavailable. For a custom density, supply control$rng; do not set control$store_callables = FALSE."
    )
  }
  draws <- .fdb_scalar_count(draws, "draws", lower = 1L)
  size <- .fdb_scalar_count(size, "size", lower = 1L)
  all_draws <- object$draws
  answer <- .fdb_with_seed(seed, {
    index <- sample.int(nrow(all_draws), draws, replace = draws > nrow(all_draws))
    out <- matrix(NA_real_, nrow = draws, ncol = size)
    for (i in seq_len(draws)) {
      theta <- unlist(all_draws[index[i], object$model$parameters,
                                drop = FALSE], use.names = TRUE)
      out[i, ] <- object$.rng(theta, size)
    }
    out
  })
  if (anyNA(answer) || any(is.nan(answer))) {
    .fdb_stop("Posterior prediction produced undefined numeric values.")
  }
  infinite <- is.infinite(answer)
  if (any(infinite)) {
    n_infinite <- sum(infinite)
    .fdb_stop(
      paste0(
        "%d posterior predictive value%s exceeded the floating-point range; ",
        "prediction stopped and no clipping was applied."
      ),
      n_infinite, if (n_infinite == 1L) "" else "s"
    )
  }
  colnames(answer) <- paste0("y_rep[", seq_len(size), "]")
  if (size == 1L) as.numeric(answer[, 1L]) else answer
}

#' Extract a pointwise log-likelihood matrix
#' @export
log_lik.fitdistrBayes <- function(object, draws = NULL, seed = NULL, ...) {
  if (is.null(object$.loglik)) {
    .fdb_stop("Pointwise log-likelihood storage was disabled.")
  }
  all_draws <- object$draws
  if (is.null(draws)) {
    index <- seq_len(nrow(all_draws))
  } else {
    draws <- .fdb_scalar_count(draws, "draws", lower = 1L)
    index <- .fdb_with_seed(
      seed,
      sample.int(nrow(all_draws), draws, replace = draws > nrow(all_draws))
    )
  }
  values <- vapply(index, function(i) {
    theta <- unlist(all_draws[i, object$model$parameters, drop = FALSE],
                    use.names = TRUE)
    object$.loglik(theta)
  }, numeric(object$model$n))
  # vapply() simplifies to a vector when n = 1. Rebuild the documented
  # draws-by-observations matrix explicitly before transposing.
  ans <- t(matrix(
    values, nrow = object$model$n, ncol = length(index)
  ))
  colnames(ans) <- paste0("log_lik[", seq_len(object$model$n), "]")
  ans
}

log_lik <- function(object, ...) UseMethod("log_lik")
