# Separate censored-data interface; sampling does not call the complete-data fit.
# Censoring is independent/non-informative; status 0 means T > x (strictly,
# including for discrete distributions). Priors are those of the complete model.

.cens_control <- function(control) {
  defaults <- list(rhat_threshold = 1.01, ess_threshold = 400,
    slice_width = 1, slice_steps = 40L, max_shrink = 1000L,
    init_jitter = 2, max_init_tries = 100L, optimize_start = TRUE,
    warn_convergence = TRUE, store_callables = TRUE,
    entropy_tol = 1e-10, entropy_exact_limit = 2000L)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(!nzchar(names(control))) ||
       anyDuplicated(names(control))))) stop("'control' must be a uniquely named list.")
  bad <- setdiff(names(control), names(defaults))
  if (length(bad)) stop("Unknown control option: ", paste(bad, collapse = ", "))
  defaults[names(control)] <- control
  for (nm in c("rhat_threshold", "ess_threshold", "slice_width", "init_jitter", "entropy_tol")) {
    z <- defaults[[nm]]
    if (!is.numeric(z) || length(z) != 1L || !is.finite(z) || z <= 0)
      stop("control$", nm, " must be a positive finite scalar.")
  }
  if (defaults$rhat_threshold < 1 || defaults$entropy_tol >= 0.1)
    stop("Invalid R-hat threshold or entropy tolerance.")
  for (nm in c("slice_steps", "max_shrink", "max_init_tries", "entropy_exact_limit"))
    defaults[[nm]] <- .fdb_scalar_count(defaults[[nm]], nm)
  for (nm in c("optimize_start", "warn_convergence", "store_callables"))
    if (!is.logical(defaults[[nm]]) || length(defaults[[nm]]) != 1L || is.na(defaults[[nm]]))
      stop("control$", nm, " must be TRUE or FALSE.")
  defaults
}

.cens_models <- function() {
  z <- fitdistrBayes_routes()
  names(z)[names(z) == "engine"] <- "complete_data_engine"
  names(z)[names(z) == "posterior_condition"] <- "exact_subset_sufficient_condition"
  z$exact_subset_sufficient_condition <- paste0(z$exact_subset_sufficient_condition,
    "; apply n to the exact-event subset, not to total/imputed n")
  ix <- z$model == "t" & z$prior == "independence-jeffreys"
  z$exact_subset_sufficient_condition[ix] <- "at least 2 pairwise distinct exact events; ties are not certified"
  z$exact_subset_sufficient_condition[z$model == "dhillon"] <- paste0(
    "exact observed-data condition: at least one event and either distinct event times ",
    "or a censoring limit strictly above the common event time")
  z$censored_methods <- "direct, augmentation"
  z$censoring <- "right; status=1 event, status=0 means T>x"
  z$prior_origin <- "complete-data sampling model"
  z
}

.cens_data <- function(x, status, model, na.action) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || !length(x))
    stop("'x' must be a nonempty numeric vector.")
  if (!(is.numeric(status) || is.logical(status)) || is.complex(status) || !is.null(dim(status)) ||
      length(status) != length(x)) stop("'status' must have the same length as x and contain 0/1.")
  missing <- is.na(x) | is.na(status)
  if (any(missing) && na.action == "fail") stop("Missing x/status values; use na.action='omit' to omit pairs.")
  omitted <- which(missing)
  x <- x[!missing]; status <- status[!missing]
  if (!length(x) || any(!is.finite(x)) || any(!is.finite(status)) ||
      any(!status %in% c(0, 1))) stop("Data must be finite, with status exactly 0 or 1.")
  discrete <- model %in% c("geometric", "Poisson", "negative binomial")
  positive <- model %in% c("gamma", "weibull", "lognormal", "chi-squared", "frechet",
                          "lomax", "nakagami", "rician", "weighted lindley", "dhillon",
                          "exponential-logarithmic", "exponential")
  if (discrete && (any(x < 0) || any(x != floor(x))))
    stop("Discrete observations/censoring limits must be nonnegative integers; status=0 means T>x.")
  if (positive && any(x < 0)) stop("This model requires nonnegative observed times.")
  if (positive && model != "exponential" && any(x[status == 1] <= 0))
    stop("Exact observations must be strictly positive for this model.")
  if (model == "beta" && (any(x < 0 | x >= 1) || any(x[status == 1] <= 0)))
    stop("Beta events must lie in (0,1), and censoring limits in [0,1).")
  list(x = x, status = as.integer(status), omitted = omitted)
}

# log(1-exp(-exp(a))), without underflow in the outer log.
.cens_log1mexp_exp <- function(a) {
  out <- a
  middle <- is.finite(a) & a > -35 & a < 7
  out[middle] <- log(-expm1(-exp(a[middle])))
  out[a >= 7] <- 0
  out
}
.cens_logadd <- function(a, b) {
  z <- pmax(a, b)
  ans <- z + log1p(exp(-abs(a - b)))
  ans[is.infinite(z) & z < 0] <- -Inf
  ans
}
# log(-log(1-exp(a))), a<=0; needed for EL survival and inverse tails.
.cens_logneglog1m <- function(a) {
  out <- a
  use <- a > -35
  out[use] <- log(-log(-expm1(a[use])))
  out
}

.cens_log_i0_scaled <- function(log_argument) {
  out <- rep(NA_real_, length(log_argument))
  ordinary <- !is.na(log_argument) & log_argument <= log(1000)
  out[ordinary] <- log(besselI(exp(log_argument[ordinary]), 0, expon.scaled = TRUE))
  large <- !is.na(log_argument) & log_argument > log(1000)
  r <- exp(-log_argument[large])
  # DLMF 10.40.1, nu=0: positive asymptotic coefficients. At x>1000
  # the omitted sixth term is below double-precision relative accuracy.
  correction <- r/8 + 9*r^2/128 + 225*r^3/3072 +
    11025*r^4/98304 + 893025*r^5/3932160
  out[large] <- -0.5 * (log(2*pi) + log_argument[large]) + log1p(correction)
  out
}

.cens_density <- function(model, x, th, fixed) {
  switch(model,
    dhillon = .fdb_dhillon_log_density(x, th),
    beta = stats::dbeta(x, th["shape1"], th["shape2"], log = TRUE),
    cauchy = stats::dcauchy(x, th["location"], th["scale"], log = TRUE),
    `chi-squared` = stats::dchisq(x, th["df"], log = TRUE),
    exponential = stats::dexp(x, th["rate"], log = TRUE),
    gamma = stats::dgamma(x, th["shape"], rate = th["rate"], log = TRUE),
    geometric = stats::dgeom(x, th["prob"], log = TRUE),
    lognormal = stats::dlnorm(x, th["meanlog"], th["sdlog"], log = TRUE),
    logistic = stats::dlogis(x, th["location"], th["scale"], log = TRUE),
    `negative binomial` = stats::dnbinom(x, size = fixed$size, mu = th["mu"], log = TRUE),
    normal = stats::dnorm(x, th["mean"], th["sd"], log = TRUE),
    Poisson = stats::dpois(x, th["lambda"], log = TRUE),
    t = stats::dt((x - th["location"]) / th["scale"],
                 if (is.null(fixed$df)) th["df"] else fixed$df, log = TRUE) - log(th["scale"]),
    weibull = stats::dweibull(x, th["shape"], th["scale"], log = TRUE),
    frechet = {
      a <- th["shape"]; b <- th["scale"]
      log(a) + log(b) - (a + 1) * log(x) - exp(log(b) - a * log(x))
    },
    gumbel = {
      z <- (x - th["location"]) / th["scale"]
      -log(th["scale"]) - z - exp(-z)
    },
    lomax = log(th["shape"]) - log(th["scale"]) -
      (th["shape"] + 1) * .fdb_softplus(log(x) - log(th["scale"])),
    nakagami = stats::dgamma(x^2, th["shape"], rate = th["shape"] / th["spread"],
                            log = TRUE) + log(2) + log(x),
    `exponential-logarithmic` = {
      lp <- log(th["theta"]); lq <- log(-expm1(lp))
      log(th["rate"]) + lq - log(-lp) - th["rate"] * x -
        log(-expm1(lq - th["rate"] * x))
    },
    rician = {
      s <- th["scale"]; a <- th["noncentrality"]
      log(x) - 2 * log(s) + .cens_log_i0_scaled(log(a) + log(x) - 2*log(s)) -
        0.5 * ((x - a) / s)^2
    },
    `weighted lindley` = {
      a <- th["phi"]; b <- th["lambda"]
      (a + 1) * log(b) - .fdb_logsumexp(log(c(a, b))) - lgamma(a) +
        (a - 1) * log(x) + log1p(x) - b * x
    })
}

.cens_survival <- function(model, x, th, fixed) {
  # All calculations use the upper tail directly, not 1 - p(...).
  ans <- switch(model,
    dhillon = .fdb_dhillon_survival(x, th),
    beta = stats::pbeta(x, th["shape1"], th["shape2"], lower.tail = FALSE, log.p = TRUE),
    cauchy = stats::pcauchy(x, th["location"], th["scale"], lower.tail = FALSE, log.p = TRUE),
    `chi-squared` = stats::pchisq(x, th["df"], lower.tail = FALSE, log.p = TRUE),
    exponential = stats::pexp(x, th["rate"], lower.tail = FALSE, log.p = TRUE),
    gamma = stats::pgamma(x, th["shape"], rate = th["rate"], lower.tail = FALSE, log.p = TRUE),
    geometric = stats::pgeom(x, th["prob"], lower.tail = FALSE, log.p = TRUE),
    lognormal = stats::plnorm(x, th["meanlog"], th["sdlog"], lower.tail = FALSE, log.p = TRUE),
    logistic = stats::plogis(x, th["location"], th["scale"], lower.tail = FALSE, log.p = TRUE),
    `negative binomial` = stats::pnbinom(x, size = fixed$size, mu = th["mu"], lower.tail = FALSE, log.p = TRUE),
    normal = stats::pnorm(x, th["mean"], th["sd"], lower.tail = FALSE, log.p = TRUE),
    Poisson = stats::ppois(x, th["lambda"], lower.tail = FALSE, log.p = TRUE),
    t = stats::pt((x - th["location"]) / th["scale"],
                 if (is.null(fixed$df)) th["df"] else fixed$df, lower.tail = FALSE, log.p = TRUE),
    weibull = stats::pweibull(x, th["shape"], th["scale"], lower.tail = FALSE, log.p = TRUE),
    frechet = .cens_log1mexp_exp(log(th["scale"]) - th["shape"] * log(pmax(x, 0))),
    gumbel = .cens_log1mexp_exp(-(x - th["location"]) / th["scale"]),
    lomax = -th["shape"] * .fdb_softplus(log(pmax(x, 0)) - log(th["scale"])),
    nakagami = stats::pgamma(pmax(x, 0)^2, th["shape"], rate = th["shape"] / th["spread"],
                            lower.tail = FALSE, log.p = TRUE),
    `exponential-logarithmic` = {
      lp <- log(th["theta"])
      .cens_logneglog1m(log(-expm1(lp)) - th["rate"] * pmax(x, 0)) - log(-lp)
    },
    rician = .cens_rician_survival(x, th),
    `weighted lindley` = {
      la <- log(th["phi"]); lb <- log(th["lambda"])
      denominator <- .fdb_logsumexp(c(la, lb))
      .cens_logadd(lb - denominator + stats::pgamma(x, th["phi"], rate = th["lambda"],
                       lower.tail = FALSE, log.p = TRUE),
                   la - denominator + stats::pgamma(x, th["phi"] + 1, rate = th["lambda"],
                       lower.tail = FALSE, log.p = TRUE))
    })
  # The boundary is exact, including for mixture survival formulas whose
  # floating-point weights otherwise need not sum to exactly one.
  if (model %in% c("beta", "chi-squared", "exponential", "gamma", "lognormal",
      "weibull", "frechet", "lomax", "nakagami", "exponential-logarithmic",
      "rician", "weighted lindley", "dhillon")) ans[x <= 0] <- 0
  if (any(is.nan(ans)) || any(ans > 1e-12)) return(rep(NaN, length(x)))
  pmin(ans, 0)
}

.cens_rician_survival <- function(x, th) {
  # R's upper noncentral chi-square tail may lose precision for large ncp.
  # In those cases integrate the scaled Rician density on standardized y.
  rho <- th["noncentrality"] / th["scale"]
  q <- pmax(x, 0) / th["scale"]
  if (!is.finite(rho) || any(!is.finite(q)))
    stop("Rician standardized parameters exceeded numerical range.")
  answer <- if (rho^2 <= 1e5) suppressWarnings(stats::pchisq(q^2, df = 2, ncp = rho^2,
                                         lower.tail = FALSE, log.p = TRUE)) else rep(NA_real_, length(q))
  fallback <- which((rho^2 >= 80 & answer < log(1e-8)) | !is.finite(answer))
  for (i in fallback) {
    if (q[i] <= 0) { answer[i] <- 0; next }
    lower_z <- q[i] - rho
    logd <- function(z) {
      y <- rho + z
      ans <- rep(-Inf, length(z)); valid <- y > 0 & is.finite(y)
      ans[valid] <- log(y[valid]) + .cens_log_i0_scaled(log(rho) + log(y[valid])) - z[valid]^2 / 2
      ans
    }
    offset <- logd(max(lower_z, 0))
    f <- function(z) exp(logd(z) - offset)
    integral <- function(a,b) stats::integrate(f, a, b, rel.tol = 1e-8, subdivisions = 300L)$value
    # Center around the noncentrality so a narrow peak is not missed when
    # rho is large; split the lower half at -8 rather than discarding it.
    val <- tryCatch({
      if (lower_z >= 0) {
        h <- 1 / max(1, lower_z)
        stats::integrate(function(t) f(lower_z + h*t), 0, Inf,
                         rel.tol = 1e-8, subdivisions = 300L)$value * h
      } else {
        integral(0, Inf) + integral(max(lower_z, -8), 0) +
          if (lower_z < -8) integral(lower_z, -8) else 0
      }
    }, error = function(e) NA_real_)
    if (!is.finite(val) || val <= 0)
      stop("Rician survival could not be evaluated accurately at this parameter state.")
    answer[i] <- offset + log(val)
  }
  answer[x <= 0] <- 0
  answer
}

.cens_prior <- function(model, prior, th, fixed, control, rician_logq) {
  if (model == "dhillon") return(.fdb_dhillon_log_prior(th, prior))
  if (model %in% c("normal", "lognormal", "logistic", "gumbel", "cauchy") ||
      (model == "t" && !is.null(fixed$df))) {
    scale_name <- switch(model, normal = "sd", lognormal = "sdlog", "scale")
    return(-(if (prior == "jeffreys") 2 else 1) * log(th[scale_name]))
  }
  switch(model,
    beta = {
      a <- th["shape1"]; b <- th["shape2"]; ab <- trigamma(a + b)
      det <- (trigamma(a) - ab) * (trigamma(b) - ab) - ab^2
      if (!is.finite(det) || det <= 0) -Inf else 0.5 * log(det)
    },
    `chi-squared` = 0.5 * log(trigamma(th["df"] / 2)),
    exponential = (if (prior == "mdi") 1 else -1) * log(th["rate"]),
    gamma = {
      a <- th["shape"]
      h <- switch(prior,
        jeffreys = 0.5 * log(.fdb_gamma_joint_term(a)),
        `first-rule` = -log(a),
        `reference-shape` = 0.5 * (log(.fdb_gamma_joint_term(a)) - log(a)),
        `reference-rate` = 0.5 * log(trigamma(a)))
      h - log(th["rate"])
    },
    geometric = {
      p <- th["prob"]
      if (prior == "mdi") log(p) + (1 - p) / p * log1p(-p) else -log(p) - 0.5 * log1p(-p)
    },
    `negative binomial` = {
      if (prior == "mdi") -.fdb_nbinom_entropy(fixed$size, th["mu"], control) else
        -0.5 * (log(th["mu"]) + log(fixed$size + th["mu"]))
    },
    Poisson = if (prior == "mdi") -.fdb_poisson_entropy(th["lambda"], control) else -0.5 * log(th["lambda"]),
    t = -log(th["scale"]) + 0.5 * (log(th["df"]) - log(th["df"] + 3) + log(.fdb_t_B(th["df"]))),
    weibull = -log(th["scale"]) - if (prior == "reference") log(th["shape"]) else 0,
    frechet = -sum(log(th[c("shape", "scale")])),
    lomax = -log(th["scale"]) - 0.5 * log(th["shape"]) - log1p(th["shape"]) - 0.5 * log(th["shape"] + 2),
    nakagami = 0.5 * (log(.fdb_gamma_joint_term(th["shape"])) -
      if (prior == "reference") log(th["shape"]) else 0) - log(th["spread"]),
    `exponential-logarithmic` = {
      lp <- log(th["theta"]); terms <- .fdb_el_prior_terms(lp)
      lr <- log(th["rate"])
      switch(prior,
        jeffreys = -lr + 0.5 * terms["log_zeta"],
        mdi = lr + 0.5 * (log(-expm1(lp)) - lp - log(-lp)) + exp(terms["log_dilog"]) / lp,
        reference =, `reference-theta` = -lr + 0.5 * (terms["log_zeta"] + log(-lp) - terms["log_dilog"]),
        `reference-rate` = -lr + 0.5 * terms["log_information_theta"])
    },
    rician = 0.5 * rician_logq(2 * (log(th["noncentrality"]) - log(th["scale"]))) - 2 * log(th["scale"]),
    `weighted lindley` = .fdb_weighted_lindley_log_prior(th["lambda"], th["phi"], prior))
}

.cens_start <- function(built, model, exact, start) {
  env <- environment(built$loglik)
  if (exists("sv", env, inherits = FALSE)) {
    initial <- get("sv", env)
  } else {
    initial <- switch(model,
      exponential = c(rate = 1 / mean(exact)),
      geometric = c(prob = (length(exact) + 0.5) / (length(exact) + sum(exact) + 1)),
      Poisson = c(lambda = (sum(exact) + 0.5) / length(exact)),
      `negative binomial` = c(mu = mean(exact) + 0.5 / length(exact)),
      normal = c(mean = mean(exact), sd = sqrt(mean((exact - mean(exact))^2))),
      lognormal = c(meanlog = mean(log(exact)), sdlog = sqrt(mean((log(exact) - mean(log(exact)))^2))))
  }
  initial <- initial[built$parameters]
  if (!is.null(start)) {
    z <- if (is.list(start)) unlist(start, use.names = TRUE) else start
    if (!is.numeric(z) || !length(z) || any(!is.finite(z))) stop("'start' must contain finite numeric values.")
    if (is.null(names(z))) {
      if (length(z) != length(initial)) stop("Unnamed 'start' has incorrect length.")
      names(z) <- names(initial)
    }
    if (any(!nzchar(names(z))) || anyDuplicated(names(z)) || any(!names(z) %in% names(initial)))
      stop("Invalid or duplicated starting parameter names.")
    initial[names(z)] <- z
  }
  real <- names(initial) %in% c("location", "mean", "meanlog")
  probs <- names(initial) %in% c("prob", if (model != "dhillon") "theta")
  if (any(!is.finite(initial)) || any(initial[!real] <= 0) || any(initial[probs] >= 1))
    stop("Starting values lie outside parameter support.")
  initial
}

.cens_transform <- function(model, initial) {
  pars <- names(initial)
  real <- pars %in% c("location", "mean", "meanlog")
  prob <- pars %in% c("prob", if (model != "dhillon") "theta")
  dhillon_center <- if (model == "dhillon") unname(-log(initial["theta"]) / initial["beta"]) else 0
  loc <- which(real)
  loc_center <- if (length(loc)) initial[loc] else numeric()
  scale_name <- intersect(c("scale", "sd", "sdlog"), pars)
  loc_scale <- if (length(loc)) initial[scale_name] else 1
  to <- function(th) {
    u <- numeric(length(th))
    u[!real & !prob] <- log(th[!real & !prob])
    u[real] <- (th[real] - loc_center) / loc_scale
    u[prob] <- stats::qlogis(th[prob])
    if (model == "dhillon")
      u <- c(log(th["beta"]), log(th["theta"]) + th["beta"] * dhillon_center)
    if (model == "gamma") u <- c(log(th["shape"] / th["rate"]), log(th["shape"]))
    if (model == "beta") u <- c(log(th["shape1"] / th["shape2"]), log(sum(th)))
    if (model == "lomax") u <- c(log(th["shape"]), log(th["scale"] / th["shape"]))
    if (model == "weighted lindley") {
      a <- th["phi"]; b <- th["lambda"]
      u <- c(log(a * (b + a + 1) / (b * (b + a))), log(a))
      if (!is.finite(u[1L])) u[1L] <- log(a) +
        .fdb_logsumexp(c(log(b), log(a), 0)) - log(b) -
        .fdb_logsumexp(c(log(b), log(a)))
    }
    if (model == "gamma" && !is.finite(u[1L])) u[1L] <- log(th["shape"]) - log(th["rate"])
    if (model == "beta") {
      if (!is.finite(u[1L])) u[1L] <- log(th["shape1"]) - log(th["shape2"])
      if (!is.finite(u[2L])) u[2L] <- .fdb_logsumexp(log(th))
    }
    if (model == "lomax" && !is.finite(u[2L])) u[2L] <- log(th["scale"]) - log(th["shape"])
    unname(u)
  }
  from <- function(u) {
    th <- exp(u)
    th[real] <- loc_center + loc_scale * u[real]
    th[prob] <- stats::plogis(u[prob])
    names(th) <- pars
    jac <- sum(u[!real & !prob]) +
      sum(stats::plogis(u[prob], log.p = TRUE) + stats::plogis(u[prob], lower.tail = FALSE, log.p = TRUE)) +
      length(loc) * log(loc_scale)
    if (model == "dhillon") {
      b <- exp(u[1]); lt <- u[2] - b * dhillon_center
      th <- stats::setNames(.fdb_exp_positive(c(u[1], lt)), c("beta", "theta"))
      jac <- u[1] + lt
    }
    if (model == "gamma") {
      th <- c(shape = exp(u[2]), rate = exp(u[2] - u[1])); jac <- 2 * u[2] - u[1]
    }
    if (model == "beta") {
      p <- stats::plogis(u[1]); th <- c(shape1 = p * exp(u[2]), shape2 = (1 - p) * exp(u[2]))
      if (p == 0 || p == 1 || any(!is.finite(th))) {
        th <- c(shape1 = exp(stats::plogis(u[1], log.p = TRUE) + u[2]),
                shape2 = exp(stats::plogis(u[1], lower.tail = FALSE, log.p = TRUE) + u[2]))
      }
      jac <- sum(log(th))
    }
    if (model == "lomax") {
      th <- c(shape = exp(u[1]), scale = exp(u[1] + u[2])); jac <- 2 * u[1] + u[2]
    }
    if (model == "weighted lindley") {
      phi <- exp(u[2]); mu <- exp(u[1]); lambda <- .fdb_weighted_lindley_lambda_from_mean(mu, phi)
      th <- c(lambda = lambda, phi = phi)
      ld <- .fdb_logsumexp(log(c(lambda, phi)))
      lA <- .fdb_logsumexp(c(2 * ld, log(2) + log(lambda), log(phi)))
      jac <- u[1] + u[2] - (log(phi) + lA - 2 * log(lambda) - 2 * ld)
    }
    list(theta = th, jacobian = jac)
  }
  list(to = to, from = from)
}

.cens_geometry <- function(target, init, optimize_start) {
  center <- init; p <- length(init); ok <- FALSE
  objective <- function(u) { z <- target(u); if (is.finite(z)) -z else 1e100 }
  if (!is.finite(target(init))) stop("The initial censored log posterior is not finite; supply 'start'.")
  if (optimize_start) {
    method <- if (p == 1L) "BFGS" else "Nelder-Mead"
    opt <- tryCatch(stats::optim(init, objective, method = method,
                     control = list(maxit = 1500, reltol = 1e-9)), error = function(e) NULL)
    if (!is.null(opt) && is.finite(target(opt$par)) && objective(opt$par) <= objective(init)) {
      center <- opt$par; ok <- opt$convergence == 0
    }
  }
  h <- tryCatch(stats::optimHess(center, objective), error = function(e) NULL)
  B <- diag(1, p)
  if (!is.null(h) && all(is.finite(h))) {
    es <- eigen((h + t(h)) / 2, symmetric = TRUE)
    if (all(es$values > 0)) {
      # Bounds affect proposal geometry only; parameter space is NOT truncated.
      B <- es$vectors %*% diag(pmin(10, pmax(1e-4, 1 / sqrt(es$values))), p)
    }
  }
  list(center = center, basis = B, optimizer_converged = ok)
}

.cens_slice <- function(current, target, width, steps, max_shrink) {
  level <- target(current) - stats::rexp(1)
  if (!is.finite(level)) stop("Non-finite current state in slice sampling.")
  left <- current - width * stats::runif(1); right <- left + width
  j <- floor(stats::runif(1, 0, steps)); k <- steps - 1L - j
  while (j > 0 && target(left) > level) { left <- left - width; j <- j - 1L }
  while (k > 0 && target(right) > level) { right <- right + width; k <- k - 1L }
  for (attempt in seq_len(max_shrink)) {
    proposal <- stats::runif(1, left, right)
    lp <- target(proposal)
    if (is.finite(lp) && lp >= level) return(proposal)
    if (proposal < current) left <- proposal else right <- proposal
    if (left == right) break
  }
  stop("Slice interval exhausted numerical resolution; no invalid draw was returned.")
}

.cens_truncated <- function(model, lower, th, fixed) {
  if (!length(lower)) return(numeric())
  lp <- .cens_survival(model, lower, th, fixed) - stats::rexp(length(lower))
  if (any(!is.finite(lp))) stop("Truncated sampling exceeded the representable tail range.")
  z <- switch(model,
    dhillon = .fdb_dhillon_truncated(lower, th, lp),
    beta = stats::qbeta(lp, th["shape1"], th["shape2"], lower.tail = FALSE, log.p = TRUE),
    cauchy = stats::qcauchy(lp, th["location"], th["scale"], lower.tail = FALSE, log.p = TRUE),
    `chi-squared` = stats::qchisq(lp, th["df"], lower.tail = FALSE, log.p = TRUE),
    exponential = stats::qexp(lp, th["rate"], lower.tail = FALSE, log.p = TRUE),
    gamma = stats::qgamma(lp, th["shape"], rate = th["rate"], lower.tail = FALSE, log.p = TRUE),
    geometric = stats::qgeom(lp, th["prob"], lower.tail = FALSE, log.p = TRUE),
    lognormal = stats::qlnorm(lp, th["meanlog"], th["sdlog"], lower.tail = FALSE, log.p = TRUE),
    logistic = stats::qlogis(lp, th["location"], th["scale"], lower.tail = FALSE, log.p = TRUE),
    `negative binomial` = stats::qnbinom(lp, size = fixed$size, mu = th["mu"], lower.tail = FALSE, log.p = TRUE),
    normal = stats::qnorm(lp, th["mean"], th["sd"], lower.tail = FALSE, log.p = TRUE),
    Poisson = stats::qpois(lp, th["lambda"], lower.tail = FALSE, log.p = TRUE),
    t = th["location"] + th["scale"] * stats::qt(lp,
            if (is.null(fixed$df)) th["df"] else fixed$df, lower.tail = FALSE, log.p = TRUE),
    weibull = stats::qweibull(lp, th["shape"], th["scale"], lower.tail = FALSE, log.p = TRUE),
    frechet = exp((log(th["scale"]) - .cens_logneglog1m(lp)) / th["shape"]),
    gumbel = th["location"] - th["scale"] * .cens_logneglog1m(lp),
    lomax = th["scale"] * expm1(-lp / th["shape"]),
    nakagami = sqrt(stats::qgamma(lp, th["shape"], rate = th["shape"] / th["spread"],
                                  lower.tail = FALSE, log.p = TRUE)),
    `exponential-logarithmic` = {
      ltheta <- log(th["theta"])
      (log(-expm1(ltheta)) - .cens_log1mexp_exp(lp + log(-ltheta))) / th["rate"]
    },
    rician = {
      # Fast vectorized inversion in the moderate-ncp regime, checked against
      # our survival. Only problematic quantiles need a scalar root search.
      rho2 <- (th["noncentrality"] / th["scale"])^2
      ans <- if (is.finite(rho2) && rho2 < 80) suppressWarnings(th["scale"] *
        sqrt(stats::qchisq(lp, df = 2, ncp = rho2, lower.tail = FALSE, log.p = TRUE))) else rep(NA_real_, length(lp))
      error <- rep(Inf, length(lp))
      good <- is.finite(ans) & ans > lower
      if (any(good)) error[good] <- abs(.cens_survival(model, ans[good], th, fixed) - lp[good])
      bad <- which(!is.finite(error) | error > 1e-7)
      for (i in bad) {
        hi <- max(lower[i] + th["scale"], th["noncentrality"] + 8 * th["scale"])
        fun <- function(z) .cens_survival(model, z, th, fixed) - lp[i]
        for (j in seq_len(100L)) {
          if (fun(hi) <= 0) break
          hi <- 2 * hi + th["scale"]
        }
        ans[i] <- stats::uniroot(fun, c(lower[i], hi), tol = max(1e-12 * th["scale"], .Machine$double.eps))$root
      }
      ans
    },
    `weighted lindley` = {
      a <- th["phi"]; b <- th["lambda"]
      l1 <- log(b) + stats::pgamma(lower, a, rate = b, lower.tail = FALSE, log.p = TRUE)
      l2 <- log(a) + stats::pgamma(lower, a + 1, rate = b, lower.tail = FALSE, log.p = TRUE)
      component <- stats::runif(length(lower)) > exp(l1 - .cens_logadd(l1, l2))
      shape <- a + component
      ltail <- stats::pgamma(lower, shape, rate = b, lower.tail = FALSE, log.p = TRUE)
      stats::qgamma(ltail - stats::rexp(length(lower)), shape, rate = b, lower.tail = FALSE, log.p = TRUE)
    })
  if (any(!is.finite(z)) || any(z <= lower))
    stop("Truncated draws were not representable strictly above their censoring limits; no clipping was applied.")
  unname(z)
}

.cens_mcmc <- function(model, prior, x, status, fixed, initial, iter, warmup,
                       thin, chains, control, method, logprior) {
  transform <- .cens_transform(model, initial)
  event <- status == 1L; cens <- !event
  state <- new.env(parent = emptyenv()); state$data <- x
  if (model == "dhillon") state$log_data <- log(x)
  dhillon_center <- if (model == "dhillon") unname(-log(initial["theta"]) / initial["beta"]) else 0
  evaluate <- function(u, augmented = FALSE) {
    if (any(!is.finite(u))) return(-Inf)
    if (model == "dhillon") return(.fdb_dhillon_target(u, dhillon_center,
      if (augmented) state$log_data else log(x),
      if (augmented) rep(TRUE, length(x)) else event, prior))
    tr <- suppressWarnings(transform$from(u)); th <- tr$theta
    if (any(!is.finite(th)) || any(th[!names(th) %in% c("location", "mean", "meanlog")] <= 0) ||
        any(th[names(th) %in% c("prob", if (model != "dhillon") "theta")] >= 1) || !is.finite(tr$jacobian)) return(-Inf)
    priorval <- suppressWarnings(logprior(th))
    if (length(priorval) != 1L || !is.finite(priorval)) return(-Inf)
    if (augmented) {
      ll <- if (model == "dhillon") .fdb_dhillon_logtime_density(state$log_data, th) else
        suppressWarnings(.cens_density(model, state$data, th, fixed))
    } else {
      ll <- c(suppressWarnings(.cens_density(model, x[event], th, fixed)),
              suppressWarnings(.cens_survival(model, x[cens], th, fixed)))
    }
    ans <- sum(ll) + priorval + tr$jacobian
    if (is.finite(ans)) unname(ans) else -Inf
  }
  geometry <- .cens_geometry(evaluate, transform$to(initial),
                            control$optimize_start && is.null(attr(initial, "supplied")))
  target <- function(z) evaluate(as.numeric(geometry$center + geometry$basis %*% z),
                                  augmented = method == "augmentation")
  out <- vector("list", chains); p <- length(initial)
  starts <- matrix(NA_real_, chains, p, dimnames = list(NULL, names(initial)))
  for (ch in seq_len(chains)) {
    z <- rep(0, p)
    for (attempt in seq_len(control$max_init_tries)) {
      candidate <- stats::rnorm(p, sd = control$init_jitter)
      if (is.finite(evaluate(as.numeric(geometry$center + geometry$basis %*% candidate)))) {
        z <- candidate; break
      }
    }
    theta <- transform$from(as.numeric(geometry$center + geometry$basis %*% z))$theta
    starts[ch, ] <- theta
    saved <- matrix(NA_real_, floor((iter - warmup) / thin), p,
                    dimnames = list(NULL, names(initial)))
    pos <- 0L
    for (i in seq_len(iter)) {
      if (method == "augmentation") {
        theta <- transform$from(as.numeric(geometry$center + geometry$basis %*% z))$theta
        if (model == "dhillon") {
          lp <- .fdb_dhillon_survival(x[cens], theta) - stats::rexp(sum(cens))
          state$log_data[cens] <- .fdb_dhillon_truncated_log(x[cens], theta, lp)
        } else state$data[cens] <- .cens_truncated(model, x[cens], theta, fixed)
        # Exact blocked Normal-inverse-Gamma conditional, not a restarted fit.
        if (model %in% c("normal", "lognormal")) {
          y <- if (model == "normal") state$data else log(state$data)
          n <- length(y); ybar <- mean(y); q <- sum((y - ybar)^2)
          a <- if (prior == "jeffreys") 2 else 1
          s2 <- 1 / .fdb_rgamma_positive(1, (n + a - 2) / 2, rate = q / 2)
          theta <- stats::setNames(c(stats::rnorm(1, ybar, sqrt(s2 / n)), sqrt(s2)), names(initial))
          z <- as.numeric(solve(geometry$basis, transform$to(theta) - geometry$center))
        } else {
          for (j in sample.int(p)) {
            one <- function(v) { proposed <- z; proposed[j] <- v; target(proposed) }
            z[j] <- .cens_slice(z[j], one, control$slice_width,
                                control$slice_steps, control$max_shrink)
          }
        }
      } else {
        for (j in sample.int(p)) {
          one <- function(v) { proposed <- z; proposed[j] <- v; target(proposed) }
          z[j] <- .cens_slice(z[j], one, control$slice_width,
                              control$slice_steps, control$max_shrink)
        }
      }
      if (i > warmup && (i - warmup) %% thin == 0L) {
        pos <- pos + 1L
        saved[pos, ] <- transform$from(as.numeric(geometry$center + geometry$basis %*% z))$theta
      }
    }
    out[[ch]] <- saved
  }
  list(chains = out, independent = FALSE,
       engine = if (method == "augmentation") {
         if (model %in% c("normal", "lognormal")) "data augmentation + exact Normal-inverse-Gamma blocks" else
           "data augmentation + preconditioned slice-within-Gibbs"
       } else "observed-likelihood preconditioned slice sampling",
       initialization = list(source = if (is.null(attr(initial, "supplied"))) "automatic" else "user-supplied center",
         method = "complete-event classical start; observed-posterior geometry; dispersed chains",
         classical = initial, center = transform$from(geometry$center)$theta,
         chain_starts = starts, optimizer_converged = geometry$optimizer_converged))
}

.cens_exact <- function(model, prior, x, status, initial, iter, warmup, thin, chains, control) {
  m <- sum(status); A <- sum(x); nc <- sum(status == 0L)
  if (model %in% c("exponential", "geometric") && !is.finite(A))
    stop("The exact censored posterior sufficient statistic exceeded numerical range.")
  nsave <- floor((iter - warmup) / thin)
  if (model == "exponential") {
    out <- lapply(seq_len(chains), function(ch) matrix(.fdb_rgamma_positive(nsave,
      m + if (prior == "mdi") 2 else 0, rate = A), ncol = 1, dimnames = list(NULL, "rate")))
    label <- "exact censored Gamma posterior"
  } else if (model == "geometric") {
    out <- lapply(seq_len(chains), function(ch) matrix(stats::rbeta(nsave, m, A + nc + 0.5),
      ncol = 1, dimnames = list(NULL, "prob")))
    label <- "exact censored Beta posterior"
  } else {
    lx <- log(x); event_lx <- lx[status == 1]
    power <- m - if (prior == "jeffreys") 1 else 2
    target <- function(u) {
      k <- exp(u)
      if (!is.finite(k) || k <= 0) return(-Inf)
      ans <- power * u + (k - 1) * sum(event_lx) - m * .fdb_logsumexp(k * lx) + u
      if (is.finite(ans)) ans else -Inf
    }
    raw <- .fdb_slice_chains(target, log(initial["shape"]), iter, warmup, thin, chains,
      .fdb_merge_control(list(slice_width = control$slice_width, slice_steps = control$slice_steps,
                              init_jitter = control$init_jitter)))
    out <- lapply(raw, function(u) {
      k <- exp(u)
      logA <- vapply(k, function(v) .fdb_logsumexp(v * lx), numeric(1))
      sc <- .fdb_exp_positive((logA - log(.fdb_rgamma_positive(length(k), m))) / k)
      cbind(shape = k, scale = sc)
    })
    label <- "censored marginal slice + exact conditional Weibull scale"
  }
  list(chains = out, independent = model != "weibull", engine = label,
       initialization = list(source = if (model == "weibull") "classical event-subset start" else "not applicable",
                             method = label, center = initial))
}

.cens_summarize <- function(chains, independent) {
  out <- .fdb_summarize(chains, independent)
  if (!independent) for (j in seq_len(nrow(out))) {
    split <- .fdb_split_matrix(chains, j)
    scale <- .fdb_stable_sd(as.numeric(split))
    chain_moves <- vapply(chains, function(z) length(unique(z[, j])) > 1L, logical(1))
    if (!all(chain_moves) || !is.finite(scale) || scale <= 0) {
      out$rhat[j] <- Inf
      out$ess_mean[j] <- out$ess_bulk[j] <- out$ess_tail[j] <- 0
      out$mcse_mean[j] <- NA_real_
    } else {
      # Scaling protects mean ESS/MCSE from absolute variance tolerances when
      # users change units (rank-based R-hat and bulk/tail ESS are invariant).
      standardized <- (split - mean(split)) / scale
      out$ess_mean[j] <- .fdb_ess_matrix(standardized)
      out$mcse_mean[j] <- out$sd[j] / sqrt(out$ess_mean[j])
    }
  }
  out
}

fitcensBayes <- function(x, status, distr, prior = NULL, start = NULL, fixed = NULL,
                        iter = 4000L, warmup = floor(iter / 2), thin = 1L,
                        chains = 4L, seed = NULL, method = c("auto", "direct", "augmentation"),
                        na.action = c("fail", "omit"), control = list(), ...,
                        criteria = FALSE) {
  call <- match.call()
  if (length(list(...))) stop("Unused arguments: ", paste(names(list(...)), collapse = ", "))
  if (!is.character(distr) || length(distr) != 1L || is.na(distr))
    stop("This separate extension accepts the 21 registered distribution names, not custom densities yet.")
  model <- .fdb_model_name(distr)
  if (!is.character(prior) || length(prior) != 1L || is.na(prior) || !nzchar(prior))
    stop("Specify one registered prior name.")
  prior_key <- .fdb_prior_name(prior, model)
  method <- match.arg(method); na.action <- match.arg(na.action)
  control <- .cens_control(control)
  if (!identical(criteria, FALSE)) .fdb_criteria_preflight(criteria, control)
  iter <- .fdb_scalar_count(iter, "iter", lower = 20L)
  warmup <- .fdb_scalar_count(warmup, "warmup", lower = 0L)
  if (warmup >= iter) stop("'warmup' must be smaller than 'iter'.")
  thin <- .fdb_scalar_count(thin, "thin")
  chains <- .fdb_scalar_count(chains, "chains", lower = 2L)
  nsave <- floor((iter - warmup) / thin)
  if (nsave < 10L) stop("Keep at least 10 draws per chain after warmup/thinning.")
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed))) stop("Invalid seed.")
  dat <- .cens_data(x, status, model, na.action)
  x <- dat$x; status <- dat$status; exact <- x[status == 1L]
  m <- length(exact); nc <- sum(status == 0L)
  # Check priors/fixed parameters before the event-subset propriety guard.
  # The upstream validator supplies detailed errors for aliases and fixed df.
  if (!m && model != "dhillon" && !(model == "exponential" && prior_key == "mdi" && sum(x) > 0))
    stop("Posterior propriety is not certified with zero events for this route; no sampling was performed.")
  guard_data <- if (!m || (model == "exponential" && sum(exact) == 0 && sum(x) > 0)) 1 else exact
  built <- if (model == "dhillon") .fdb_build_dhillon(x, prior_key, fixed, start, status) else
    tryCatch(.fdb_build_builtin(guard_data, model, prior_key, fixed, start),
    error = function(e) stop("Censored-model validation failed. The complete-event sufficient criterion or model/prior restriction was not met: ",
                            conditionMessage(e), ". This does not by itself establish impropriety of every censored case.", call. = FALSE))
  fixed <- built$fixed
  if (model == "exponential" && (!is.finite(sum(x)) || sum(x) <= 0))
    stop("Exponential posterior requires a finite positive total observed time.")
  initial <- .cens_start(built, model, guard_data, start)
  if (!is.null(start)) attr(initial, "supplied") <- TRUE
  rq <- if (model == "rician") .fdb_rician_logq_factory() else NULL
  logprior <- function(th) .cens_prior(model, prior_key, th, fixed, control, rq)
  exact_route <- model %in% c("exponential", "weibull") ||
    (model == "geometric" && prior_key != "mdi")
  selected <- if (method == "auto") "direct" else method
  sampled <- .fdb_with_seed(seed, {
    if (nc == 0L && method == "auto") {
      up_control <- control[intersect(names(control), names(.fdb_merge_control(list())))]
      built$sampler(iter, warmup, thin, chains, .fdb_merge_control(up_control))
    } else if (exact_route && selected == "direct") {
      .cens_exact(model, prior_key, x, status, initial, iter, warmup, thin, chains, control)
    } else {
      .cens_mcmc(model, prior_key, x, status, fixed, initial, iter, warmup, thin,
                 chains, control, selected, logprior)
    }
  })
  chain_list <- sampled$chains
  if (length(chain_list) != chains || any(vapply(chain_list, nrow, integer(1)) != nsave) ||
      any(vapply(chain_list, function(z) any(!is.finite(z)), logical(1)))) stop("Invalid posterior chains; fit not returned.")
  positive_parameters <- setdiff(built$parameters, c("location", "mean", "meanlog"))
  probability_parameters <- intersect(built$parameters, c("prob", if (model != "dhillon") "theta"))
  if (any(vapply(chain_list, function(z)
    any(z[, positive_parameters, drop=FALSE] <= 0) ||
    any(z[, probability_parameters, drop=FALSE] >= 1), logical(1))))
    stop("Posterior draws reached a nonrepresentable parameter boundary; no clipping was applied.")
  posterior <- .cens_summarize(chain_list, independent = sampled$independent)
  moments <- if (model == "dhillon") built$moment_status else
    .fdb_moment_status(model, prior_key, built$parameters, length(guard_data), fixed, guard_data)
  if (nc && model != "dhillon") {
    # Dominance proves existence, never nonexistence, of moments. Unproved
    # censored moments are NA, not asserted infinite from complete-data tests.
    for (col in c("mean_exists", "variance_exists")) moments[[col]][!moments[[col]] %in% TRUE] <- NA
    moments$note <- "Finite moments transferred from exact-event subset when TRUE; NA means not certified under censoring."
    if (model == "exponential") {
      moments$mean_exists <- moments$variance_exists <- TRUE
      moments$note <- "Exact Gamma posterior: all positive rate moments are finite."
    }
    if (model == "weibull") {
      idx <- moments$parameter == "scale"
      moments[idx, c("mean_exists", "variance_exists")] <- FALSE
      moments$note[idx] <- "Under the censored Weibull conditional Gamma form, positive scale moments diverge as shape approaches zero."
    }
  }
  posterior$mean_exists <- moments$mean_exists
  posterior$variance_exists <- moments$variance_exists
  posterior$moment_note <- moments$note
  posterior$mean[!posterior$mean_exists %in% TRUE] <- NA_real_
  posterior$sd[!posterior$variance_exists %in% TRUE] <- NA_real_
  posterior$mcse_mean[!(posterior$mean_exists %in% TRUE & posterior$variance_exists %in% TRUE)] <- NA_real_
  good <- is.finite(posterior$rhat) & posterior$rhat <= control$rhat_threshold &
    is.finite(posterior$ess_bulk) & posterior$ess_bulk >= control$ess_threshold &
    is.finite(posterior$ess_tail) & posterior$ess_tail >= control$ess_threshold
  passed <- sampled$independent || all(good)
  message <- if (sampled$independent) "Independent exact simulation; MCMC convergence is not applicable." else
    if (passed) "Rank-normalized split/folded R-hat and bulk/tail ESS targets met; this is not a proof of convergence." else
      paste0("Diagnostic targets not met for: ", paste(posterior$parameter[!good], collapse = ", "), ". Increase iter and inspect chains.")
  propriety <- if (model == "exponential") "proper: exact censored Gamma posterior" else
    if (!nc || model == "dhillon") built$propriety else paste0("proper by exact-event subset dominance; m=", m, "; ", built$propriety)
  answer <- list(call = call,
    model = list(name = model, parameters = built$parameters, fixed = fixed, n = length(x)),
    prior = list(input = prior, key = prior_key, label = built$prior_label, kernel = built$prior_kernel,
                 origin = "complete-data sampling model", posterior_propriety = propriety),
    censoring = list(type = "right", events = m, censored = nc, fraction = nc / length(x),
                     assumption = "independent/non-informative; status=0 means T>x"),
    initialization = sampled$initialization,
    engine = list(algorithm = sampled$engine, requested = method, iterations = iter, warmup = warmup,
                  thin = thin, chains = chains, saved_per_chain = nsave, seed = seed),
    estimates = stats::setNames(posterior$median, posterior$parameter), summary = posterior,
    moment_status = moments,
    diagnostics = list(converged = passed, exact_or_independent = sampled$independent,
      max_rhat = max(posterior$rhat), min_ess_bulk = min(posterior$ess_bulk), min_ess_tail = min(posterior$ess_tail),
      rhat_threshold = control$rhat_threshold, ess_threshold = control$ess_threshold, messages = message),
    chains = chain_list, draws = .fdb_long_draws(chain_list), data = x, status = status,
    omitted = dat$omitted, control = control)
  if (control$store_callables) {
    answer$.loglik <- function(th) {
      out <- numeric(length(x)); event <- status == 1
      out[event] <- .cens_density(model, x[event], th, fixed)
      out[!event] <- .cens_survival(model, x[!event], th, fixed)
      out
    }
    answer$.rng <- built$rng
    answer$.survival <- function(times, th) .cens_survival(model, times, th, fixed)
    answer$.truncated <- function(lower, th) .cens_truncated(model, lower, th, fixed)
  }
  class(answer) <- "fitcensBayes"
  if (model == "dhillon") {
    answer$predictive_moments <- .fdb_dhillon_predictive_moments()
    answer$initialization$method <- paste0("Dhillon log-moment start; observed-posterior geometry; dispersed chains")
  }
  if (!identical(criteria, FALSE)) {
    answer$criteria <- .fdb_criteria_after_fit(answer, criteria)
  }
  if (!passed && control$warn_convergence) warning(message, call. = FALSE)
  answer
}

.cens_print <- function(x, digits = 4L, ...) {
  cat("\nObjective Bayesian fit with right censoring\n")
  cat("Model:", x$model$name, " | Prior:", x$prior$label, "\n")
  cat("Events:", x$censoring$events, " | Censored:", x$censoring$censored,
      sprintf(" (%.1f%%)\n", 100 * x$censoring$fraction))
  cat("Algorithm:", x$engine$algorithm, "\n")
  cat("Diagnostic:", if (x$diagnostics$converged) "targets met" else "attention required", "\n\n")
  print(x$summary[, c("parameter", "mean", "sd", "q2.5", "median", "q97.5", "rhat", "ess_bulk", "ess_tail")],
        row.names = FALSE, digits = digits)
  if (anyNA(x$summary$mean) || anyNA(x$summary$sd))
    cat("\nNA moments are infinite or not certified. coef() returns posterior medians.\n")
  invisible(x)
}

.cens_plot <- function(x, type, pars, ...) {
  if (!is.character(pars) || !length(pars) || any(!pars %in% x$model$parameters))
    stop("'pars' must name fitted parameters.")
  if (type == "pairs") {
    if (length(pars) < 2L) stop("Pairs plots require at least two parameters.")
    graphics::pairs(x$draws[, pars, drop = FALSE], ...)
    return(invisible(x))
  }
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old))
  count <- length(pars) * if (type == "acf") length(x$chains) else 1L
  graphics::par(mfrow = grDevices::n2mfrow(count), mar = c(4.2, 5, 2.5, 1.2))
  for (p in pars) {
    if (type == "trace") {
      series <- lapply(x$chains, function(z) z[, p])
      args <- list(x = seq_along(series[[1]]), y = series[[1]], type = "l", col = 1,
                   ylim = range(unlist(series)), xlab = "Saved iteration", ylab = p,
                   main = paste("Trace:", p))
      do.call(graphics::plot, utils::modifyList(args, list(...)))
      for (ch in seq.int(2L, length(series))) graphics::lines(series[[ch]], col = ch)
    } else if (type == "density") {
      dens <- lapply(x$chains, function(z) stats::density(z[, p]))
      args <- list(x = dens[[1]], col = 1, xlab = p, main = paste("Posterior:", p),
                   xlim = range(vapply(dens, function(z) range(z$x), numeric(2))),
                   ylim = range(vapply(dens, function(z) range(z$y), numeric(2))))
      do.call(graphics::plot, utils::modifyList(args, list(...)))
      for (ch in seq.int(2L, length(dens))) graphics::lines(dens[[ch]], col = ch)
    } else {
      for (ch in seq_along(x$chains)) {
        args <- list(x = x$chains[[ch]][, p], main = sprintf("ACF: %s, chain %d", p, ch))
        do.call(stats::acf, utils::modifyList(args, list(...)))
      }
    }
  }
  invisible(x)
}

.cens_predict <- function(object, type, times, draws, size, seed) {
  draws <- .fdb_scalar_count(draws, "draws"); size <- .fdb_scalar_count(size, "size")
  if (is.null(object$.rng)) stop("Prediction callbacks were not stored.")
  if (type == "survival" && (!is.numeric(times) || !length(times) || any(!is.finite(times))))
    stop("Supply finite numeric 'times' for survival prediction.")
  if (type == "impute" && !any(object$status == 0)) stop("There are no censored observations to impute.")
  if (type != "survival" && !is.null(times)) stop("'times' is only used for type='survival'.")
  .fdb_with_seed(seed, {
    idx <- sample.int(nrow(object$draws), draws, replace = draws > nrow(object$draws))
    cols <- switch(type, response = size, survival = length(times), impute = sum(object$status == 0))
    ans <- matrix(NA_real_, draws, cols)
    for (i in seq_len(draws)) {
      th <- unlist(object$draws[idx[i], object$model$parameters, drop = FALSE], use.names = TRUE)
      ans[i, ] <- switch(type, response = object$.rng(th, size),
                        survival = exp(object$.survival(times, th)),
                        impute = object$.truncated(object$data[object$status == 0], th))
    }
    if (any(!is.finite(ans))) stop("Prediction exceeded numerical range; no clipping was applied.")
    colnames(ans) <- switch(type,
      response = paste0("replicate[", seq_len(cols), "]"),
      survival = paste0("S(", times, ")"),
      impute = paste0("observation[", which(object$status == 0), "]"))
    ans
  })
}
