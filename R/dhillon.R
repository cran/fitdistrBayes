# Dhillon parametrization: S(t) = 1 / (1 + theta * t^beta), beta, theta > 0.
# This is log-logistic, with log(T) logistic(location=-log(theta)/beta,
# scale=1/beta). theta is NOT a probability or an ordinary scale parameter.

.fdb_dhillon_log_density <- function(x, th) {
  b <- unname(th["beta"]); lt <- log(unname(th["theta"]))
  out <- rep(-Inf, length(x))
  keep <- is.finite(x) & x > 0
  lx <- log(x[keep]); z <- lt + b * lx
  # log logistic density, avoiding exp(z) and cancellation in z-2*softplus(z).
  out[keep] <- log(b) - lx - abs(z) - 2 * log1p(exp(-abs(z)))
  out
}

.fdb_dhillon_survival <- function(x, th) {
  out <- rep(0, length(x)); keep <- x > 0
  out[keep] <- -.fdb_softplus(log(unname(th["theta"])) +
                              unname(th["beta"]) * log(x[keep]))
  out
}

.fdb_dhillon_log_prior <- function(th, prior) {
  lt <- log(unname(th["theta"]))
  ans <- -log(unname(th["beta"])) - lt
  if (prior == "reference-theta") {
    # log(sqrt(pi^2+3+3*lt^2)) without squaring a potentially large number.
    a <- sqrt((pi^2 + 3) / 3)
    s <- max(a, abs(lt))
    ans <- ans - (0.5 * log(3) + log(s) +
                    0.5 * log((a / s)^2 + (lt / s)^2))
  }
  ans
}

.fdb_dhillon_rng <- function(th, size) {
  .fdb_exp_positive((stats::rlogis(size) - log(unname(th["theta"]))) /
                     unname(th["beta"]))
}

.fdb_dhillon_truncated <- function(lower, th, log_survival) {
  .fdb_exp_positive(.fdb_dhillon_truncated_log(lower, th, log_survival))
}

.fdb_dhillon_truncated_log <- function(lower, th, log_survival) {
  # log_survival = log(S(lower)) - Exp(1); inverse upper tail in log space.
  log_odds <- log(-expm1(log_survival)) - log_survival
  out <- (log_odds - log(unname(th["theta"]))) / unname(th["beta"])
  if (any(!is.finite(out)) || any(out <= log(lower)))
    .fdb_stop("Dhillon latent log times were not representable above their censoring limits; no clipping was applied.")
  out
}

.fdb_dhillon_logtime_density <- function(lx, th) {
  # Density of log(T). Its parameter-free Jacobian can be omitted in Gibbs
  # updates. Keeping latent times on this scale avoids exp(large log time).
  z <- log(unname(th["theta"])) + unname(th["beta"]) * lx
  log(unname(th["beta"])) - abs(z) - 2 * log1p(exp(-abs(z)))
}

.fdb_dhillon_target <- function(u, center, lx, event, prior) {
  # Posterior density in (log(beta), log(theta)+beta*center), up to a
  # data-only constant. Do NOT exponentiate log(theta) while evaluating it:
  # floating-point limits must not silently truncate the posterior target.
  b <- exp(u[1]); lt <- u[2] - b * center
  if (!is.finite(b) || b <= 0 || !is.finite(lt)) return(-Inf)
  z <- u[2] + b * (lx - center)
  a <- abs(z[event])
  ans <- sum(event) * u[1] + sum(-a - 2 * log1p(exp(-a))) -
    sum(.fdb_softplus(z[!event]))
  if (prior == "reference-theta") {
    c0 <- sqrt((pi^2 + 3) / 3); s <- max(c0, abs(lt))
    ans <- ans - (0.5 * log(3) + log(s) +
      0.5 * log((c0 / s)^2 + (lt / s)^2))
  }
  if (is.finite(ans)) unname(ans) else -Inf
}

.fdb_dhillon_proper <- function(x, status) {
  events <- x[status == 1L]
  if (!length(events)) return(FALSE)
  min(events) < max(events) || any(x[status == 0L] > max(events))
}

.fdb_dhillon_moment <- function(x, status, order) {
  # Exact criterion under either ordered reference prior or Jeffreys.
  # A_q(z) = sum_events |z+log(x)| + sum_censored (z+log(x))_+ - q*z.
  # Zero censoring limits have survival one and contribute nothing.
  keep <- x > 0; lx <- log(x[keep]); event <- status[keep] == 1L
  N <- length(lx)
  if (order >= N) return(FALSE)
  ix <- order(lx); weights <- 1L + event[ix]
  # A minimizing knot is a weighted quantile. O(n log n), not O(n^2).
  knot <- lx[ix[which(cumsum(weights) >= N - order)[1L]]]
  terms <- c(abs(lx[event] - knot), pmax(lx[!event] - knot, 0))
  value <- sum(terms) + order * knot
  if (value == 0) return(FALSE)
  tolerance <- 64 * .Machine$double.eps * (sum(terms) + abs(order * knot))
  # Do not certify a sign determined only by floating-point cancellation.
  if (abs(value) <= tolerance) return(NA)
  value > 0
}

.fdb_dhillon_moments <- function(x, status = rep(1L, length(x))) {
  data.frame(parameter = c("beta", "theta"),
    mean_exists = c(TRUE, .fdb_dhillon_moment(x, status, 1)),
    variance_exists = c(TRUE, .fdb_dhillon_moment(x, status, 2)),
    note = c("All positive posterior moments of beta exist under the registered propriety condition.",
      paste0("The q-th theta moment exists iff q < N and min A_q > 0; N counts positive observed times, ",
             "A_q(z)=sum_events abs(z+log(x))+sum_censored pmax(z+log(x),0)-q*z. ",
             "This depends on measurement units; NA denotes an unresolved numerical boundary.")),
    stringsAsFactors = FALSE)
}

.fdb_dhillon_predictive_moments <- function() {
  list(mean_exists = FALSE, variance_exists = FALSE,
    note = paste0("All positive predictive lifetime moments are infinite under these full-support posteriors. ",
                  "The predictive mean is infinite; variance and standard deviation are undefined. ",
                  "Use predictive quantiles and survival probabilities."))
}

.fdb_dhillon_start <- function(x, status, start) {
  events <- x[status == 1L]
  event_start <- length(events) >= 2L && min(events) < max(events)
  y <- if (event_start) events else x[x > 0]
  ly <- log(y); center <- mean(ly)
  scale <- sqrt(mean((ly - center)^2))
  beta <- pi / (sqrt(3) * scale)
  automatic <- c(beta = beta, theta = exp(-beta * center))
  method <- if (event_start) "closed-form moments of log exact-event times" else
    "log-moment starting heuristic from positive observed times; refined by observed-posterior geometry"
  supplied <- NULL
  if (!is.null(start)) {
    supplied <- if (is.list(start)) unlist(start, use.names = TRUE) else start
    if (!is.numeric(supplied) || is.complex(supplied) || !length(supplied) ||
        !is.null(dim(supplied)) || any(!is.finite(supplied)))
      .fdb_stop("Dhillon 'start' must be a finite named numeric vector or list.")
    if (is.null(names(supplied))) {
      if (length(supplied) != 2L) .fdb_stop("Unnamed Dhillon 'start' must have length two: beta, theta.")
      names(supplied) <- c("beta", "theta")
    }
    if (any(!nzchar(names(supplied))) || anyDuplicated(names(supplied)) ||
        any(!names(supplied) %in% c("beta", "theta")) || any(supplied <= 0))
      .fdb_stop("Dhillon starting parameters must be uniquely named beta/theta and strictly positive.")
  }
  values <- automatic
  if (!is.null(supplied)) values[names(supplied)] <- supplied
  if (any(!is.finite(values)) || any(values <= 0))
    .fdb_stop("Dhillon automatic starting values exceed numerical range; supply positive 'start' or rescale the time unit.")
  list(values = values, automatic = automatic, method = method, supplied = supplied)
}

.fdb_build_dhillon <- function(x, prior, fixed, start = NULL,
                               status = rep(1L, length(x))) {
  fixed <- .fdb_check_fixed(fixed, character())
  if (prior == "mdi")
    .fdb_stop("The MDI prior for the Dhillon model is disabled: its posterior is improper for every finite complete or right-censored sample.")
  if (prior == "reference") prior <- "reference-beta"
  if (!prior %in% c("jeffreys", "reference-beta", "reference-theta"))
    .fdb_stop("Dhillon supports jeffreys, reference-beta (alias reference), and reference-theta.")
  if (!length(x) || any(!is.finite(x)) || any(x < 0) || any(x[status == 1L] <= 0))
    .fdb_stop("Dhillon requires strictly positive exact observations and nonnegative censoring limits.")
  if (!.fdb_dhillon_proper(x, status))
    .fdb_stop(paste0("Dhillon posterior is improper: require at least one event and either two distinct ",
                     "event times or a censoring limit strictly above the common event time."))
  classical <- .fdb_dhillon_start(x, status, start)
  sv <- classical$values
  if (!is.null(start)) attr(sv, "supplied") <- TRUE
  loglik <- function(theta) .fdb_dhillon_log_density(x, theta)
  sampler <- function(iter, warmup, thin, chains, control) {
    fields <- intersect(names(control), names(.cens_control(list())))
    cc <- .cens_control(control[fields])
    ans <- .cens_mcmc("dhillon", prior, x, status, fixed, sv,
      iter, warmup, thin, chains, cc, "direct",
      function(th) .fdb_dhillon_log_prior(th, prior))
    ans$engine <- "preconditioned slice sampling in centered log coordinates"
    ans$initialization$method <- paste(classical$method,
      "observed-posterior geometry and dispersed starts", sep = "; ")
    ans$initialization$automatic <- classical$automatic
    ans$initialization$supplied <- classical$supplied
    ans
  }
  list(sampler = sampler, parameters = c("beta", "theta"), fixed = fixed,
    prior_label = if (prior == "jeffreys") "Jeffreys" else
      if (prior == "reference-beta") "ordered reference (beta interest; equal to Jeffreys)" else
        "ordered reference (theta interest)",
    prior_kernel = if (prior == "reference-theta")
      "1/(beta*theta*sqrt(pi^2+3+3*log(theta)^2))" else "1/(beta*theta)",
    propriety = paste0("proper by the exact Dhillon observed-data condition: at least one event and ",
      "either distinct event times or a censoring limit above the common event time"),
    moment_status = .fdb_dhillon_moments(x, status),
    loglik = loglik, rng = .fdb_dhillon_rng)
}

.fdb_dhillon_loo <- function(x, event) {
  n <- length(x)
  if (n < 2L) return(rep(FALSE, n))
  without_extreme <- function(z, fun, identity) {
    left <- fun(z); right <- rev(fun(rev(z)))
    before <- c(identity, left[-n]); after <- c(right[-1L], identity)
    if (identical(fun, cummin)) pmin(before, after) else pmax(before, after)
  }
  emin <- without_extreme(ifelse(event, x, Inf), cummin, Inf)
  emax <- without_extreme(ifelse(event, x, -Inf), cummax, -Inf)
  cmax <- without_extreme(ifelse(!event, x, -Inf), cummax, -Inf)
  (sum(event) - event > 0L) & (emin < emax | cmax > emax)
}
