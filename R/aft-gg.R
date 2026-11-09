# Prentice generalized gamma kernels and the ordered reference measure.
# Adapted from the authors' GG AFT manuscript; see inst/doc/gg-aft-theory.txt.
.aft_logcosh <- function(z) abs(z) + log1p(exp(-2 * abs(z))) - log(2)

.aft_logomega <- function(q) {
  s <- q * q
  if (!is.finite(s)) return(-Inf)
  if (s == 0) return(-0.5 * log(24))
  if (s < 0.0225) {
    t <- 0.5 + s * (1/6 + s^2 * (-1/30 + s^2 * (1/42 + s^2 * (-1/30 + s^2 * 5/66))))
    r <- 1/3 + s^2 * (-1/15 + s^2 * (1/21 + s^2 * (-1/15 + s^2 * 5/33)))
    return(0.5 * log((r + t*t - t)/(2 + s*(t-1))))
  }
  u <- trigamma(1 + 1/s)
  # Divide numerator and denominator by s^2 and s to avoid q^4 overflow.
  -3 * log(abs(q)) + 0.5 * (log(s) + log(u-1+(u/s)^2) - log1p(u/s))
}

.aft_logdensity <- function(w, logsigma, q) {
  if (q*q == 0) return(stats::dnorm(w, log = TRUE) - logsigma)
  s <- q*q; a <- 1/s; x <- q*w
  norm <- if (s < .02) -.5*log(2*pi)-s/12+s^3/360-s^5/1260+s^7/1680 else
    log(abs(q)) + a*log(a) - lgamma(a) - a
  h <- expm1(x) - x
  small <- is.finite(x) & abs(x) < .001
  z <- x[small]
  penalty <- h/s
  # Divide symbolically before floating-point evaluation: x^2/q^2 = w^2.
  # Computing x^2 first loses the lognormal limit when q^2 underflows.
  penalty[small] <- w[small]^2*(.5 + z*(1/6 + z*(1/24 + z*(1/120 + z/720))))
  ans <- norm - logsigma - penalty
  ans[!is.finite(w)] <- -Inf
  ans
}

.aft_logsurvival <- function(w, q) {
  if (q == 0) return(stats::pnorm(w, lower.tail = FALSE, log.p = TRUE))
  ans <- rep(NA_real_, length(w))
  ans[w == -Inf] <- 0; ans[w == Inf] <- -Inf
  use <- which(is.finite(w))
  if (!length(use)) return(ans)
  v <- w[use]
  if (abs(q) < 1e-5) {
    # The previous |w|<10 switch left small-Q tails to pgamma with a shape
    # so large that forming its argument lost the standardized deviation.
    extreme <- abs(v)>1e6
    if (any(extreme)) {
      z <- v[extreme]
      log_slope <- log(abs(expm1(q*z)))-log(abs(q))
      tail <- .aft_logdensity(z,0,q)-log_slope
      ans[use[extreme]] <- ifelse(z>0,tail,.aft_log1mexp(tail))
      use <- use[!extreme];v <- v[!extreme]
    }
    near <- abs(v)<10 | log(abs(q))+3*log(pmax(1,abs(v)))<log(.01)
    l0 <- stats::pnorm(v[near], lower.tail = FALSE, log.p = TRUE)
    z <- v[near]
    mills <- exp(stats::dnorm(z, log = TRUE)-l0)
    large <- z>50
    zz <- z[large]
    mills[large] <- zz+1/zz-2/zz^3+10/zz^5-74/zz^7+706/zz^9
    # Integrate the density expansion through q^3 against the normal tail.
    correction <- mills * (q*(z*z+2)/6-q*q*(z^5+2*z^3+6*z)/72 +
      q^3*(5*z^8-5*z^6+24*z^4+6*z*z+12)/6480)
    ans[use[near]] <- l0 + log1p(-correction)
    use <- use[!near]; v <- v[!near]
  }
  if (length(use)) {
    a <- q^-2; lg <- log(a) + q*v
    limit <- log(.Machine$double.xmax)
    argument <- exp(pmin(lg, limit))
    central <- abs(q*v)<.5
    # Avoid adding the tiny standardized displacement to a large log(a).
    argument[central] <- a+a*expm1(q*v[central])
    val <- stats::pgamma(argument, a, lower.tail = q < 0, log.p = TRUE)
    small <- lg < -30
    lgamma1 <- if (a<1e-5) -0.5772156649015329*a+(pi^2/12)*a*a else lgamma(a+1)
    lo <- a*lg[small] - lgamma1
    val[small] <- if (q < 0) lo else ifelse(lo < -log(2), log1p(-exp(lo)), log(-expm1(lo)))
    val[lg > limit] <- if (q > 0) -Inf else 0
    ans[use] <- val
  }
  ans
}

.aft_loglik <- function(theta, y, status, X) {
  p <- ncol(X); sigma <- theta[p+1L]; q <- theta[p+2L]
  if (any(!is.finite(theta)) || sigma <= 0 || !is.finite(q*q)) return(rep(-Inf, length(y)))
  w <- (y - as.vector(X %*% theta[seq_len(p)]))/sigma
  out <- numeric(length(y)); event <- status == 1L
  out[event] <- .aft_logdensity(w[event], log(sigma), q) - y[event]
  out[!event] <- .aft_logsurvival(w[!event], q)
  out
}

.aft_propriety <- function(y, status, X, tol = 1e-10) {
  event <- which(status == 1L); p <- ncol(X)
  if (length(event) <= p) return(list(ok = FALSE, reason = "The number of failures must exceed the number of regression coefficients (including the intercept)."))
  Z <- X[event, , drop = FALSE]; yy <- y[event]
  # Normalize design columns for a relative rank assessment.
  Z <- .aft_unit_columns(Z)
  qrZ <- qr(Z, tol = tol)
  if (qrZ$rank < p) return(list(ok = FALSE, reason = "The design matrix of the observed failures must have full column rank."))
  yc <- yy - mean(yy)
  res <- qr.resid(qrZ, yc)
  scale <- max(abs(yc))
  ok <- scale > 0 && sqrt(sum((res/scale)^2)) > tol * max(1, sqrt(sum((yc/scale)^2)))
  list(ok = ok, reason = if (ok) "Failure design has full rank, failures exceed coefficients, and log failure times are not exactly interpolated." else
    "Log failure times are exactly (or numerically indistinguishably) interpolated by the failure design.")
}

.aft_prior <- function(prior) {
  if (!is.character(prior) || length(prior) != 1L || is.na(prior)) stop("Specify one AFT prior name.")
  prior <- tolower(prior)
  if (prior %in% c("jeffreys", "joint-jeffreys")) stop("The joint Jeffreys GG AFT posterior is improper; no sampling was performed.")
  if (!prior %in% c("reference", "simple")) stop("AFT supports only the objective 'reference' and 'simple' priors; no proper-prior hyperparameters are exposed.")
  list(key = prior, label = prior)
}

.aft_physical <- function(v, p) c(v[seq_len(p)], exp(v[p+1L]-.aft_logcosh(v[p+2L])), sinh(v[p+2L]))

.aft_logtarget <- function(v, y, status, X, prior) {
  p <- ncol(X); th <- .aft_physical(v, p)
  ll <- .aft_loglik(th, y, status, X)
  if (any(!is.finite(ll))) return(-Inf)
  q <- th[p+2L]; lc <- .aft_logcosh(v[p+2L])
  lp <- switch(prior$key,
    reference = .aft_logomega(q) + lc,
    simple = -log1p(q*q) + lc)
  value <- sum(ll) + lp
  if (is.finite(value)) value else -Inf
}

.aft_loo_guard <- function(object) {
  n <- length(object$data)
  model <- object$model$aft_family
  check <- function(i) {
    if (is.null(model)) .aft_propriety(log(object$data[-i]),
      object$status[-i],object$design[-i,,drop=FALSE]) else
      .aft_standard_propriety(log(object$data[-i]),object$status[-i],
        object$design[-i,,drop=FALSE],model,object$prior)
  }
  ok <- vapply(seq_len(n), function(i) check(i)$ok, TRUE)
  list(certified = ok, reason = if (all(ok)) "Every AFT training sample satisfies the sufficient failure-design propriety conditions." else
    "At least one AFT training posterior is not certified by the sufficient failure-design conditions; PSIS-LOO is withheld.")
}
