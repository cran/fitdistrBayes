# A fixed reparametrization of the SAME posterior, not a new prior.
# h(Q) is a standardized log-time quantile divided by sqrt(1+Q^2).
.aft_quantile_offset <- function(q, prob) {
  if (!is.finite(q) || !is.finite(q*q)) return(NA_real_)
  if (abs(q) < 1e-4) {
    z <- stats::qnorm(prob)
    return((z-q*(z*z+2)/6+q*q*(z^3+5*z)/36)/sqrt(1+q*q))
  }
  a <- q^-2
  lp <- if (q > 0) log(prob) else log1p(-prob)
  logg <- (lp+lgamma(a+1))/a
  # Gamma quantiles can underflow even though the log-time offset is finite.
  # The small-gamma-argument expansion has relative error O(exp(logg)).
  if (logg < -30)
    return(sign(q)*(lp+lgamma(a+1)-a*log(a))/sqrt(1+a))
  g <- stats::qgamma(lp, shape = a, log.p = TRUE)
  (log(g)-log(a))/(q*sqrt(1+q*q))
}

.aft_uncenter <- function(v, prob) {
  d <- length(v)
  v[1L] <- v[1L]-exp(v[d-1L])*.aft_quantile_offset(sinh(v[d]), prob)
  v
}

.aft_sample <- function(target, mode, p, iter, warmup, thin, chains, control,
                        parameters, center_prob) {
  d <- p+2L
  C <- t(chol(mode$V))
  D <- diag(sqrt(diag(mode$V)), nrow = d)
  nsave <- floor((iter-warmup)/thin)
  output <- starts <- vector("list", chains)
  movement <- evaluations <- numeric(chains)
  for (k in seq_len(chains)) {
    calls <- 0L
    counted_target <- function(v) {calls <<- calls+1L; target(v)}
    initialized <- FALSE
    for (j in seq_len(100L)) {
      v <- mode$center + as.vector(C %*% stats::rnorm(d))*(0.5+k/chains)
      if (is.finite(counted_target(v))) {initialized <- TRUE; break}
    }
    if (!initialized) stop("Could not generate dispersed finite AFT chain starts.")
    starts[[k]] <- v
    out <- matrix(NA_real_, nsave, d, dimnames = list(NULL, parameters))
    moved <- save <- 0L
    for (i in seq_len(iter)) {
      previous <- v
      # Both bases are fixed and the choice is state-independent. Each line
      # update is Neal's stepping-out/shrinkage slice kernel. Mixing bases
      # avoids imposing the mode's correlations throughout the shape tails.
      basis <- if (stats::runif(1) < .5) C else D
      for (j in sample.int(d)) {
        direction <- basis[, j]
        one <- function(s) counted_target(v+s*direction)
        step <- .cens_slice(0, one, control$slice_width,
                            control$slice_steps, control$max_shrink)
        v <- v+step*direction
      }
      if (i > warmup) {
        moved <- moved+as.integer(any(v != previous))
        if ((i-warmup) %% thin == 0L) {
          save <- save+1L
          th <- .aft_physical(.aft_uncenter(v, center_prob), p)
          if (any(!is.finite(th)) || th[p+1L] <= 0)
            stop("AFT draw is outside the representable range; no clipping was applied.")
          out[save, ] <- th
        }
      }
    }
    output[[k]] <- out
    movement[k] <- moved/(iter-warmup)
    evaluations[k] <- calls
  }
  list(chains = output, movement = movement, evaluations = evaluations, starts = starts)
}
