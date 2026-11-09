# Approximate quantile MCSE from the ESS of the empirical-CDF indicator.
# Inverting a one-standard-error probability interval avoids a density
# estimate and does not require finite parameter means or variances.
.aft_quantile_precision <- function(chains) {
  probs <- c(.025, .5, .975)
  parameters <- colnames(chains[[1L]])
  do.call(rbind, lapply(seq_along(parameters), function(j) {
    values <- unlist(lapply(chains, function(ch) ch[, j]), use.names = FALSE)
    split <- .fdb_split_matrix(chains, j)
    q <- stats::quantile(values, probs, names = FALSE, type = 8)
    ess <- mcse <- rep(NA_real_, 3L)
    moving <- all(vapply(chains, function(ch) any(ch[, j] != ch[1L, j]), TRUE))
    if (moving && !is.null(split)) for (k in seq_along(probs)) {
      ess[k] <- .fdb_ess_matrix(1*(split <= q[k]))
      # A few effective tail observations, or a stuck chain, cannot certify
      # precision. min(p,1-p)*ESS>=5 also keeps both CDF bounds inside (0,1).
      if (is.finite(ess[k]) && ess[k]*min(probs[k], 1-probs[k]) >= 5) {
        se <- sqrt(probs[k]*(1-probs[k])/ess[k])
        bounds <- stats::quantile(values, probs[k]+c(-se,se), names = FALSE, type = 8)
        mcse[k] <- diff(bounds)/2
      }
    }
    width <- q[3L]-q[1L]
    ratio <- if (is.finite(width) && width > 0 && all(is.finite(mcse[c(1L,3L)])))
      max(mcse[c(1L,3L)])/width else NA_real_
    data.frame(parameter = parameters[j], ess_q025 = ess[1L], ess_median = ess[2L],
      ess_q975 = ess[3L], mcse_q025 = mcse[1L], mcse_median = mcse[2L],
      mcse_q975 = mcse[3L], max_mcse_interval_ratio = ratio, row.names = NULL)
  }))
}
