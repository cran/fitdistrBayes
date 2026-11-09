# Keep callback environments independent of the complete fitting frame.
.aft_loglik_factory <- function(X, y, status) {
  force(X); force(y); force(status)
  function(theta) .aft_loglik(theta, y, status, X)
}
