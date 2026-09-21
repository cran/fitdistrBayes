\name{fitdistrBayes-methods}
\alias{print.fitdistrBayes}
\alias{summary.fitdistrBayes}
\alias{print.summary.fitdistrBayes}
\alias{coef.fitdistrBayes}
\alias{confint.fitdistrBayes}
\alias{as.data.frame.fitdistrBayes}
\alias{plot.fitdistrBayes}
\alias{predict.fitdistrBayes}
\alias{log_lik}
\alias{log_lik.fitdistrBayes}
\title{Methods for Objective Bayesian Distribution Fits}
\usage{
\method{print}{fitdistrBayes}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{summary}{fitdistrBayes}(object, ...)
\method{print}{summary.fitdistrBayes}(x,
  digits = max(3L, getOption("digits") - 3L), ...)
\method{coef}{fitdistrBayes}(object, ...)
\method{confint}{fitdistrBayes}(object, parm = object$model$parameters,
  level = 0.95, ...)
\method{as.data.frame}{fitdistrBayes}(x, row.names = NULL,
  optional = FALSE, ...)
\method{plot}{fitdistrBayes}(x,
  type = c("trace", "density", "acf", "pairs"),
  pars = x$model$parameters, ...)
\method{predict}{fitdistrBayes}(object, draws = 1000L, size = 1L,
  seed = NULL, ...)
log_lik(object, ...)
\method{log_lik}{fitdistrBayes}(object, draws = NULL, seed = NULL, ...)
}
\arguments{
  \item{x, object}{A fitted \code{"fitdistrBayes"} object. The
  \code{log_lik} generic also accepts \code{"fitcensBayes"} objects;
  see \code{\link{fitcensBayes-methods}} for their separate methods.}
  \item{digits}{Number of printed significant digits.}
  \item{parm}{Parameter names for interval extraction.}
  \item{level}{Credible level.}
  \item{row.names, optional}{Arguments for data-frame conversion.}
  \item{type}{Diagnostic plot type.}
  \item{pars}{Optional subset of parameters.}
  \item{draws}{Number of posterior or predictive draws.}
  \item{size}{Number of observations in each predictive data set.}
  \item{seed}{Optional reproducibility seed passed to \code{set.seed()} before
  sampling posterior draws or posterior predictive observations.}
  \item{...}{Additional arguments.}
}
\value{
The returned value depends on the method:
\itemize{
  \item \code{print.fitdistrBayes()} returns \code{x}, invisibly, and prints
  the model, prior, computational engine, posterior summaries, and diagnostic
  status.
  \item \code{summary.fitdistrBayes()} returns an object of class
  \code{"summary.fitdistrBayes"}. It is a list containing the original call;
  model, prior, initialization, and computational-engine records; the
  posterior summary table; the posterior-moment audit; convergence
  diagnostics; and the available prediction and log-likelihood capabilities.
  \item \code{print.summary.fitdistrBayes()} returns \code{x}, invisibly, and
  prints the contents of the summary object.
  \item \code{coef.fitdistrBayes()} returns a named numeric vector containing
  the posterior median of each estimated parameter. These medians are the
  default point estimates used by the package.
  \item \code{confint.fitdistrBayes()} returns a numeric matrix with one row
  per requested parameter and two columns containing the lower and upper
  equal-tail posterior credible limits.
  \item \code{as.data.frame.fitdistrBayes()} returns a data frame of the
  combined posterior draws. The columns \code{.chain}, \code{.iteration}, and
  \code{.draw} identify the simulation draw; the remaining columns contain
  parameter values.
  \item \code{plot.fitdistrBayes()} returns \code{x}, invisibly, and produces
  the requested diagnostic plot as a side effect.
  \item \code{predict.fitdistrBayes()} returns posterior predictive
  observations. It is a numeric vector of length \code{draws} when
  \code{size = 1}, and otherwise a numeric matrix with \code{draws} rows and
  \code{size} columns. Each row is one replicated data set generated after
  selecting a posterior draw.
  \item \code{log_lik()} and \code{log_lik.fitdistrBayes()} return a numeric
  matrix whose rows correspond to posterior draws and whose columns
  correspond to observations. Each entry is that observation's
  log-likelihood contribution at the selected posterior draw.
}
}
\description{
Methods for inspecting posterior summaries and draws, extracting posterior
medians and credible intervals, producing standard diagnostic plots,
simulating posterior predictive observations, and computing pointwise
log-likelihood matrices.
Prediction and pointwise log-likelihood are computed on demand. They are
unavailable when the fitted object was created with
\code{control = list(store_callables = FALSE)}; inspect
\code{object$capabilities} before calling them in reusable workflows.
}
