\name{fitcensBayes-methods}
\alias{fitcensBayes-methods}
\alias{print.fitcensBayes}
\alias{summary.fitcensBayes}
\alias{print.summary.fitcensBayes}
\alias{coef.fitcensBayes}
\alias{confint.fitcensBayes}
\alias{as.data.frame.fitcensBayes}
\alias{plot.fitcensBayes}
\alias{predict.fitcensBayes}
\alias{log_lik.fitcensBayes}
\alias{log_lik_cens}
\title{Methods for Objective Bayesian Fits with Right Censoring}
\description{
Summarizes posterior draws, assesses simulation diagnostics, generates
uncensored future observations, computes survival probabilities, imputes
censored lifetimes, and extracts the observed-data log likelihood.
}
\usage{
\method{print}{fitcensBayes}(x, digits = 4L, ...)
\method{summary}{fitcensBayes}(object, ...)
\method{print}{summary.fitcensBayes}(x, digits = 4L, ...)
\method{coef}{fitcensBayes}(object, ...)
\method{confint}{fitcensBayes}(object,
  parm = object$model$parameters, level = 0.95, ...)
\method{as.data.frame}{fitcensBayes}(x,
  row.names = NULL, optional = FALSE, ...)
\method{plot}{fitcensBayes}(x,
  type = c("trace", "density", "acf", "pairs"),
  pars = x$model$parameters, ...)
\method{predict}{fitcensBayes}(object,
  type = c("response", "survival", "impute"),
  times = NULL, draws = 1000L, size = 1L, seed = NULL, ...)
\method{log_lik}{fitcensBayes}(object, draws = NULL, seed = NULL, ...)
log_lik_cens(object, draws = NULL, seed = NULL)
}
\arguments{
  \item{x, object}{A \code{"fitcensBayes"} object, or for its summary
  print method an object of class \code{"summary.fitcensBayes"}.}
  \item{digits}{Number of significant digits printed.}
  \item{parm}{Parameter names for credible intervals.}
  \item{level}{Credible probability, strictly between zero and one.}
  \item{row.names, optional}{Compatibility arguments, currently ignored.
  Stored draw identifiers are retained.}
  \item{type}{Plot or prediction type. See Details for prediction.}
  \item{pars}{Parameter names to plot. Pairs plots need at least two.}
  \item{times}{Finite numeric evaluation times for survival prediction.
  Must be \code{NULL} for other prediction types.}
  \item{draws}{Number of posterior samples to select. For
  \code{log_lik}, \code{NULL} uses all stored draws. Prediction samples
  posterior rows without replacement unless more rows are requested
  than are stored.}
  \item{size}{Number of future observations per draw for response
  prediction. A positive integer; irrelevant to the other prediction types.}
  \item{seed}{Optional seed passed to \code{set.seed()} before simulation
  or posterior row selection. RNG state advances normally.}
  \item{...}{Additional graphical arguments for \code{plot}. Unused
  arguments to prediction or log-likelihood extraction are errors;
  other methods accept them for generic compatibility.}
}
\details{
The methods dispatch separately from those for uncensored fits.
\code{predict(type = "response")} generates future uncensored values
from the lifetime distribution, not an observed value/indicator pair.
\code{predict(type = "survival")} returns
\eqn{P(T>t\mid\theta)} at each requested time for sampled parameter values.
\code{predict(type = "impute")} samples latent lifetimes conditionally on
each censored record. Every imputed value is strictly above its observed
limit, including in discrete models. Imputation does not change the fit.
No imputation is possible without censored records.

\code{\link{log_lik}} dispatches to the censored method.
\code{log_lik_cens} is an explicit convenience wrapper with the same
result. Each observation contributes either its log density or its
log survival, according to the observed indicator, not a completed-data
log density. The censoring mechanism's own likelihood is excluded under
the non-informative censoring assumption.

Prediction and likelihood extraction require callbacks stored by
\code{control = list(store_callables = TRUE)}, the fitting default.
}
\value{
\describe{
  \item{\code{print}, \code{print.summary}}{The supplied object,
  invisibly, after printing model, censoring, posterior and diagnostic
  information. The summary print method also displays the moment audit.}
  \item{\code{summary}}{A list of class \code{"summary.fitcensBayes"}
  containing \code{model}, \code{prior}, \code{summary},
  \code{moment_status}, \code{diagnostics}, \code{initialization},
  \code{engine}, and \code{censoring} from the fit.}
  \item{\code{coef}}{A named numeric vector of posterior medians.}
  \item{\code{confint}}{A numeric matrix with one row per requested
  parameter and two columns of equal-tail posterior credible limits.
  These are credible intervals, not frequentist confidence intervals.}
  \item{\code{as.data.frame}}{A data frame of retained posterior draws
  with chain, iteration and draw identifiers and parameter columns.}
  \item{\code{plot}}{The original fitted object, invisibly. Produces
  trace, density, autocorrelation or pairs plots as a side effect.}
  \item{\code{predict}}{Always a numeric matrix with \code{draws} rows.
  For response prediction it has \code{size} columns of future values;
  for survival it has \code{length(times)} columns of conditional
  survival probabilities; for imputation it has one column per censored
  individual, in their order in the retained data. Rows incorporate
  posterior parameter uncertainty.}
  \item{\code{log_lik}, \code{log_lik_cens}}{A numeric matrix with
  selected posterior draws in rows and retained observations in columns.
  Entries are observed-data log-likelihood contributions.}
}
}
\seealso{\code{\link{fitcensBayes}}, \code{\link{fitcensBayes_models}}}
\examples{
fit <- fitcensBayes(c(0.3, 0.7, 1.2, 1.5, 2),
                   c(1, 1, 0, 1, 0), "exponential", "jeffreys",
                   seed = 15)
summary(fit)
coef(fit)
confint(fit, level = 0.9)
head(as.data.frame(fit))
predict(fit, type = "survival", times = c(1, 2), draws = 3, seed = 16)
predict(fit, type = "impute", draws = 3, seed = 17)
predict(fit, type = "response", draws = 3, size = 2, seed = 18)
log_lik(fit, draws = 3, seed = 19)
plot(fit, type = "density")
}
\keyword{models}
