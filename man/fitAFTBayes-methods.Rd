\name{fitAFTBayes-methods}
\alias{fitAFTBayes-methods}
\alias{print.fitAFTBayes}
\alias{summary.fitAFTBayes}
\alias{print.summary.fitAFTBayes}
\alias{coef.fitAFTBayes}
\alias{confint.fitAFTBayes}
\alias{as.data.frame.fitAFTBayes}
\alias{log_lik.fitAFTBayes}
\alias{predict.fitAFTBayes}
\alias{plot.fitAFTBayes}
\alias{nobs.fitAFTBayes}
\alias{formula.fitAFTBayes}
\alias{model.matrix.fitAFTBayes}
\title{Methods for Objective Bayesian AFT Fits}
\description{Summaries, time ratios, prediction, diagnostics and observed-data
log likelihoods for \code{\link{fitAFTBayes}} objects.}
\usage{
\method{print}{fitAFTBayes}(x, digits = 4L, ...)
\method{summary}{fitAFTBayes}(object, ...)
\method{print}{summary.fitAFTBayes}(x, digits = 4L, ...)
\method{coef}{fitAFTBayes}(object,
  type = c("regression", "all", "time_ratio"), ...)
\method{confint}{fitAFTBayes}(object, parm = NULL, level = 0.95,
  type = c("regression", "all", "time_ratio"), ...)
\method{as.data.frame}{fitAFTBayes}(x, row.names = NULL, optional = FALSE, ...)
\method{log_lik}{fitAFTBayes}(object, draws = NULL, seed = NULL, ...)
\method{predict}{fitAFTBayes}(object, newdata = NULL,
  type = c("survival", "linear"), times = NULL, summary = TRUE,
  level = 0.95, draws = NULL, seed = NULL, ...)
\method{plot}{fitAFTBayes}(x, type = c("trace", "density", "acf", "pairs"),
  pars = x$model$parameters, ...)
\method{nobs}{fitAFTBayes}(object, ...)
\method{formula}{fitAFTBayes}(x, ...)
\method{model.matrix}{fitAFTBayes}(object, ...)
}
\arguments{
\item{x, object}{An AFT fit, or its summary for the summary print method.}
\item{digits}{Printing precision.}
\item{type}{For coefficients/intervals: regression beta only, all parameters,
or exp(beta) for non-intercept design columns. For prediction: survival
probabilities or the linear log-time location \eqn{x^T\beta}. For plots:
trace, density, autocorrelation or pairs.}
\item{parm}{Parameter names or indices in the selected coefficient type;
\code{NULL} selects all. Regression names match design columns; \code{"all"}
uses the internal \code{beta:} prefix.}
\item{level}{Equal-tailed credible interval probability, strictly between zero and one.}
\item{row.names, optional}{Accepted for data-frame method compatibility.}
\item{draws}{\code{NULL} uses all saved draws. A positive integer randomly
selects that many draws, with replacement only if greater than available.}
\item{seed}{Optional seed used when selecting a subset of draws.}
\item{newdata}{Prediction data frame. \code{NULL} uses retained training rows.
Factor levels and transformed formula bases use the fitted metadata.
Missing values and new factor levels are errors.}
\item{times}{Nonnegative survival prediction times; zero and infinity are
allowed. Required for survival, unused for linear prediction. Units must
match the fitted response.}
\item{summary}{Return summaries rather than individual posterior draws.}
\item{pars}{Internal parameter names to plot.}
\item{...}{Printing/plotting arguments where appropriate. Extra prediction
and log-likelihood arguments are rejected.}
}
\details{
The term \emph{credible interval} is intended even though the standard generic
is named \code{confint}. Coefficients are posterior medians, not posterior
means. Time ratios describe a unit change in a design column with other
columns fixed. With interactions or nonlinear transformations, such a change
may not correspond to a realizable covariate contrast; construct the desired
linear contrast explicitly from the posterior regression draws.

Survival summaries contain \code{row}, \code{time}, \code{estimate},
\code{median}, \code{lower}, and \code{upper}. \code{estimate} is the posterior
mean survival probability, equivalently the predictive survival probability.
The median and interval summarize the posterior distribution of
\eqn{S(t\mid x,\theta)} and are not prediction intervals for a lifetime.
Linear summaries use the posterior median as \code{estimate}, avoiding
uncertified parameter means. Predictions do not repair poor sampling.

Observed-data log likelihoods are log density for events and log survival for
censoring, on the original response-time scale. No imputed lifetime density
is substituted. \code{WAIC}, \code{LOOIC}, \code{DIC} and
\code{compare_models} support this class. Each leave-one-out training design
must satisfy the sufficient propriety conditions. DIC is unavailable unless
every necessary parameter mean is certified finite, and is therefore
unavailable under GG reference/simple. Other families use their own
moment and deletion-propriety checks. No automatic leave-one-out refitting is
performed when Pareto diagnostics fail.
}
\value{Print and plot methods return their input invisibly. Summary returns
\code{summary.fitAFTBayes}; coefficient and interval methods return a vector
and matrix. The data-frame method returns all draws with chain/iteration
identifiers. Log likelihood is a draws-by-observations matrix. Unsummarized
linear predictions are draws-by-rows matrices; unsummarized survival
predictions are arrays indexed by draw, newdata row and time. Other prediction
results are data frames. Metadata methods return the retained sample size,
formula or design matrix.}
\seealso{\code{\link{fitAFTBayes}}, \code{\link{criteria}}}
