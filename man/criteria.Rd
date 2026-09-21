\name{criteria}
\alias{criteria}
\alias{WAIC}
\alias{LOOIC}
\alias{DIC}
\alias{compare_models}
\alias{print.fitdistrBayes_criteria}
\alias{as.data.frame.fitdistrBayes_criteria}
\alias{print.fitdistrBayes_comparison}
\title{Optional Bayesian Information Criteria and Model Comparison}
\description{
Computes WAIC, Pareto-smoothed importance-sampling leave-one-out cross-validation
(PSIS-LOO), and DIC from existing posterior draws. Both complete and independently
right-censored fits use the pointwise observed-data likelihood.
No sampling or automatic refitting is performed.
}
\usage{
criteria(object, methods = c("waic", "looic", "dic"), cores = 1L)
WAIC(object)
LOOIC(object, cores = 1L)
DIC(object)
compare_models(..., criterion = "looic", cores = 1L)
\method{print}{fitdistrBayes_criteria}(x, digits = 4L, ...)
\method{as.data.frame}{fitdistrBayes_criteria}(x, row.names = NULL,
  optional = FALSE, ...)
\method{print}{fitdistrBayes_comparison}(x, digits = 4L, ...)
}
\arguments{
  \item{object}{A fitted \code{fitdistrBayes} or \code{fitcensBayes} object
  retaining its log-likelihood callback (\code{store_callables = TRUE}, the default).}
  \item{methods}{Character vector selecting \code{"waic"}, \code{"looic"},
  and/or \code{"dic"}. Names ignore case; \code{"loo"} aliases \code{"looic"}.}
  \item{cores}{Number of PSIS-LOO computation cores: one (default) or two.
  This does not change the number of sampling chains.}
  \item{criterion}{One criterion used to compare all supplied models.}
  \item{...}{For \code{compare_models}, two or more preferably named fits or
  previously computed \code{fitdistrBayes_criteria} objects. In methods,
  additional arguments passed to the corresponding data-frame method.}
  \item{x}{A criteria or comparison result, as appropriate.}
  \item{digits}{Printing precision.}
  \item{row.names, optional}{Passed to \code{as.data.frame}.}
}
\details{
The fitting functions default to \code{criteria = FALSE}. They do not evaluate
the extra pointwise likelihood matrix, load \pkg{loo}, or calculate these
criteria unless requested. \code{criteria = TRUE} adds post-processing after
the chains are complete; a vector such as \code{criteria = c("waic", "dic")}
selects a subset. The chains and random-number state are unaffected by the
built-in criterion computations. Existing convergence diagnostics are separate
and remain available regardless of this option.

Post-processing uses all saved draws and temporarily requires a matrix with
one row per draw and one column per observation, approximately
\eqn{8 S n} bytes, plus working copies. A combined call evaluates this matrix
only once; the result does not retain it. Each on-demand call recomputes its
criteria, so use the combined call for efficiency. PSIS-LOO alone requires
the suggested package \pkg{loo}; WAIC and DIC use base \R.
An automatic request for LOO without \pkg{loo}, or for any criterion without
stored callbacks, fails before sampling starts.
If a numerical error occurs in optional post-processing after sampling,
the fitting functions retain the fitted object and chains, issue a warning,
and mark the requested criteria unavailable. The error is recorded in
\code{fit$criteria$diagnostics$computation_error}. Direct calls to
\code{criteria()} still report invalid likelihood matrices as errors.
Numerically unrepresentable aggregate scores are never labelled reliable.

For exact observations, \eqn{l_{si}=\log f(x_i\mid\theta_s)}. For censored
observations, \eqn{l_{si}=\log S(x_i\mid\theta_s)} instead. The censoring
mechanism is assumed independent and identical across compared models;
its parameter-free factors are omitted. Imputed lifetimes never replace
the observed likelihood in these criteria, including for augmentation fits.
For discrete data, censoring means the strict event \eqn{T>x_i}.

WAIC uses
\deqn{p_{WAIC}=\sum_i \mathrm{Var}_s(l_{si}),\qquad
 WAIC=-2\sum_i\left[\log\left(S^{-1}\sum_s e^{l_{si}}\right)
                         -\mathrm{Var}_s(l_{si})\right].}
The variance uses denominator \eqn{S-1}. Pointwise variances above 0.4
produce a warning. PSIS-LOO uses \code{loo::loo}, with chain-aware relative
effective sample sizes (one for independent exact draws). Its result includes
Pareto-k and effective-sample-size diagnostics. The diagnostic threshold is
\eqn{\min(0.7,1-1/\log_{10}(S))}. Large Pareto-k values indicate unreliable
importance sampling; exact refits and moment matching are not performed.
\eqn{LOOIC=-2\,elpd_{loo}}.

For a positive continuous distribution, a record \code{x = 0, status = 0}
is the sure event \eqn{T>0}. Its pointwise log predictive value, penalty,
and Monte Carlo error are exactly zero. Such rows are handled analytically
and restored to the full observation order after PSIS on the informative rows.
Their Pareto-k entry is \code{-Inf} as a sentinel, not an estimated tail
shape. Their indices are stored in \code{loo_zero_information_observations}
and in an attribute of the extended \code{psis_loo} result. This exception
does not apply to discrete censoring at zero.

With improper objective priors, a proper full-data posterior does not imply
proper leave-one-out training posteriors. The registered sufficient conditions
are checked for every training set before PSIS-LOO. Under censoring these are
applied to the exact-event subset; the exponential case is checked analytically.
If any training posterior is not certified, LOOIC is returned as \code{NA}
with an explanation, not an unchecked numerical score. These sufficient
checks can be conservative. LOOIC is not certified automatically for models
or priors defined by the user; WAIC and DIC rely on the user's correctly
normalized densities, proper posterior, and moment declarations.

DIC uses \eqn{D(\theta)=-2\sum_i\log p(y_i\mid\theta)},
\eqn{\bar D=E[D(\theta)\mid y]}, \eqn{p_D=\bar D-D(E[\theta\mid y])},
and \eqn{DIC=\bar D+p_D}. It is parameterization-dependent and additionally
assumes a finite expected deviance; the moment audit checks parameter means,
not this latter expectation. All parameter means must be certified finite.
Otherwise DIC is \code{NA}, including the registered Weibull scale and
unrestricted Student-t degrees of freedom. Medians are not substituted.
Even when parameter means exist, the likelihood at their joint mean can
be undefined, for example in a disconnected parameter space. In this case
DIC is unavailable, but other computable criteria and the chains are retained.
Negative \eqn{p_D} is flagged. DIC has no reported predictive standard error.

\code{compare_models} requires identical retained data, order, missing-row
indices, censoring indicators, and likelihood measures; continuous and
discrete likelihoods cannot be mixed. Keep measurement scales and censoring
definitions identical. For custom models, automatic identity checks cannot
verify the common likelihood measure; this remains the user's responsibility.
Lower criteria are better; differences do not represent posterior model
probabilities. For WAIC/LOOIC, \code{delta} is relative to the lowest score
and \code{se_delta} is \eqn{\sqrt{n\,\mathrm{Var}_i(d_i)}} using paired
pointwise IC differences. This is uncertainty across observations, not MCMC
error. With one observation it is undefined. Diagnostic failures are retained
and warned about; a ranking is not evidence that those failures are harmless.
}
\value{
\code{criteria}, \code{WAIC}, \code{LOOIC}, and \code{DIC} return a list of
class \code{"fitdistrBayes_criteria"}, with:
\describe{
  \item{\code{estimates}}{Data frame with criterion name, IC \code{estimate},
  its observation-level \code{se}, effective parameter count \code{p_eff},
  \code{elpd}, \code{available}, \code{reliable}, and explanatory \code{reason}.
  Unavailable criteria and undefined quantities are \code{NA}. The reliability
  flag only means the implemented diagnostic checks passed, not a guarantee.
  Scores from custom fits without a posterior-propriety declaration are
  not labelled reliable.}
  \item{\code{details}}{Pointwise WAIC components; the native \code{psis_loo}
  result (extended with analytic zero-information rows when necessary),
  relative efficiencies and Pareto diagnostics; and DIC components
  with the posterior-mean parameter vector, as applicable.}
  \item{\code{diagnostics}}{Sampling diagnostics, per-observation LOO
  propriety certification, and captured warning messages.}
  \item{\code{model}, \code{prior}, \code{nobs}, \code{ndraws}}{Model/prior
  identifiers and numbers of observations and saved posterior draws.}
  \item{\code{observation}}{Retained data, status, omitted-row indices,
  and likelihood measure used to prevent incompatible comparisons.}
}
\code{compare_models} returns class \code{"fitdistrBayes_comparison"}, a list
with sorted data frame \code{table}, named criterion \code{results},
\code{criterion}, and the lowest-score \code{reference} label. The table
contains model labels, estimates, differences, paired difference SEs, individual
SEs, effective parameter counts, and diagnostic reliability flags.
The print methods return their argument invisibly after printing.
The data-frame method returns the \code{estimates} data frame.
}
\references{
Gelman A, Hwang J, Vehtari A (2014). Understanding predictive information
criteria for Bayesian models. \emph{Statistics and Computing}, 24, 997--1016.
\doi{10.1007/s11222-013-9416-2}.

Vehtari A, Gelman A, Gabry J (2017). Practical Bayesian model evaluation
using leave-one-out cross-validation and WAIC.
\emph{Statistics and Computing}, 27, 1413--1432.
\doi{10.1007/s11222-016-9696-4}.

Spiegelhalter DJ, Best NG, Carlin BP, van der Linde A (2002).
Bayesian measures of model complexity and fit.
\emph{Journal of the Royal Statistical Society B}, 64, 583--639.
\doi{10.1111/1467-9868.00353}.
}
\seealso{\code{\link{fitdistrBayes}}, \code{\link{fitcensBayes}},
\code{\link{log_lik}}.}
\examples{
set.seed(62)
x <- rexp(40, rate = 0.8)
fit <- fitdistrBayes(x, "exponential", "reference", seed = 63,
  iter = 1000, warmup = 500, chains = 2, criteria = c("waic", "dic"))
fit$criteria
DIC(fit)

status <- as.integer(x <= 2)
observed <- pmin(x, 2)
cens <- fitcensBayes(observed, status, "exponential", "reference",
  seed = 64, iter = 1000, warmup = 500, chains = 2)
WAIC(cens)

if (requireNamespace("loo", quietly = TRUE)) {
  LOOIC(cens)
  fit_mdi <- fitdistrBayes(x, "exponential", "mdi", seed = 65,
    iter = 1000, warmup = 500, chains = 2)
  compare_models(Reference = fit, MDI = fit_mdi, criterion = "looic")
}
}
\keyword{models}
