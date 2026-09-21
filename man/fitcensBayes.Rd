\name{fitcensBayes}
\alias{fitcensBayes}
\title{Objective Bayesian Fitting with Independent Right Censoring}
\description{
Fits registered univariate distributions to exact and independently
right-censored observations using the objective priors of the complete-data
model. This is separate from the complete-data interface \code{\link{fitdistrBayes}}.
}
\usage{
fitcensBayes(x, status, distr, prior = NULL, start = NULL, fixed = NULL,
  iter = 4000L, warmup = floor(iter / 2), thin = 1L,
  chains = 4L, seed = NULL,
  method = c("auto", "direct", "augmentation"),
  na.action = c("fail", "omit"), control = list(), ...,
  criteria = FALSE)
}
\arguments{
  \item{x}{Nonempty numeric vector of observed values: the event time when
  \code{status = 1}, or the lower censoring limit when \code{status = 0}.
  Values must satisfy the support of the selected model.}
  \item{status}{Numeric or logical vector of the same length as \code{x}.
  One (\code{TRUE}) denotes an exact event; zero (\code{FALSE}) denotes the
  strict event \eqn{T > x}. For a discrete model, a record meaning
  \eqn{T \ge 3} must be entered as \code{x = 2, status = 0}.}
  \item{distr}{One of the 20 registered distribution names, ignoring case.
  See \code{\link{fitcensBayes_models}} and the Details section.
  A density function or model specification defined by the user is not
  supported by this censored interface.}
  \item{prior}{A registered prior name. It must be supplied explicitly.
  See \code{\link{fitcensBayes_models}} for the available names by model.}
  \item{start}{Optional named list or numeric vector of initial parameter
  values on the original scale. Missing starting values use the classical
  initialization of the complete-data model applied to the exact events.
  Starts initialize computation; they do not determine propriety.}
  \item{fixed}{Optional named list of fixed parameters. The negative binomial
  requires positive \code{size}. The usual Student-t routes require positive
  \code{df}; its \code{"independence-jeffreys"} route estimates \code{df}.}
  \item{iter}{Total iterations per chain, including warmup; at least 20.}
  \item{warmup}{Number of initial iterations discarded; nonnegative and
  smaller than \code{iter}.}
  \item{thin}{Retain every \code{thin}-th iteration after warmup.
  At least ten draws per chain must remain.}
  \item{chains}{Number of independently initialized chains; at least two.}
  \item{seed}{Optional nonnegative integer seed passed to \code{set.seed()}
  before posterior simulation. The RNG state advances normally.}
  \item{method}{Sampling method. \code{"auto"} uses the original sampler
  without censoring and the observed likelihood with censoring.
  \code{"direct"} uses the observed likelihood (not necessarily independent
  simulation); \code{"augmentation"} samples latent censored lifetimes
  and updates the parameters conditionally.}
  \item{na.action}{Whether to reject missing values or omit the entire
  \code{x, status} pair. Infinite values are always rejected.}
  \item{control}{Named list of controls, described below.}
  \item{...}{Reserved for future extensions. Unused arguments are errors.}
  \item{criteria}{\code{FALSE} by default: no model-comparison computation.
  \code{TRUE} requests WAIC, PSIS-LOO, and DIC after sampling; a character
  vector selects \code{"waic"}, \code{"looic"}, or \code{"dic"}.
  Must be named in full. PSIS-LOO requires the optional \pkg{loo} package.
  See \code{\link{criteria}}.}
}
\details{
Censoring is assumed independent and non-informative about the lifetime
parameters. The observed likelihood is proportional to
\deqn{
  L(\theta) =
    \prod_{i:\delta_i=1} f(x_i\mid\theta)
    \prod_{i:\delta_i=0} S(x_i\mid\theta),
}
where \eqn{S(t\mid\theta)=P(T>t\mid\theta)}. The censoring distribution
is not estimated. Only right censoring is implemented.

The supported distributions are beta, Cauchy, chi-squared, exponential,
gamma, geometric, lognormal, logistic, negative binomial, normal, Poisson,
Student-t, Weibull, Frechet, Gumbel, Lomax, Nakagami-m,
exponential-logarithmic, Rician, and weighted Lindley. Their parameterizations
and aliases are those of \code{\link{fitdistrBayes}}. There are 56 registered
model-prior routes. Not every prior is available for every model.

\strong{Interpretation of the priors.}
The function reuses Jeffreys, reference, and maximal data information priors
derived for the complete-data sampling model, including the registered
first-rule and parameter-ordered variants. It does not claim that these
are the corresponding objective priors derived from the information matrix
of a particular censoring design.

\strong{Posterior propriety and moments.}
Since every survival factor is at most one, a proper posterior based on the
exact-event subset implies a proper censored posterior. The complete-data
sufficient conditions are therefore applied to this subset, not to the
total or imputed sample size. A finite nonnegative posterior moment also
transfers by the same bound. Failure of this sufficient criterion is not
claimed to prove impropriety; the function stops without sampling.
Uncertified moments are labelled \code{NA} and are not estimated from a
finite chain. Known infinite moments are distinguished by \code{FALSE}.
In particular, positive moments of the Weibull scale diverge under the
registered priors, even though medians and credible intervals exist.

For the exponential model, let \eqn{m=\sum_i\delta_i} and
\eqn{A=\sum_i x_i}. The rate posterior is
\eqn{\mathrm{Gamma}(m,A)} under Jeffreys/reference and
\eqn{\mathrm{Gamma}(m+2,A)} under MDI, using the shape-rate convention.
It requires \eqn{A>0} finite and positive shape. Thus all-censored
exponential data are accepted under MDI; other all-censored routes are
not certified. For geometric Jeffreys/reference, the probability posterior
is \eqn{\mathrm{Beta}(m,\sum_i x_i+n_c+1/2)}, where \eqn{n_c} is the
number censored and counts start at zero.

\strong{Computation.}
Observed-likelihood sampling uses exact simulation for the analytic cases,
marginal slice sampling with conditional scale simulation for Weibull,
and preconditioned slice sampling otherwise. Augmentation generates
\eqn{T_i\mid T_i>x_i,\theta} and uses conditional parameter updates;
normal and lognormal models have exact blocked updates. Initial values
use classical estimates from the exact events, optional optimization of
the observed log posterior, and dispersed chain starts. The geometry is
fixed during retained simulation. Imputed data are never used to certify
propriety. Difficult tails and heavy censoring can require longer runs;
inspect diagnostics rather than relying on default iteration counts.

\strong{Controls.}
\describe{
  \item{\code{rhat_threshold}, \code{ess_threshold}}{Diagnostic targets,
  default 1.01 and 400 for rank-normalized split/folded R-hat and bulk/tail ESS.}
  \item{\code{slice_width}, \code{slice_steps}, \code{max_shrink}}{Slice
  width (1), maximum stepping-out budget (40), and shrinkage budget (1000).}
  \item{\code{init_jitter}, \code{max_init_tries}}{Initial dispersion (2)
  and maximum attempts to find a finite initial target (100).}
  \item{\code{optimize_start}}{Use observed-posterior optimization for
  automatic starts (default \code{TRUE}). A supplied start keeps its center.}
  \item{\code{warn_convergence}}{Warn if MCMC diagnostic targets fail
  (default \code{TRUE}). Turning off the warning does not alter diagnostics.}
  \item{\code{store_callables}}{Store callbacks for prediction and
  log-likelihood extraction (default \code{TRUE}).}
  \item{\code{entropy_tol}, \code{entropy_exact_limit}}{Tolerance
  (\code{1e-10}) and exact-sum limit (2000) for applicable MDI calculations.}
}
Unknown controls are errors. No new mandatory packages are required.
}
\value{
A list of class \code{"fitcensBayes"}, separate from \code{"fitdistrBayes"}:
\describe{
  \item{\code{call}, \code{model}, \code{prior}}{The matched call;
  model name, parameter names, fixed values and sample size; and prior
  label, kernel, complete-data origin and posterior-propriety justification.}
  \item{\code{censoring}}{Censoring type, exact/censored counts, censored
  fraction and the non-informative censoring assumption.}
  \item{\code{initialization}, \code{engine}}{Starting procedure and
  values when applicable; selected algorithm, iteration settings and seed.}
  \item{\code{estimates}}{Named numeric vector of posterior medians.}
  \item{\code{summary}}{Data frame with one row per parameter: justified
  posterior mean and standard deviation, median, equal-tail 95 percent
  limits, R-hat, ESS, mean Monte Carlo standard error, and moment status.}
  \item{\code{moment_status}}{Data frame recording whether each mean
  and variance is finite (\code{TRUE}), infinite (\code{FALSE}), or not
  certified (\code{NA}), with explanations. Uncertified or infinite means,
  standard deviations and associated MCSEs are suppressed in the summary.}
  \item{\code{diagnostics}}{Diagnostic targets, extrema, explanatory
  messages and a \code{converged} flag. This flag denotes meeting the
  diagnostic targets, not a mathematical proof of convergence.
  \code{exact_or_independent} identifies independent posterior simulation,
  for which MCMC convergence checks do not apply.}
  \item{\code{chains}}{List of numeric matrices, one per chain, with
  retained iterations in rows and model parameters in columns.}
  \item{\code{draws}}{Combined data frame with columns \code{.chain},
  \code{.iteration}, \code{.draw}, and the parameter values.}
  \item{\code{data}, \code{status}, \code{omitted}}{Observed values and
  indicators after missing-pair removal, and omitted original indices.}
  \item{\code{control}}{Validated computational settings.}
  \item{\code{criteria}}{Present only when requested: a
  \code{"fitdistrBayes_criteria"} object based on the observed likelihood,
  not the density of imputed lifetimes. Does not change the chains.}
}
When \code{store_callables = TRUE}, internal callbacks are stored for the
observed log likelihood, predictive sampling, survival and imputation.
See \code{\link{fitcensBayes-methods}} for extraction and prediction methods.
}
\seealso{
\code{\link{fitdistrBayes}}, \code{\link{fitcensBayes_models}},
\code{\link{fitcensBayes-methods}}.
}
\examples{
set.seed(42)
lifetime <- rexp(80, rate = 0.7)
censoring <- rexp(80, rate = 0.3)
x <- pmin(lifetime, censoring)
status <- as.integer(lifetime <= censoring)
fit <- fitcensBayes(x, status, "exponential", "reference", seed = 43)
fit
confint(fit)

# Discrete censoring means strictly greater than the reported count.
fit_geom <- fitcensBayes(c(0, 1, 2, 3, 4, 5), c(1, 1, 0, 1, 0, 0),
                        "geometric", "jeffreys", seed = 44)
coef(fit_geom)

# Inspect the installed tutorial without executing long MCMC runs.
system.file("examples", "tutorial_fitcensBayes.R", package = "fitdistrBayes")
}
\keyword{models}
\keyword{survival}
