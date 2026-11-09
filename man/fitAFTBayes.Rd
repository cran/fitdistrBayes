\name{fitAFTBayes}
\alias{fitAFTBayes}
\title{Objective Bayesian Accelerated Failure Time Regression}
\description{Fits nine parametric AFT lifetime families with complete or
independently right-censored times, objective priors and model-specific
propriety and moment checks. See \code{\link{fitAFTBayes_models}} for the catalogue.}
\usage{
fitAFTBayes(formula, data, distr = "gg", prior = "reference",
  iter = 4000L, warmup = floor(iter/2),
  thin = 1L, chains = 4L, seed = NULL,
  na.action = c("fail", "omit"), control = list(), ..., criteria = FALSE,
  df = NULL, shape = NULL)
}
\arguments{
\item{formula}{A two-sided formula. The response is either positive complete
times or \code{survival::Surv(time, status)} with ordinary right censoring.
An intercept is required. Factors, interactions and standard formula
transformations are supported. Offsets are not supported.}
\item{data}{Data frame containing the formula variables.}
\item{distr}{\code{"gg"} (default, Prentice generalized gamma),
\code{"weibull"}, \code{"lognormal"}, \code{"loglogistic"},
\code{"loglaplace"}, \code{"logt"}, \code{"exponential"},
\code{"rayleigh"}, or \code{"gamma"}. Names ignore case, spaces and hyphens.
Aliases include \code{"gengamma"}, \code{"generalized gamma"},
\code{"loggaussian"}, \code{"log-student-t"} and \code{"t"}.
Here \code{"t"} means Student errors in log time, not Student lifetimes.}
\item{prior}{\code{"reference"} (default). GG additionally accepts
\code{"simple"}; all other families accept \code{"jeffreys"} for the joint
Jeffreys measure. For these other families, \code{"independence-jeffreys"}
and \code{"joint-jeffreys"} are aliases for reference and Jeffreys.
Jeffreys is rejected for GG because that joint posterior is improper.
There is no MDI route for this regression interface.}
\item{iter}{Total sweeps per chain, including warmup; at least 20. Each
sweep updates every sampling coordinate once.}
\item{warmup}{Iterations discarded per chain.}
\item{thin}{Keep every \code{thin}-th post-warmup iteration; default one.}
\item{chains}{Number of chains, at least two.}
\item{seed}{Optional integer random-number seed. As in the other fitting
interfaces, setting this resets R's random-number generator.}
\item{na.action}{Fail on missing response/covariates, or omit the corresponding
whole model-frame rows. Retained observations and omitted indices are stored.}
\item{control}{Named list. Defaults: \code{rhat_threshold=1.01},
\code{ess_threshold=400}, \code{warn_convergence=TRUE},
\code{store_callables=TRUE}, \code{slice_width=1},
\code{slice_steps=40L}, \code{max_shrink=1000L},
\code{optimize_maxit=1000L}. Unknown options are errors.}
\item{...}{Reserved; unused arguments are errors.}
\item{criteria}{\code{FALSE} by default. Optionally \code{TRUE} or a subset of
\code{c("waic","looic","dic")}; see \code{\link{criteria}}.}
\item{df}{Positive fixed degrees of freedom, required for \code{"logt"}.
Degrees of freedom are not estimated. Unused for other models and rejected
if supplied to them. This is a model parameter, not a prior hyperparameter.}
\item{shape}{Optional positive fixed shape for \code{"gamma"}. The default
\code{NULL} estimates the shape over its unrestricted positive support.
Rejected for other models.}
}
\details{
\strong{Fixed-error and gamma models.}
Weibull, lognormal, log-logistic, log-Laplace and log-Student-t use
\eqn{\log T_i=x_i^T\beta+\sigma Z_i}. With p coefficients including the
intercept, reference is \eqn{1/\sigma} and joint Jeffreys is
\eqn{\sigma^{-(p+1)}} with respect to \eqn{d\beta\,d\sigma}.
Exponential and Rayleigh fix sigma at 1 and 1/2, respectively; both prior
names then give the same flat coefficient measure. Weibull lifetime shape
is \eqn{1/\sigma}, and its lifetime scale is \eqn{\exp(x_i^T\beta)}.

Gamma uses shape k and rate \eqn{k\exp(-x_i^T\beta)}, so the conditional
mean is \eqn{\exp(x_i^T\beta)}. There is no additional free sigma.
The unrestricted reference and joint Jeffreys priors are respectively
\eqn{\sqrt{\psi_1(k)-1/k}} and
\eqn{k^{p/2}\sqrt{\psi_1(k)-1/k}}, flat in beta.
The gamma model estimates k unless \code{shape} is supplied; a fixed k
makes both priors flat in beta. No truncation or proper shape hyperprior
is silently introduced.

Unknown-scale exponential-tail families and unknown-shape gamma require
more failures than coefficients, full failure-design rank and no exact
interpolation of log failure times. Fixed-scale families need only full
failure-design rank. For Student-t, an additional conservative certificate
partitions failures into an independent p-row design and disjoint
noninterpolable blocks of p+1 rows. With B blocks and prior power a, it
requires \eqn{B(\mathrm{df}+1)>d+a-p-1}, where d is the number of failures.
This protects against tied and partially interpolated samples. Failure to
pass a sufficient check is not a declaration that every such posterior is
improper. Details and boundary proofs are in \file{doc/aft-models-theory.txt}.

All families use the complete-data prior with the observed censored likelihood.
For the four exponentially bounded unknown-scale errors, absolute parameter
moments of positive order r exist exactly when \eqn{r<d+a-p-1} under the
certified design. Student coefficient moments additionally require
\eqn{r<\mathrm{df}} for this sufficient certificate; others may be unknown.
Gamma coefficient moments have thresholds d-p (reference) and d-p/2
(Jeffreys), while all positive shape moments are finite. Fixed-scale
families have all polynomial coefficient moments. Reported means/SDs and
DIC respect these flags; \code{coef()} continues to return medians.
No posterior mean of a time ratio or predictive lifetime is inferred from
the existence of a conditional lifetime mean.

Complete lognormal regression uses exact inverse-gamma/normal simulation.
Intercept-only exponential/Rayleigh with right censoring and complete
intercept-only known-shape gamma use exact transformed-gamma simulation.
For these routes, iter/warmup/thin determine the number of independent
saved draws, without an MCMC burn-in requirement. Other new routes use
quantile-centered, mixed-fixed-basis slice sweeps, initialized from failure
least squares and numerical optimization. A fixed positive design-based
geometry replaces an unusable Hessian, especially for nonsmooth Laplace
targets; it does not change the posterior. Covariates are internally scaled
and all saved coefficients are returned on the original scale.

\strong{Generalized gamma.}
For GG the model is \eqn{\log T_i=x_i^T\beta+\sigma W_i}, with \eqn{\sigma>0}.
For \eqn{Q\ne0}, \eqn{G\sim\mathrm{Gamma}(Q^{-2},1)} and
\eqn{W=\log(Q^2G)/Q}. At \eqn{Q=0}, \eqn{W} is standard normal.
Thus \eqn{Q=0} gives lognormal AFT, and \eqn{Q=1} gives Weibull AFT.
This route estimates Q. The lognormal and Weibull routes are separately
fitted submodels with their own priors, not restrictions of the GG prior.

The reference prior is \eqn{\omega(Q)/\sigma}, flat in beta, from the
complete-data reference construction. It is applied to the observed
right-censored likelihood; it is not asserted to be the reference prior
of every censoring experiment. With \eqn{s=Q^2} and
\eqn{u=\psi_1(1+1/s)},
\deqn{\omega(Q)=|Q|^{-3}\{[s^2(u-1)+u^2]/(s+u)\}^{1/2},}
and \eqn{\omega(0)=1/\sqrt{24}}. Simple uses
\eqn{1/[\sigma(1+Q^2)]}, also flat in beta. Both joint priors are improper
and contain no user-specified hyperparameters. The slice controls govern
numerical computation, not the prior.

Before sampling, both GG routes require more failures than regression
coefficients (including the intercept), a full-rank failure design, and log
failure times that are not exactly interpolated by that design. These are
sufficient conditions; their failure is reported as lack of certification,
not a proof of impropriety for every censored sample.
The numerical rank/interpolation tolerance is relative, \eqn{10^{-10}}.
Failure times must be positive; censoring at zero is allowed and contributes
exactly one to the likelihood. Left truncation, interval/left censoring,
multi-state responses and observation weights are not implemented.
Censoring must be conditionally independent and non-informative given the
covariates; this assumption is not established by the numerical checks.

The sampler uses stepping-out and shrinkage slice updates in coordinates
\eqn{m=\beta_0+\sigma w_c(Q)}, the regression slopes,
\eqn{\log\tau=\log\sigma+\log\cosh\zeta}, and \eqn{\zeta=\mathrm{asinh}(Q)}.
Here \eqn{w_c(Q)} is a numerically stable approximation to the standardized
log-time quantile at \eqn{c=\max(0.05,\min(0.5,r/(2n)))}, with r failures.
The intercept translation has unit Jacobian; both target evaluation and
saved draws use its same inverse. Thus centering, including its numerical
approximation, does not change the prior or likelihood.

Each sweep randomly chooses either the full posterior-mode Cholesky basis
or its marginal-scale diagonal basis, with equal probabilities, then updates
all coordinates in random order. Both bases are fixed before sampling.
The Hessian supplies computational geometry only, never a Wald covariance.
The width of each slice bracket adapts locally through valid stepping-out
and shrinkage updates, without changing the target. \code{slice_steps}
limits bracket expansion and \code{max_shrink} limits rejected proposals;
exhaustion of numerical resolution or shrinkage attempts is an error.

Warmup sweeps are discarded. There is no automatic refitting if diagnostics
fail. Increase iterations and inspect trace plots, especially for weakly
identified shapes and heavy censoring. An iteration now contains a full
slice sweep and has a different computational cost from the Metropolis
iteration in the local 0.7.0 development release.

For GG reference and simple, \eqn{E[|Q|]} is infinite: posterior medians and
quantiles are reported and the mean of Q is suppressed. Other uncertified
moments are also suppressed. DIC is unavailable under both priors because
the required parameter means are not certified finite. No finite predictive
lifetime mean is asserted. Survival probabilities are bounded quantities.
The mathematical specification is provided in
the installed \file{doc/gg-aft-theory.txt}.

The \code{quantile_precision} table reports approximate Monte Carlo standard
errors for the 2.5\%, 50\% and 97.5\% quantiles and their empirical-CDF
indicator ESS. CDF uncertainty \eqn{\{a(1-a)/\mathrm{ESS}_a\}^{1/2}} is
mapped back through empirical quantiles; no parameter mean or variance is
needed. Estimates are withheld for stuck chains or fewer than five effective
observations in the relevant tail. The endpoint MCSE divided by the 95\%
interval width helps assess accuracy, especially for Q. These approximations
are not convergence proofs and may be unreliable before mixing is adequate.
The \code{converged} flag refers to the stated R-hat/bulk/tail ESS thresholds;
it does not certify a user-specific number of accurate interval digits.
}
\value{An object of class \code{fitAFTBayes}, containing the formula and
prediction terms, design and factor metadata, retained times/status, posterior
chains and a draws data frame, medians and quantile summaries, moment status,
propriety assessment, initialization, sampler settings and convergence
diagnostics. The engine records the centering probability and number of
coordinate updates per iteration. Diagnostics record target evaluations
and movement fractions per chain; the latter are not Metropolis acceptance
rates or convergence certificates. Parameter columns are named \code{beta:}
followed by the design column name, then the free parameters: \code{sigma}
and \code{Q} for GG, \code{sigma} for other unknown-scale errors, or
\code{shape} for gamma with estimated shape. Fixed-scale/shape routes
contain only beta columns. Model metadata retain fixed arguments.
Criteria are added only if requested.}
\references{
Neal, R. M. (2003). Slice sampling. \emph{The Annals of Statistics},
31, 705--767. \doi{10.1214/aos/1056562461}.

Prentice, R. L. (1974). A log gamma model and its maximum likelihood estimation.
\emph{Biometrika}, 61, 539--544. \doi{10.1093/biomet/61.3.539}.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., and Buerkner, P.-C. (2021).
Rank-Normalization, Folding, and Localization: An Improved R-hat for Assessing
Convergence of MCMC. \emph{Bayesian Analysis}, 16, 667--718.
\doi{10.1214/20-BA1221}.

Ramos, P. L., Bravo Palma, G., Moala, F. A., and Dey, D. K. (2026).
Generalized Gamma Accelerated Failure Time Regression for Right-Censored
Survival Data: An Objective Bayesian Analysis. Unpublished manuscript.
}
\seealso{\code{\link{fitAFTBayes-methods}}, \code{\link{fitcensBayes}},
\code{\link{fitAFTBayes_models}}, \code{\link{criteria}}, \code{\link[survival]{Surv}}}
\examples{
library(survival)
set.seed(7)
d <- data.frame(age10 = rnorm(60), treatment = rep(0:1, 30))
t <- exp(0.5 - 0.3*d$age10 + 0.4*d$treatment + rnorm(60, sd = 0.7))
d$time <- pmin(t, 3)
d$status <- as.integer(t <= 3)
# Short illustration; use longer chains and inspect diagnostics for analysis.
fit <- fitAFTBayes(Surv(time, status) ~ age10 + treatment, d,
  iter = 100, warmup = 50, chains = 2, seed = 8,
  control = list(warn_convergence = FALSE))
coef(fit, type = "time_ratio")
confint(fit, type = "time_ratio")
predict(fit, data.frame(age10 = 0, treatment = 0:1), times = c(1, 3), draws = 50)
# Fixed objective alternative, without prior hyperparameters.
simple <- fitAFTBayes(Surv(time, status) ~ age10 + treatment, d,
  prior = "simple",
  iter = 100, warmup = 50, chains = 2, seed = 8,
  control = list(warn_convergence = FALSE))
simple$diagnostics

# Other families share the formula/method interface.
fitAFTBayes_models()
d$complete_time <- t
normal <- fitAFTBayes(complete_time ~ age10 + treatment, d, distr = "lognormal",
  iter = 200, warmup = 0, chains = 2, seed = 9,
  control = list(warn_convergence = FALSE))
coef(normal)
gamma <- fitAFTBayes(Surv(time, status) ~ age10 + treatment, d,
  distr = "gamma", prior = "jeffreys", iter = 100, warmup = 50,
  chains = 2, seed = 10, control = list(warn_convergence = FALSE))
coef(gamma, type = "all")
}
