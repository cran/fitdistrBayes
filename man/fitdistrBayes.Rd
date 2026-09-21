\name{fitdistrBayes}
\alias{fitdistrBayes}
\title{Objective Bayesian Distribution Fitting}
\usage{
fitdistrBayes(x, distr, prior = NULL, start = NULL, fixed = NULL,
  iter = 4000L, warmup = floor(iter / 2), thin = 1L,
  chains = 4L, seed = NULL, na.action = c("fail", "omit"),
  control = list(), ..., criteria = FALSE)
}
\arguments{
  \item{x}{Numeric vector of observations.}
  \item{distr}{A recognized distribution name, matched without regard to
  case, a density function evaluated at its first argument, or an object made
  by \code{fitdistrBayes_model()}.}
  \item{prior}{A registered objective-prior name or, for a custom density
  function, a prior function on the documented parameterization. Omit this
  argument when \code{distr} is a \code{fitdistrBayes_model} object.}
  \item{start}{Optional named starting values. If omitted for a built-in
  non-exact route, closed-form classical estimates are used: ordinary
  moments where available, quantile matching for models without the required
  moments, and L-moments for Weibull, Frechet, and Lomax. The
  Exponential-Logarithmic start solves one scalar moment equation, and the
  weighted Lindley start uses its closed-form likelihood estimator with a
  numerical maximum-likelihood fallback. Required for custom densities and
  ignored with a warning for exact independent posterior simulation.}
  \item{fixed}{Optional named list of fixed parameters.}
  \item{iter}{Total iterations per chain, including warmup.}
  \item{warmup}{Warmup iterations per chain.}
  \item{thin}{Positive thinning interval.}
  \item{chains}{Number of chains.}
  \item{seed}{Optional reproducibility seed passed to \code{set.seed()} before
  posterior simulation.}
  \item{na.action}{Either \code{"fail"} or \code{"omit"}.}
  \item{control}{Named list of computational, diagnostic, storage, and custom
  model controls described in Details.}
  \item{...}{Additional fixed arguments forwarded to a custom density and
  predictive generator. Unused arguments are errors for built-in distributions.}
  \item{criteria}{Optional post-processing: \code{FALSE} (default) does no
  information-criterion computation; \code{TRUE} requests WAIC, PSIS-LOO,
  and DIC. A character vector selects \code{"waic"}, \code{"looic"}, or
  \code{"dic"}. Must be named in full. PSIS-LOO requires the optional
  \pkg{loo} package. See \code{\link{criteria}} for assumptions and diagnostics.}
}
\value{
An object of class \code{"fitdistrBayes"}. It is a list with the following
principal components:
\describe{
  \item{\code{call}}{The matched fitting call.}
  \item{\code{model}}{The canonical model name, estimated parameter names,
  fixed parameters, sample size, and support information.}
  \item{\code{prior}}{The normalized prior identifier, its displayed label and
  kernel, and the posterior-propriety condition verified before fitting.}
  \item{\code{initialization}}{The automatic or user-supplied initialization
  rule, its classical estimation method, and the center used to initialize
  non-exact algorithms.}
  \item{\code{engine}}{The computational algorithm and the requested chain,
  iteration, warmup, thinning, saved-draw, and seed settings.}
  \item{\code{estimates}}{A named numeric vector of posterior medians, used as
  the default point estimates.}
  \item{\code{summary}}{A data frame containing posterior means and standard
  deviations when certified to exist, medians, equal-tail 95 percent credible
  limits, Monte Carlo standard errors, rank-normalized split/folded
  \eqn{\widehat R}, bulk and tail effective sample sizes, and moment-existence
  indicators.}
  \item{\code{moment_status}}{A data frame recording whether the posterior
  mean and variance of each parameter are known to exist and explaining the
  corresponding analytical condition.}
  \item{\code{diagnostics}}{A list containing the overall convergence
  decision, thresholds, extrema of \eqn{\widehat R} and effective sample
  sizes, acceptance information when applicable, and explanatory messages.
  A non-moving MCMC chain sets \eqn{\widehat R} to \code{Inf} and ESS to
  zero for the affected parameter; it cannot certify convergence.
  Independent exact draws do not require MCMC convergence diagnostics.}
  \item{\code{capabilities}}{Logical indicators describing whether prediction
  and pointwise log-likelihood evaluation are available.}
  \item{\code{chains}}{A list of posterior-draw matrices, one matrix per
  chain.}
  \item{\code{draws}}{A data frame combining all posterior draws. Columns
  \code{.chain}, \code{.iteration}, and \code{.draw} identify each draw and
  the remaining columns contain parameter values.}
  \item{\code{data}, \code{omitted}, \code{control}}{The analyzed data, the
  original indices of omitted missing observations, and the validated computational
  controls.}
  \item{\code{criteria}}{Present only when requested: an object of class
  \code{"fitdistrBayes_criteria"}. Computed after sampling, without
  changing the chains. Criteria can also be computed later from the fit.}
}
Density and random-generation closures are retained internally by default to
support \code{log_lik()} and \code{predict()}; they are omitted when
\code{control = list(store_callables = FALSE)}.
}
\description{
Fits a univariate parametric distribution under a registered objective
Bayesian prior or a user-supplied model and prior. For built-in models,
data-dependent posterior-propriety checks
are executed before any numerical integration or simulation. Posterior means,
variances, and their Monte Carlo errors are reported only when the
corresponding moments have been certified to exist.
}
\details{
Recognized distributions are beta, Cauchy, chi-squared, exponential,
Exponential-Logarithmic, Frechet, Gamma, geometric, Gumbel, lognormal,
logistic, Lomax, Nakagami-m, negative binomial, Normal, Poisson, Rician,
Student t, Weibull, and weighted Lindley. Available priors depend on the model and parameter of interest;
unsupported or known-improper model-prior combinations stop with an
informative error.

Use \code{fitdistrBayes_routes()} to obtain the machine-readable catalogue of
all enabled combinations, their estimated and required fixed parameters,
computational engines, and concise propriety conditions. The fitting function
performs the authoritative sample-dependent check.

For the Exponential-Logarithmic model, \code{"reference"} and
\code{"reference-theta"} select the ordering in which theta is the parameter
of interest; \code{"reference-rate"} selects the rate parameter. The Lomax
reference posterior and Nakagami-m MDI posterior are disabled because they
are improper. Rician currently has only the proper joint-Jeffreys route. For
weighted Lindley, \code{"reference"} is the one-group reference prior and is
equal to the Fisher-information Jeffreys prior; \code{"reference-lambda"} and
\code{"reference-phi"} select the exact ordered reference prior for the stated
scalar target. Its MDI posterior is improper and therefore disabled.

For numerical routes, the automatic center and any user overrides are
recorded in \code{fit$initialization}. Each MCMC chain is dispersed from this
center on the unconstrained scale using \code{control$init_jitter}. For
conditionally exact Gamma, Frechet, Nakagami-m, and Weibull algorithms, only the marginal shape
parameter requires initialization; rate or scale is drawn from its exact
conditional distribution.

The negative-binomial routes estimate \code{mu} and require a known positive
\code{fixed$size}. For the Frechet model, the implemented distribution function
is \eqn{F(x)=\exp\{-\mathit{scale}\,x^{-\mathit{shape}}\}}. Thus, the argument
named \code{scale} is the coefficient in the exponent; the conventional
quantile scale is \eqn{\mathit{scale}^{1/\mathit{shape}}}.

Sampler controls are \code{target_accept}, \code{adapt_interval},
\code{proposal_scale}, \code{init_jitter}, \code{slice_width},
\code{slice_steps}, and \code{max_init_tries}. Diagnostic controls are
\code{rhat_threshold}, \code{ess_threshold}, and \code{warn_convergence}; both
bulk and tail ESS must meet the ESS threshold. The \code{overall} acceptance
column uses post-warmup iterations, with warmup and all-iteration rates stored
separately. Entropy evaluation uses \code{entropy_tol} and the positive integer
\code{entropy_exact_limit}. Set \code{store_callables = FALSE} to omit the
closures needed by \code{predict()} and \code{log_lik()}.

For a custom density function, \code{lower} and \code{upper} define componentwise bounds;
the package supplies the corresponding transforms and Jacobian. Log mode is
inferred only from an explicit \code{log} formal. It can be forced with
\code{density_is_log} and \code{prior_is_log}; an initial ordinary/log
consistency check is made when possible. \code{prior_style} is one of
\code{"auto"}, \code{"scalar"}, or \code{"vector"}. A predictive function can
be supplied as \code{rng}, and \code{rng_validator} may check its returned
support. Posterior propriety and moment existence remain the user's
responsibility for custom targets. Use \code{fitdistrBayes_model()} for an
explicit reusable specification that can additionally select the generic
slice sampler or a user-supplied posterior sampler and carry executable
support, propriety, and moment information.
}
\section{Built-in models, parameters, and priors}{
The table below gives the canonical model name, the parameters returned by the
fit, and every registered objective-prior label. Distribution names and prior
labels are matched without regard to case. Hyphens, spaces, and common aliases
such as \code{"log-normal"}, \code{"negative binomial"}, \code{"student-t"},
\code{"el"}, and \code{"nakagami"} are normalized internally.

\tabular{lll}{
\strong{Model} \tab \strong{Estimated parameters} \tab \strong{Available priors} \cr
\code{beta} \tab \code{shape1, shape2} \tab \code{jeffreys, reference} \cr
\code{cauchy} \tab \code{location, scale} \tab \code{jeffreys, reference, mdi} \cr
\code{chi-squared} \tab \code{df} \tab \code{jeffreys, reference} \cr
\code{exponential} \tab \code{rate} \tab \code{jeffreys, reference, mdi} \cr
\code{exponential-logarithmic} \tab \code{theta, rate} \tab \code{jeffreys, mdi} \cr
 \tab  \tab \code{reference-theta, reference-rate} \cr
\code{frechet} \tab \code{shape, scale} \tab \code{jeffreys, reference} \cr
\code{gamma} \tab \code{shape, rate} \tab \code{jeffreys, first-rule} \cr
 \tab  \tab \code{reference-shape, reference-rate} \cr
\code{geometric} \tab \code{prob} \tab \code{jeffreys, reference, mdi} \cr
\code{gumbel} \tab \code{location, scale} \tab \code{jeffreys, reference, mdi} \cr
\code{lognormal} \tab \code{meanlog, sdlog} \tab \code{jeffreys, reference} \cr
\code{logistic} \tab \code{location, scale} \tab \code{jeffreys, reference, mdi} \cr
\code{lomax} \tab \code{shape, scale} \tab \code{jeffreys} \cr
\code{nakagami-m} \tab \code{shape, spread} \tab \code{jeffreys, reference} \cr
\code{negative binomial} \tab \code{mu} \tab \code{jeffreys, reference, mdi} \cr
\code{normal} \tab \code{mean, sd} \tab \code{jeffreys, reference, mdi} \cr
\code{Poisson} \tab \code{lambda} \tab \code{jeffreys, reference, mdi} \cr
\code{rician} \tab \code{noncentrality, scale} \tab \code{jeffreys} \cr
\code{t} \tab \code{location, scale, df} \tab \code{jeffreys, reference, mdi} \cr
 \tab  \tab \code{independence-jeffreys} \cr
\code{weibull} \tab \code{shape, scale} \tab \code{jeffreys, reference} \cr
\code{weighted lindley} \tab \code{lambda, phi} \tab \code{jeffreys, reference, first-rule} \cr
 \tab  \tab \code{independence-jeffreys} \cr
 \tab  \tab \code{reference-lambda, reference-phi}
}

For Student t, \code{independence-jeffreys} estimates \code{df} and requires
at least two pairwise distinct observations. Tied observations produce an
improper posterior when \code{df} is unrestricted and are rejected before
sampling. The other three priors require \code{fixed = list(df = ...)};
their propriety checks account for the largest number of repeated observations
and the fixed value of \code{df}. The generic label
\code{"reference"} is used only where its parameter ordering is unambiguous.
For Gamma and Exponential-Logarithmic models, the parameter of interest is
therefore stated in the label. Unsupported combinations, including objective
priors known to yield an improper posterior, are rejected before sampling.
Run \code{fitdistrBayes_routes()} for the precise data condition and
computational engine attached to every route.

Beta observations must lie strictly between 0 and 1; chi-squared and Rician
observations are positive; Poisson and geometric observations are
nonnegative integers. Exponential, Gamma, geometric, lognormal, logistic,
Normal, Poisson, and Weibull follow the corresponding base-R
parameterizations. Geometric observations count failures before the first
success. The Gumbel model is the location-scale distribution for maxima. For
Lomax, the support is \eqn{x >= 0}; for Nakagami-m, \code{spread} is
\eqn{E(X^2)}. Negative-binomial \code{size} is known and supplied in
\code{fixed}, whereas \code{mu} is estimated. Weighted Lindley observations
are strictly positive, and its density is
\deqn{f(x\mid\lambda,\phi)=\frac{\lambda^{\phi+1}}
{(\lambda+\phi)\Gamma(\phi)}x^{\phi-1}(1+x)e^{-\lambda x}.}
Its MCMC is carried out in the exactly Fisher-orthogonal mean/shape
parameterization while the returned parameters remain \code{lambda} and
\code{phi}.
}
\section{Posterior propriety and reported moments}{
Most objective priors are improper as prior measures. A fit is attempted only
after the built-in route verifies its sample-dependent posterior-propriety
condition. Typical requirements include positivity, a nonconstant sample, a
minimum sample size, or multiplicity restrictions for location-scale models.
The complete concise conditions are returned by
\code{fitdistrBayes_routes()}, while the authoritative check is performed on
the supplied \code{x} by \code{fitdistrBayes()}.

Posterior propriety does not imply that every posterior moment exists. Thus,
\code{mean} or \code{sd} can legitimately be \code{NA} even when medians,
credible intervals, draws, and diagnostics are available. Inspect
\code{fit$moment_status} and its explanatory \code{note} column before treating
an \code{NA} as a computational failure.
}
\section{Output and diagnostics}{
Printing a fit reports the model, prior, computational engine, initialization,
propriety decision, posterior summaries, and diagnostic status. The full draws
are in \code{fit$draws}; chain-specific matrices are in \code{fit$chains}. The
summary contains posterior means when they exist, medians, central credible
limits, MCSE, rank-normalized split/folded \eqn{\widehat R}, and bulk and tail
effective sample sizes.

Use \code{summary()}, \code{coef()}, \code{confint()}, \code{as.data.frame()},
and \code{plot()} for inspection. Use \code{predict()} for posterior predictive
simulation and \code{log_lik()} for a pointwise log-likelihood matrix. The last
two methods are computed on demand and require retained callables.
}
\section{Complete console tutorial}{
The installed file
\code{system.file("examples", "tutorial_fitdistrBayes_all_models.R",
package = "fitdistrBayes")} contains simulated data for all 20 built-in
families, all 56 enabled model-prior routes, comparisons with
\code{MASS::fitdistr()} where directly available, and a final diagnostic table.
It prints to the console and does not save analysis results. Set
\code{quick_mode <- TRUE} for a short demonstration or \code{FALSE} for the
longer teaching run.
}
\references{
Hosking JRM (1990). ``L-moments: Analysis and estimation of distributions
using linear combinations of order statistics.'' \emph{Journal of the Royal
Statistical Society: Series B}, 52(1), 105--124.
\doi{10.1111/j.2517-6161.1990.tb01775.x}.

Ramos PL, Louzada F, Ramos E, Dey S (2020). ``The Frechet distribution:
Estimation and application---an overview.'' \emph{Journal of Statistics and
Management Systems}, 23(3), 549--578. \doi{10.1080/09720510.2019.1645400}.

Ferreira PH, Ramos E, Ramos PL, et al. (2020). ``Objective Bayesian analysis
for the Lomax distribution.'' \emph{Statistics & Probability Letters}, 159,
108677. \doi{10.1016/j.spl.2019.108677}.

Ramos PL, Louzada F, Ramos E (2018). ``Posterior properties of the
Nakagami-m distribution using noninformative priors and applications in
reliability.'' \emph{IEEE Transactions on Reliability}, 67(1), 105--117.
\doi{10.1109/TR.2017.2778139}.

Moala FA, Achire Quispe E, Ramos PL (2026). ``Objective Bayesian inference
for the Exponential-Logarithmic distribution.'' \emph{Journal of Statistical
Computation and Simulation}, 96(8), 1773--1801.
\doi{10.1080/00949655.2025.2608789}.

Achire E, Ramos E, Ramos PL (2025). ``On the posterior property of the
Rician distribution.'' \emph{Statistics}, 59(1), 167--186.
\doi{10.1080/02331888.2024.2425688}.

Mota AL, Ramos PL, Ferreira PH, Tomazella VLD, Louzada F (2021).
``A reparameterized weighted Lindley distribution: Properties, estimation
and applications.'' \emph{Revista Colombiana de Estadistica}, 44(1), 65--90.
\doi{10.15446/rce.v44n1.86566}.
}
\seealso{
\code{\link{fitcensBayes}} provides a separate interface for independently
right-censored observations using the registered complete-data priors.
}
\examples{
# A first exact-posterior example.
set.seed(1)
x <- rexp(30, rate = 2)
fit <- fitdistrBayes(x, "exponential", "jeffreys",
                     iter = 400, warmup = 100, chains = 2, seed = 2)
coef(fit)
confint(fit)
fitdistrBayes_routes("exponential")

# Group different objective priors for the same Gamma data.
set.seed(3)
x_gamma <- rgamma(40, shape = 2.5, rate = 1.3)
gamma_fits <- list(
  Jeffreys = fitdistrBayes(x_gamma, "gamma", "jeffreys",
    iter = 400, warmup = 100, chains = 2, seed = 4,
    control = list(warn_convergence = FALSE)),
  Reference_shape = fitdistrBayes(x_gamma, "gamma", "reference-shape",
    iter = 400, warmup = 100, chains = 2, seed = 5,
    control = list(warn_convergence = FALSE)),
  Reference_rate = fitdistrBayes(x_gamma, "gamma", "reference-rate",
    iter = 400, warmup = 100, chains = 2, seed = 6,
    control = list(warn_convergence = FALSE))
)
lapply(gamma_fits, coef)

# A fixed parameter: negative-binomial size is known and mu is estimated.
set.seed(7)
x_nb <- rnbinom(40, size = 5, mu = 8)
fit_nb <- fitdistrBayes(x_nb, "negative binomial", "reference",
  fixed = list(size = 5), iter = 400, warmup = 100,
  chains = 2, seed = 8)
coef(fit_nb)

# Student t: fixed df versus estimated df.
set.seed(9)
x_t <- 1 + 2 * rt(50, df = 7)
fit_t_fixed <- fitdistrBayes(x_t, "t", "jeffreys",
  fixed = list(df = 7), iter = 400, warmup = 100,
  chains = 2, seed = 10, control = list(warn_convergence = FALSE))
fit_t_unknown <- fitdistrBayes(x_t, "t", "independence-jeffreys",
  iter = 400, warmup = 100, chains = 2, seed = 11,
  control = list(warn_convergence = FALSE))
coef(fit_t_fixed)
coef(fit_t_unknown)

# One short fit for every built-in family. The installed tutorial additionally
# runs every available prior for each family.
quick <- list(warn_convergence = FALSE)

set.seed(101)
fit_beta <- fitdistrBayes(rbeta(30, 2, 5), "beta", "jeffreys",
  iter = 120, warmup = 40, chains = 2, seed = 102, control = quick)
fit_cauchy <- fitdistrBayes(rcauchy(40, 1, 2), "cauchy", "reference",
  iter = 120, warmup = 40, chains = 2, seed = 103, control = quick)
fit_chisq <- fitdistrBayes(rchisq(30, 6), "chi-squared", "jeffreys",
  iter = 120, warmup = 40, chains = 2, seed = 104, control = quick)
fit_exponential <- fitdistrBayes(rexp(30, 1.5), "exponential", "mdi",
  iter = 120, warmup = 40, chains = 2, seed = 105, control = quick)
fit_gamma <- fitdistrBayes(rgamma(30, 2.5, rate = 1.3),
  "gamma", "reference-shape", iter = 120, warmup = 40,
  chains = 2, seed = 106, control = quick)
fit_geometric <- fitdistrBayes(rgeom(30, 0.35), "geometric", "reference",
  iter = 120, warmup = 40, chains = 2, seed = 107, control = quick)
fit_lognormal <- fitdistrBayes(rlnorm(30, 1, 0.6), "lognormal", "jeffreys",
  iter = 120, warmup = 40, chains = 2, seed = 108, control = quick)
fit_logistic <- fitdistrBayes(rlogis(40, 1, 2), "logistic", "mdi",
  iter = 120, warmup = 40, chains = 2, seed = 109, control = quick)
fit_nb <- fitdistrBayes(rnbinom(30, size = 5, mu = 8),
  "negative binomial", "jeffreys", fixed = list(size = 5),
  iter = 120, warmup = 40, chains = 2, seed = 110, control = quick)
fit_normal <- fitdistrBayes(rnorm(30, 3, 2), "normal", "reference",
  iter = 120, warmup = 40, chains = 2, seed = 111, control = quick)
fit_poisson <- fitdistrBayes(rpois(30, 4), "Poisson", "mdi",
  iter = 120, warmup = 40, chains = 2, seed = 112, control = quick)
fit_t <- fitdistrBayes(rt(40, 7), "t", "jeffreys", fixed = list(df = 7),
  iter = 120, warmup = 40, chains = 2, seed = 113, control = quick)
fit_weibull <- fitdistrBayes(rweibull(30, 1.6, 2), "weibull", "reference",
  iter = 120, warmup = 40, chains = 2, seed = 114, control = quick)

x_gumbel <- 1 - 2 * log(-log(runif(40)))
fit_gumbel <- fitdistrBayes(x_gumbel, "gumbel", "reference",
  iter = 120, warmup = 40, chains = 2, seed = 115, control = quick)

x_frechet <- (4 / rexp(40))^(1 / 2.5)
fit_frechet <- fitdistrBayes(x_frechet, "frechet", "jeffreys",
  iter = 120, warmup = 40, chains = 2, seed = 116, control = quick)

x_lomax <- 2 * expm1(-log(runif(40)) / 3)
fit_lomax <- fitdistrBayes(x_lomax, "lomax", "jeffreys",
  iter = 120, warmup = 40, chains = 2, seed = 117, control = quick)

x_nakagami <- sqrt(rgamma(40, 2.5, rate = 2.5 / 4))
fit_nakagami <- fitdistrBayes(x_nakagami, "nakagami-m", "reference",
  iter = 120, warmup = 40, chains = 2, seed = 118, control = quick)

theta_el <- 0.4
rate_el <- 1.3
u_el <- runif(40)
x_el <- -(log(-expm1((1 - u_el) * log(theta_el))) -
  log1p(-theta_el)) / rate_el
fit_el <- fitdistrBayes(x_el, "exponential-logarithmic", "reference-rate",
  iter = 120, warmup = 40, chains = 2, seed = 119, control = quick)

x_rician <- sqrt(rnorm(40, 5, 2)^2 + rnorm(40, 0, 2)^2)
fit_rician <- fitdistrBayes(x_rician, "rician", "jeffreys",
  iter = 120, warmup = 40, chains = 2, seed = 120, control = quick)

lambda_wl <- 2.5
phi_wl <- 0.8
component_wl <- runif(40) < lambda_wl / (lambda_wl + phi_wl)
x_wl <- rgamma(40, shape = phi_wl + as.numeric(!component_wl),
  rate = lambda_wl)
fit_wl <- fitdistrBayes(x_wl, "weighted lindley", "reference-lambda",
  iter = 120, warmup = 40, chains = 2, seed = 121, control = quick)

tutorial <- system.file(
  "examples", "tutorial_fitdistrBayes_all_models.R",
  package = "fitdistrBayes"
)
tutorial
if (interactive()) file.show(tutorial)
}
