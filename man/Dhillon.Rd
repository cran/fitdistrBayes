\name{Dhillon}
\alias{Dhillon}
\alias{dhillon}
\title{Dhillon Distribution: Priors, Censoring and Moment Conditions}
\description{
The Dhillon model is available through \code{fitdistrBayes(x, "dhillon", prior)}
and \code{fitcensBayes(x, status, "dhillon", prior)}. The aliases
\code{"log-logistic"} and \code{"loglogistic"} use the same parameterization.
}
\details{
\strong{Parameterization.} Both \code{beta} and \code{theta} are positive:
\deqn{f(t\mid\beta,\theta)=\frac{\beta\theta t^{\beta-1}}{(1+\theta t^\beta)^2},
\qquad S(t\mid\beta,\theta)=\frac{1}{1+\theta t^\beta},\quad t>0.}
The conventional log-logistic scale is \eqn{\theta^{-1/\beta}}, not
\eqn{\theta}. Equivalently, \eqn{\log T} is logistic with location
\eqn{-\log(\theta)/\beta} and scale \eqn{1/\beta}. Simulate using
\code{exp((rlogis(n) - log(theta)) / beta)}.

\strong{Priors.} All densities below are with respect to
\eqn{d\beta\,d\theta} in this parameterization:
\deqn{\pi_J(\beta,\theta)=\pi_{R,\beta}(\beta,\theta)
\propto\frac{1}{\beta\theta},\qquad
\pi_{R,\theta}(\beta,\theta)\propto
\frac{1}{\beta\theta\sqrt{\pi^2+3+3\log^2\theta}}.}
Select \code{"jeffreys"}, \code{"reference-beta"}, or
\code{"reference-theta"}. The alias \code{"reference"} selects
\code{"reference-beta"}. The two ordered reference priors must not be
confused. The MDI kernel is \eqn{\beta\theta^{1/\beta}} and yields an
improper posterior for every finite complete or right-censored sample;
\code{"mdi"} is therefore rejected before sampling. Fixed beta or theta
submodels are not registered in this version.

\strong{Right censoring.} Priors are those of the complete-data model,
not priors rederived for a particular censoring design. Censoring must be
independent and non-informative. With event indicators \eqn{\delta_i},
the observed likelihood is \eqn{\prod_i f(x_i)^{\delta_i}S(x_i)^{1-\delta_i}}.
Exact times must be positive; censoring limits may be zero. Zero limits
have survival one and add no information.

\strong{Posterior propriety.} For complete data, at least two distinct
positive observations are necessary and sufficient under each registered
prior. For right-censored data the exact condition is: at least one event,
and either two distinct event times or a censoring limit strictly greater
than the common event time. Thus a sample with one event and a later
censoring limit can be fitted. All-censored samples, and samples with
one common event time and only earlier or equal censoring limits, are rejected.

\strong{Parameter moments.} Whenever the posterior is proper, all positive
moments of \eqn{\beta} exist. Remove zero censoring limits and let \eqn{N}
be the number of remaining observations. Define
\deqn{A_q(z)=\sum_{i:\delta_i=1}|z+\log x_i|
+\sum_{i:\delta_i=0}\max(z+\log x_i,0)-qz.}
For each registered prior, \eqn{E(\theta^q\mid x,\delta)} is finite if and
only if \eqn{0<q<N} and \eqn{\min_z A_q(z)>0}. The minimum is attained at
a knot \eqn{-\log x_i}, obtainable as a weighted quantile in
\eqn{O(N\log N)} operations. The package checks orders one and two before
reporting means, SDs, and mean Monte Carlo errors. A sign too close to a
floating-point cancellation boundary is not certified and is labelled
\code{NA}. Known divergence is labelled \code{FALSE}. DIC is withheld
when either required posterior mean is not certified finite.

These theta-moment conditions depend on the measurement unit: rescaling
times by \eqn{c} changes theta to \eqn{\theta c^{-\beta}}. The reference
prior for theta is also tied to this stated parameterization and compact-set
construction; it should not be silently replaced after a change of units.

\strong{Predictive moments.} All positive predictive lifetime moments are
infinite, because the posterior has positive mass at arbitrarily small
positive beta. The predictive mean is infinite; variance and SD are
undefined. Finite simulated averages are not estimates of a finite
predictive mean. Use predictive medians, quantiles, and survival probabilities.
This information is recorded in \code{fit$predictive_moments}.

\strong{Computation.} Automatic starts use closed-form moments of log exact
times: \eqn{\beta_0=\pi/(\sqrt{3}s_{\log x})} and
\eqn{\log\theta_0=-\beta_0\overline{\log x}}, where the variance divisor is
the number of exact events. This rule does not require finite raw lifetime
moments. If the exact events have no variation but the censored posterior
is proper, positive observed times provide a starting heuristic only.
Observed-posterior optimization refines automatic starts. A supplied start
is respected. Sampling uses preconditioned slice updates in centered
log coordinates, with the correct Jacobian and geometry fixed during
retained simulation. The censored interface also supports
\code{method="augmentation"}, using exact inverse-survival truncated draws.
Latent log times are retained internally, so a lifetime beyond the natural
floating-point range does not prevent parameter sampling. Predictions requested
on the original time scale still stop explicitly if a draw is not representable;
they are not clipped. Survival probabilities avoid this conversion.
The posterior target is evaluated without exponentiating log theta. If an
initial or saved parameter value itself is outside the numerical range,
fitting stops explicitly rather than silently bounding the parameter space.
WAIC and PSIS-LOO always use the observed likelihood. Each LOO training
posterior is checked separately. The extended censored-data conditions are
derived in the installed file \code{doc/dhillon-theory.txt}.
}
\value{
Fitting returns a \code{"fitdistrBayes"} or \code{"fitcensBayes"} list,
respectively, with parameter names \code{beta, theta}, posterior draws,
medians, intervals, convergence diagnostics and \code{moment_status}.
The additional list \code{predictive_moments} contains logical
\code{mean_exists}, \code{variance_exists} and an explanatory \code{note}.
This help topic itself does not define a separate fitting function.
}
\references{
Ramos PL, Achire Quispe E, de Oliveira RP, Achcar JA (2026).
``Objective Bayesian inference for the Dhillon distribution.''
Unpublished manuscript.
}
\seealso{\code{\link{fitdistrBayes}}, \code{\link{fitcensBayes}},
\code{\link{fitdistrBayes_routes}}, \code{\link{criteria}}}
\examples{
set.seed(19)
x <- exp((rlogis(30) - log(2)) / 3)
# Short runs demonstrate the interface, not convergence.
fit <- fitdistrBayes(x, "dhillon", "reference-beta",
  iter = 120, warmup = 40, chains = 2, seed = 20,
  control = list(warn_convergence = FALSE))
coef(fit)
fit$moment_status
fit$predictive_moments

observed <- pmin(x, 1)
status <- as.integer(x <= 1)
fit_cens <- fitcensBayes(observed, status, "dhillon", "reference-theta",
  iter = 120, warmup = 40, chains = 2, seed = 21,
  control = list(warn_convergence = FALSE))
coef(fit_cens)
predict(fit_cens, type = "survival", times = c(0, 1, 2), draws = 3)
}
\keyword{models}
