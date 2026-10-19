# fitdistrBayes

`fitdistrBayes` 0.6.0 fits common univariate distributions under registered
objective Bayesian priors through a deliberately small interface:

```r
fitdistrBayes(x, distr, prior)
```

For a built-in route, the package checks the sampling support and the known
posterior-propriety condition before computation. It also records whether the
posterior mean and variance of each parameter are mathematically certified to
exist, so a finite chain is not used to justify a divergent moment.

## Installation

Install the source archive supplied with this research release and load it:

```r
install.packages("fitdistrBayes_0.6.0.tar.gz", repos = NULL,
                 type = "source")
library(fitdistrBayes)
```

## Dhillon (new in 0.6.0)

The complete-data and censored interfaces accept `"dhillon"`, with survival
`1/(1 + theta*t^beta)` and positive parameters `beta, theta`. The aliases
`"log-logistic"` and `"loglogistic"` use exactly this parameterization.
Jeffreys and `"reference-beta"` coincide; `"reference-theta"` is a distinct
ordered prior. `"reference"` means `"reference-beta"`. MDI is rejected.

```r
set.seed(19)
x <- exp((rlogis(100) - log(2)) / 3)
fit <- fitdistrBayes(x, "dhillon", "reference-beta", seed = 20)
observed <- pmin(x, 1)
status <- as.integer(x <= 1)
fit_cens <- fitcensBayes(observed, status, "dhillon", "reference-theta", seed = 21)
fit$moment_status
fit$predictive_moments
help("Dhillon")
```

Both observed-likelihood and augmentation methods are available. Exact
observed-data propriety and parameter-moment checks run before sampling.
Predictive lifetime means are infinite: use quantiles and survival probabilities.
The console tutorial is `inst/examples/tutorial_dhillon.R`; the censored-data
derivation is in `inst/doc/dhillon-theory.txt` (installed as `doc/`).

## Optional model comparison

Both fitting functions default to `criteria = FALSE`: no additional likelihood
matrix or information criterion is computed. Existing R-hat/ESS diagnostics
are independent of this option. Use `criteria = TRUE` to compute all criteria
after sampling, or choose `criteria = c("waic", "dic")` without extra dependencies.
PSIS-LOO requires the optional `loo` package. The fitting chains are unchanged.

```r
set.seed(20)
x <- rexp(60)
fit <- fitdistrBayes(x, "exponential", "reference", seed = 21,
                    criteria = c("waic", "dic"))
fit$criteria
WAIC(fit)                       # May also be computed later
# install.packages("loo")       # Once, if PSIS-LOO is desired
# LOOIC(fit)
```

The same functions accept `fitcensBayes` fits and use log survival for censored
observations, never imputed-data densities. DIC is unavailable if necessary
posterior means are not certified finite. LOO is withheld if a training posterior
is not certified proper. Inspect diagnostic warnings before model comparison.
`compare_models()` requires the same observations, censoring status, and scale.
See `help("criteria")` for definitions, output components, and limitations.

## Independent right censoring

The separate function `fitcensBayes(x, status, distr, prior)` uses `status = 1`
for exact observations and `status = 0` for the strict event `T > x`.
For discrete distributions this distinction is important: a record meaning
`T >= 3` must be supplied as `x = 2, status = 0`.

```r
set.seed(20)
lifetime <- rexp(80, rate = 0.7)
censoring <- rexp(80, rate = 0.3)
x <- pmin(lifetime, censoring)
status <- as.integer(lifetime <= censoring)
fit <- fitcensBayes(x, status, "exponential", "reference", seed = 21)
summary(fit)
confint(fit)
predict(fit, type = "survival", times = c(1, 2, 3), draws = 5)
log_lik(fit, draws = 5)
```

The catalogue `fitcensBayes_models()` lists 59 registered model-prior routes
for 21 distributions. The priors are inherited from the **complete-data model**;
they are not claimed to be Jeffreys or reference priors rederived for a
particular censoring design. Posterior propriety is certified by sufficient
conditions on the exact-event subset, with an additional analytic Exponential
case and an exact observed-data condition for Dhillon. An uncertified case is rejected before sampling. Imputed observations
are never used to satisfy these conditions.

`method = "auto"` uses the observed likelihood when censoring is present;
`method = "augmentation"` alternates truncated lifetime simulation with
parameter updates. Without censoring, `method = "auto"` uses the original
complete-data sampler. Posterior medians, credible intervals, rank-based
diagnostics and justified moment summaries are returned. Heavy censoring may
require substantially longer chains. A diagnostic pass is not a convergence
proof. Densities or priors defined by the user remain available through the
original `fitdistrBayes()` interface, but not the new censored interface.

See `help("fitcensBayes")` and the installed English tutorial:

```r
system.file("examples", "tutorial_fitcensBayes.R", package = "fitdistrBayes")
```

## Quick start

An exact-posterior example:

```r
set.seed(10)
x <- rexp(40, rate = 2)
fit_exp <- fitdistrBayes(x, "exponential", "jeffreys", seed = 11)
fit_exp
confint(fit_exp)
```

An MCMC example with a parameter-specific reference prior:

```r
set.seed(20)
x <- rgamma(50, shape = 2.5, rate = 1.3)
fit_gamma <- fitdistrBayes(
  x, "gamma", "reference-shape",
  iter = 6000, warmup = 1000, chains = 4, seed = 21
)
summary(fit_gamma)
plot(fit_gamma, type = "trace")
plot(fit_gamma, type = "acf")
```

`fit$summary` contains posterior medians, equal-tail intervals, R-hat, bulk and
tail ESS, and MCSE when the relevant moments exist. Acceptance rates labelled
`overall` use post-warmup iterations; warmup and all-iteration rates are stored
separately.

Posterior prediction and pointwise log likelihood are computed on demand:

```r
yrep <- predict(fit_gamma, draws = 500, size = length(x), seed = 22)
ll <- log_lik(fit_gamma, draws = 500, seed = 23)
dim(yrep)
dim(ll)
```

## Supported routes

The package has 59 enabled model--prior combinations for 21 distributions.
For Student-t with unknown degrees of freedom, the independence Jeffreys prior
requires at least two pairwise distinct observations; samples with ties are
rejected because the posterior is improper. Fixed-df Student-t models use
their separate, multiplicity-dependent conditions.
The public machine-readable catalogue is the concise way to inspect them:

```r
fitdistrBayes_routes()
fitdistrBayes_routes("gamma")
fitdistrBayes_routes("t")
fitdistrBayes_routes("weighted lindley")
```

The catalogue reports the estimated parameters, any required fixed parameter,
the computational engine, and a concise propriety condition. The fitting
function performs the authoritative sample-dependent check. In particular,
the negative-binomial routes estimate `mu` with known positive `size`:

```r
x <- rnbinom(40, size = 5, mu = 3)
fit_nb <- fitdistrBayes(
  x, "negative binomial", "reference", fixed = list(size = 5), seed = 30
)
```

## Automatic initialization

When `start = NULL` and numerical initialization is required, the package uses
classical estimators: ordinary moments where they exist, quantile matching for
Cauchy and heavy-tailed Student-t cases, and closed-form L-moments for Weibull,
Frechet, and Lomax. The Exponential-Logarithmic model uses a stable scalar
moment equation. Exact independent posterior simulation does not require a
starting point. The selected rule and center are stored in
`fit$initialization`.

For weighted Lindley data, the package uses the model's closed-form likelihood
estimator (with a numerical maximum-likelihood fallback) and samples in the
exactly Fisher-orthogonal mean/shape parameterization. The public output remains
on the original `lambda, phi` scale. The one-group `reference` prior equals
Fisher-information Jeffreys; `reference-lambda` and `reference-phi` select the
two exact ordered reference priors.

For the implemented Frechet parameterization,
`F(x) = exp(-scale * x^(-shape))`; hence the argument named `scale` is the
positive coefficient in the exponent, while the conventional quantile scale
is `scale^(1 / shape)`.

## Custom densities

A custom density requires named starting values and a prior function. Log mode
is inferred only from an explicit `log` argument; the presence of `...` alone
does not imply that `log = TRUE` is honored. The contract can be set explicitly
with `control$density_is_log` and `control$prior_is_log`. Use
`control$prior_style = "scalar"` or `"vector"` when automatic prior calling is
ambiguous. Bounds are supplied through `control$lower` and `control$upper`; the
package constructs the componentwise transform and Jacobian. A predictive RNG
and its optional support validator can be supplied through `control$rng` and
`control$rng_validator`.

Posterior propriety and posterior-moment existence remain the user's
responsibility for every custom target.

## Scope and limitations

The package handles complete iid univariate samples and independent,
noninformative right censoring through `fitcensBayes()`. Optional WAIC, DIC,
and PSIS-LOO calculations use the observed-data likelihood. It refuses
routes known to be improper, including the Lomax independent
Jeffreys/reference posterior and the Nakagami-m MDI posterior, and it does not
substitute an approximation when full-model propriety has not been established.
Numerical overflow or underflow stops with an explicit error; values are not
silently clipped into the floating-point range.

## Teaching example, citation, and support

A console-oriented walkthrough is installed at `inst/examples/teaching.R`.
For the complete mathematical catalogue and computational validation, cite the
accompanying manuscript, *fitdistrBayes: Objective Bayesian Distribution
Fitting in R*. Until a public issue tracker is announced, reproducible problem
reports can be sent to the maintainer address in `DESCRIPTION`.

This research version is licensed under GPL-3.
