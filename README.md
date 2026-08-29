# fitdistrBayes

`fitdistrBayes` 0.2.2 fits common univariate distributions under registered
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
install.packages("fitdistrBayes_0.2.2.tar.gz", repos = NULL,
                 type = "source")
library(fitdistrBayes)
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

Version 0.2.2 has 56 enabled model--prior combinations for 20 distributions.
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

The current package is for complete iid univariate samples. It does not
implement censoring, truncation, observation weights, regression, hierarchical
models, or automatic comparison among candidate distributions. It refuses
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
