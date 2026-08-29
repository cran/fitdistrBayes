# fitdistrBayes 0.2.2

- Added the weighted Lindley distribution with Fisher-information Jeffreys,
  one-group reference, Jeffreys' first-rule, independence-Jeffreys, and exact
  ordered reference priors for lambda and phi. All six routes enforce the proved condition that the
  sample is positive, nonconstant, and has at least two observations.
- Disabled weighted Lindley MDI fitting because the corresponding posterior
  is improper for every finite sample.
- Added the closed-form weighted Lindley likelihood initialization with a
  numerical maximum-likelihood fallback, stable log-prior evaluation, a
  Fisher-orthogonal mean/shape Metropolis parameterization, mixture-based posterior prediction,
  and pointwise log-likelihood evaluation.
- Expanded the public catalogue to 20 distributions and 56 enabled
  model--prior routes, with analytic, initialization, predictive, guard, and
  registry tests for the new family.

# fitdistrBayes 0.2.1

- Added `fitdistrBayes_model()`, a public model specification for combining a
  user-defined density and prior with parameter support, starting values,
  validation, posterior-propriety and moment declarations, posterior
  prediction, and an adaptive Metropolis, slice, or user-supplied sampler.
- Extended `fitdistrBayes()` so the distribution argument can be a
  `fitdistrBayes_model` object. The fitted object records which extension
  components and computational engine were used.
- Added executable tests for a Laplace model with a proper non-objective prior,
  a slice-sampled Exponential model, and an exact Poisson--Gamma posterior
  simulator.

# fitdistrBayes 0.2.0

- Expanded the help pages to document the class, structure, and statistical
  meaning of every exported function and S3 method return value.
- Replaced the non-executable tutorial wrapper with an interactive file-view
  example and left all fitting examples unwrapped after verifying that the
  complete help example runs in less than five seconds.
- Removed explicit access to `.GlobalEnv` and all superassignment operators
  from the package code and validation scripts.
- Added Gumbel, Frechet, Lomax, Nakagami-m, Exponential-Logarithmic, and
  Rician models with documented propriety conditions enforced on the supported
  data domain.
- Added closed-form moment starts for Gumbel, Nakagami-m, and Rician,
  closed-form L-moment starts for Frechet and Lomax, and a stable
  one-dimensional moment start for the Exponential-Logarithmic model.
- Added marginal/conditional samplers for Frechet and Nakagami-m, stable
  transformed Metropolis samplers for the other new models, and an accelerated
  Gauss-Laguerre evaluation of the Rician Jeffreys factor.
- Reparameterized the Lomax Metropolis sampler by shape and scale/shape to
  reduce posterior correlation while retaining the documented output scale.
- Added explicit failures for the improper Lomax reference posterior and
  Nakagami-m MDI posterior, as well as for unproved Rician reference/MDI
  routes.
- Expanded the analytic, numerical, initialization, route, and stress tests to
  19 built-in distributions and 50 enabled model-prior routes.
- Added `fitdistrBayes_routes()` as a machine-readable catalogue of all 50
  enabled routes and made the route test execute directly from this catalogue.
- Expanded `help("fitdistrBayes")` with all built-in models, parameterizations,
  priors, fixed-parameter requirements, propriety guidance, and one example for
  every family; the English-language 50-route console tutorial is installed as
  `examples/tutorial_fitdistrBayes_all_models.R`.
- Tightened the Nakagami-m propriety check to reject constant positive samples.
- Made the custom density/prior log-mode contract explicit and reject undefined
  or positive-infinite log values; custom RNG outputs can now be validated.
- Reported post-warmup, warmup, and all-iteration acceptance rates separately,
  included tail ESS in the convergence decision, and corrected constant-chain
  and per-chain ACF diagnostics.
- Removed silent floating-point clipping: unrepresentable natural-scale draws
  now produce transparent numerical-range errors.

# fitdistrBayes 0.1.2

- Added classical automatic initialization, including the closed-form
  Weibull L-moment estimator.
- Made non-constant-sample validation invariant to translation.
- Replaced costly truncated negative-binomial entropy sums with an exact
  numerical integral for MDI fitting.
- Stabilized natural-scale Gamma and Weibull draws at floating-point limits.
- Preserved the matrix shape of `log_lik()` for one-observation samples.
- Scaled Metropolis location proposals to the model's initial scale, fixing
  poor mixing for very small- and very large-scale data.
- Added broad stress, entropy, initialization, edge-case, and difficult-chain
  regression tests.
