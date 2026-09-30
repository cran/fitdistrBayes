# fitdistrBayes 0.5.1

- Fix WAIC portability on macOS ARM and R builds without extended-precision
  long double arithmetic. Center and scale log-likelihood contributions before
  computing their dispersion, avoiding overflow in intermediate sums.
- Use the same stable dispersion calculation for paired information-criterion
  standard errors. Return an explicitly unavailable WAIC when finite inputs
  produce pointwise or aggregate scores outside the representable range.
- Retain the extreme-value regression that exposed the problem, and add tests
  for double-only arithmetic, large constant offsets, genuine overflow and
  agreement with ordinary-range WAIC calculations. Priors and posterior
  sampling code are unchanged from 0.5.0.

# fitdistrBayes 0.5.0

- Harden the pre-release numerical audit: use exceptional-range log-scale
  fallbacks for the weighted Lindley mean inversion and censored transforms;
  avoid rounding a positive beta shape to zero; stabilize SD and mean ESS
  under changes of units; flag non-moving MCMC chains explicitly.
- Handle sure censoring events at zero analytically in PSIS-LOO instead of
  fitting degenerate Pareto tails. Guard aggregate criterion/SE overflow,
  withhold DIC when the likelihood at the posterior mean is undefined, and
  retain completed posterior draws if optional post-processing fails.
- Add permanent regression tests for these audit findings. Priors and
  sampling algorithms are unchanged; ordinary-range regression outputs
  remain identical for all 56 complete-data model-prior routes.
- Reject unused built-in fitting arguments and non-vector/complex data rather
  than silently ignoring arguments or flattening arrays. Do not clip boundary
  Beta draws in the complete-data geometric posterior sampler.
- Add opt-in criteria = FALSE to both fitting interfaces. TRUE computes
  WAIC, PSIS-LOO, and DIC after sampling; character vectors select a subset.
  Add criteria(), WAIC(), LOOIC(), DIC(), and compare_models(). Default
  fitting does not evaluate additional likelihoods or load 'loo'.
- Use the observed likelihood under censoring, guard leave-one-out posterior
  propriety, respect certified parameter moments for DIC, and report
  pointwise diagnostics and paired comparison uncertainty. PSIS-LOO uses
  the optional 'loo' package; WAIC and DIC need no additional dependency.
- Add the separate function fitcensBayes(x, status, distr, prior) for
  independent right censoring, with status 1 for exact events and status 0
  for strictly greater lifetimes, including discrete distributions.
- Support the 20 registered distributions and 56 complete-data model-prior
  routes. These priors are reused, not rederived for a censoring design.
  Check sufficient posterior-propriety conditions on the exact-event subset
  before sampling; handle the Exponential posterior analytically.
- Add exact censored Gamma and Beta posteriors where available, collapsed
  Weibull sampling, observed-likelihood slice sampling, and optional data
  augmentation with truncated lifetimes. Improve numerical evaluation of
  extreme Rician tails in the new censored interface.
- Add the fitcensBayes class, posterior summaries and diagnostics, survival
  prediction, imputation, and observed pointwise log-likelihood extraction.
  Only certified posterior moments are reported.
- Add fitcensBayes_models(), help pages, examples and regression tests.
- Preserve sampling algorithms and ordinary-range regression results; the new
  criteria argument is named-only to preserve existing positional calls.

# fitdistrBayes 0.2.3

- Reject tied observations for the Student-t model with unknown degrees of
  freedom under the independence Jeffreys prior. Such samples have an improper
  posterior because the prior gives positive mass to degrees of freedom for
  which the location-scale posterior is not integrable near zero scale.
- Document the pairwise-distinct sample requirement in the help page and public
  catalogue, and add regression tests for tied and distinct observations.
- The fixed-degree-of-freedom Student-t routes retain their existing,
  multiplicity-dependent propriety conditions. Sampling algorithms and the
  remaining model-prior routes are unchanged.

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
