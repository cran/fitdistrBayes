library("fitdistrBayes")

set.seed(270719)
dat <- list(
  beta = rbeta(25, 2, 4),
  cauchy = rcauchy(25, 0, 1),
  chisq = rchisq(25, 5),
  exponential = rexp(25, 2),
  el = NULL,
  frechet = (4 / rexp(25))^(1 / 2.5),
  gamma = rgamma(25, 3, rate = 2),
  geometric = rgeom(25, 0.4),
  gumbel = 1 - 2 * log(-log(runif(25))),
  lognormal = rlnorm(25, 0.2, 0.5),
  logistic = rlogis(25, 0, 1),
  lomax = 3 * (runif(25)^(-1 / 3) - 1),
  nakagami = sqrt(rgamma(25, 2.5, rate = 2.5 / 4)),
  nbinom = rnbinom(25, size = 5, mu = 4),
  normal = rnorm(25),
  poisson = rpois(25, 3),
  rician = sqrt(rnorm(25, 5, 2)^2 + rnorm(25, 0, 2)^2),
  t = rt(30, 7),
  weibull = rweibull(25, 2, 3)
)
theta_el <- 0.4
rate_el <- 1.3
u_el <- runif(25)
dat$el <- -(
  log(-expm1((1 - u_el) * log(theta_el))) - log1p(-theta_el)
) / rate_el

registry <- fitdistrBayes_routes()
if (nrow(registry) != 50L || anyDuplicated(registry)) {
  stop("The public route catalogue must contain 50 unique entries.",
       call. = FALSE)
}

data_key <- c(
  beta = "beta", cauchy = "cauchy", `chi-squared` = "chisq",
  exponential = "exponential", `exponential-logarithmic` = "el",
  frechet = "frechet", gamma = "gamma", geometric = "geometric",
  gumbel = "gumbel", lognormal = "lognormal", logistic = "logistic",
  lomax = "lomax", `nakagami-m` = "nakagami",
  `negative binomial` = "nbinom", normal = "normal", Poisson = "poisson",
  rician = "rician", t = "t", weibull = "weibull"
)

results <- vector("list", nrow(registry))
for (position in seq_len(nrow(registry))) {
  route <- registry[position, ]
  model <- route$model
  prior <- route$prior
  fixed <- if (identical(route$required_fixed, "size")) list(size = 5) else
    if (identical(route$required_fixed, "df")) list(df = 7) else NULL
  fit <- fitdistrBayes(
    dat[[data_key[[model]]]], model, prior, fixed = fixed,
    iter = 240, warmup = 120, chains = 2, seed = 1000 + position,
    control = list(rhat_threshold = 2, ess_threshold = 1,
                   warn_convergence = FALSE)
  )
  if (!inherits(fit, "fitdistrBayes") ||
      any(!is.finite(as.matrix(fit$draws[, fit$model$parameters,
                                         drop = FALSE]))) ||
      !identical(fit$engine$algorithm, route$engine)) {
    stop(sprintf("Registry test failed for %s / %s.", model, prior),
         call. = FALSE)
  }
  results[[position]] <- data.frame(
    model = model, prior = prior, engine = fit$engine$algorithm,
    draws = nrow(fit$draws), row.names = NULL
  )
}

results <- do.call(rbind, results)
print(results, row.names = FALSE)
cat(sprintf("All %d enabled model-prior registry entries executed successfully.\n",
            nrow(results)))
