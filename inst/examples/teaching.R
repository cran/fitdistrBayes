# fitdistrBayes 0.2.2: newer built-in models
# Run one section at a time. Results are printed in the console only.

library(fitdistrBayes)

shown <- c("parameter", "median", "q2.5", "q97.5", "rhat", "ess_bulk")

# Gumbel ----------------------------------------------------------------------
set.seed(101)
x_gumbel <- 1 - 2 * log(-log(runif(80)))
gumbel_jeffreys <- fitdistrBayes(x_gumbel, "gumbel", "jeffreys", seed = 102)
gumbel_reference <- fitdistrBayes(x_gumbel, "gumbel", "reference", seed = 103)
gumbel_mdi <- fitdistrBayes(x_gumbel, "gumbel", "mdi", seed = 104)
list(
  Jeffreys = gumbel_jeffreys$summary[, shown],
  Reference = gumbel_reference$summary[, shown],
  MDI = gumbel_mdi$summary[, shown]
)

# Frechet: F(x)=exp(-scale*x^(-shape)) ---------------------------------------
set.seed(201)
x_frechet <- (4 / rexp(80))^(1 / 2.5)
frechet_jeffreys <- fitdistrBayes(x_frechet, "frechet", "jeffreys", seed = 202)
frechet_reference <- fitdistrBayes(x_frechet, "frechet", "reference", seed = 203)
list(
  Jeffreys = frechet_jeffreys$summary[, shown],
  Reference = frechet_reference$summary[, shown]
)

# Lomax -----------------------------------------------------------------------
set.seed(301)
x_lomax <- 3 * expm1(-log(runif(80)) / 3)
lomax_jeffreys <- fitdistrBayes(x_lomax, "lomax", "jeffreys", seed = 302)
list(Jeffreys = lomax_jeffreys$summary[, shown])
# The reference posterior is improper and is deliberately unavailable.

# Nakagami-m: spread=E(X^2) ---------------------------------------------------
set.seed(401)
x_nakagami <- sqrt(rgamma(80, shape = 2.5, rate = 2.5 / 4))
nakagami_jeffreys <- fitdistrBayes(
  x_nakagami, "nakagami-m", "jeffreys", seed = 402
)
nakagami_reference <- fitdistrBayes(
  x_nakagami, "nakagami-m", "reference", seed = 403
)
list(
  Jeffreys = nakagami_jeffreys$summary[, shown],
  Reference = nakagami_reference$summary[, shown]
)
# The MDI posterior is improper and is deliberately unavailable.

# Exponential-Logarithmic -----------------------------------------------------
set.seed(501)
theta <- 0.4
rate <- 1.3
u <- runif(80)
x_el <- -(log(-expm1((1 - u) * log(theta))) - log1p(-theta)) / rate
el_jeffreys <- fitdistrBayes(x_el, "el", "jeffreys", seed = 502)
el_mdi <- fitdistrBayes(x_el, "el", "mdi", seed = 503)
el_reference_theta <- fitdistrBayes(x_el, "el", "reference-theta", seed = 504)
el_reference_rate <- fitdistrBayes(x_el, "el", "reference-rate", seed = 505)
list(
  Jeffreys = el_jeffreys$summary[, shown],
  MDI = el_mdi$summary[, shown],
  Reference_theta = el_reference_theta$summary[, shown],
  Reference_rate = el_reference_rate$summary[, shown]
)

# Rician ----------------------------------------------------------------------
set.seed(601)
x_rician <- sqrt(rnorm(80, 5, 2)^2 + rnorm(80, 0, 2)^2)
rician_jeffreys <- fitdistrBayes(x_rician, "rician", "jeffreys", seed = 602)
list(Jeffreys = rician_jeffreys$summary[, shown])

# Weighted Lindley ------------------------------------------------------------
set.seed(701)
lambda_wl <- 2.5
phi_wl <- 0.8
first_component <- runif(80) < lambda_wl / (lambda_wl + phi_wl)
x_wl <- rgamma(80, shape = phi_wl + as.numeric(!first_component),
               rate = lambda_wl)
wl_reference_lambda <- fitdistrBayes(
  x_wl, "weighted lindley", "reference-lambda", seed = 702
)
wl_reference_phi <- fitdistrBayes(
  x_wl, "weighted lindley", "reference-phi", seed = 703
)
list(
  Reference_lambda = wl_reference_lambda$summary[, shown],
  Reference_phi = wl_reference_phi$summary[, shown]
)

# Inspect the automatic classical starts used by the samplers.
gumbel_jeffreys$initialization
frechet_jeffreys$initialization
lomax_jeffreys$initialization
nakagami_jeffreys$initialization
el_jeffreys$initialization
rician_jeffreys$initialization
wl_reference_lambda$initialization

# Standard methods remain available for every fit.
summary(rician_jeffreys)
confint(rician_jeffreys)
head(as.data.frame(rician_jeffreys))
predict(rician_jeffreys, draws = 10, size = 5, seed = 603)

if (interactive()) {
  plot(rician_jeffreys, type = "trace")
  plot(rician_jeffreys, type = "density")
}
