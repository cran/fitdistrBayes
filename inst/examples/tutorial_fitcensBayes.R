# fitcensBayes: run one block at a time in the console.
# Install fitdistrBayes 0.5.0 first. This tutorial does not install anything
# and does not save outputs. Long MCMC examples are opt-in, not help examples.
# status = 1: exact event; status = 0: strictly T > x, also for counts.
# The first example uses random independent censoring; the others use fixed
# upper limits. On the real line these can represent detection limits.
# Priors within each model are grouped for comparison. The priors are those
# of the complete-data model, not derived for a particular censoring design.
# Missing mean/SD entries indicate infinite or uncertified posterior moments.
# coef() returns medians; confint() returns equal-tail credible intervals.
# If diagnostics fail, increase iterations and inspect trace/ACF plots.
# Twenty distributions and all 56 registered model-prior routes are shown;
# two additional calls illustrate data augmentation.

library(fitdistrBayes)

set.seed(101)
lifetime <- rexp(150, rate = 0.7)
censoring <- rexp(150, rate = 0.3)
x <- pmin(lifetime, censoring)
status <- as.integer(lifetime <= censoring)
table(status)

exp_J <- fitcensBayes(x, status, "exponential", "jeffreys", seed = 1)
exp_R <- fitcensBayes(x, status, "exponential", "reference", seed = 1)
exp_MDI <- fitcensBayes(x, status, "exponential", "mdi", seed = 1)
print(list(Jeffreys = exp_J, Reference = exp_R, MDI = exp_MDI))

summary(exp_R)
coef(exp_R)
confint(exp_R)
head(exp_R$draws)
exp_R$diagnostics
head(log_lik_cens(exp_R, draws = 5))

S <- predict(exp_R, type = "survival", times = c(0, 1, 2, 5), draws = 1000, seed = 2)
apply(S, 2, quantile, probs = c(0.025, 0.5, 0.975))
imputed_lifetimes <- predict(exp_R, type = "impute", draws = 5, seed = 3)
imputed_lifetimes[, 1:min(5, ncol(imputed_lifetimes)), drop = FALSE]

set.seed(102)
lifetime <- rgamma(150, shape = 2.5, rate = 1.3)
x <- pmin(lifetime, 2.5)
status <- as.integer(lifetime <= 2.5)
gamma_J <- fitcensBayes(x, status, "gamma", "jeffreys", seed = 2)
gamma_FR <- fitcensBayes(x, status, "gamma", "first-rule", seed = 2)
gamma_RS <- fitcensBayes(x, status, "gamma", "reference-shape", seed = 2)
gamma_RR <- fitcensBayes(x, status, "gamma", "reference-rate", seed = 2)
print(list(Jeffreys = gamma_J, First_rule = gamma_FR,
           Reference_shape = gamma_RS, Reference_rate = gamma_RR))

gamma_DA <- fitcensBayes(x, status, "gamma", "jeffreys",
                        method = "augmentation", seed = 3)
print(list(Direct = gamma_J, Augmentation = gamma_DA))

set.seed(103)
lifetime <- rweibull(150, shape = 1.7, scale = 2)
x <- pmin(lifetime, 2.3)
status <- as.integer(lifetime <= 2.3)
wei_J <- fitcensBayes(x, status, "weibull", "jeffreys", seed = 3)
wei_R <- fitcensBayes(x, status, "weibull", "reference", seed = 3)
print(list(Jeffreys = wei_J, Reference = wei_R))

set.seed(104)
lifetime <- rlnorm(150, meanlog = 0.3, sdlog = 0.7)
x <- pmin(lifetime, 2)
status <- as.integer(lifetime <= 2)
ln_J <- fitcensBayes(x, status, "lognormal", "jeffreys", seed = 4)
ln_R <- fitcensBayes(x, status, "lognormal", "reference", seed = 4)
print(list(Jeffreys = ln_J, Reference = ln_R))
ln_DA <- fitcensBayes(x, status, "lognormal", "reference",
                     method = "augmentation", seed = 4)
print(ln_DA)

set.seed(105)
lifetime <- rchisq(150, df = 4)
x <- pmin(lifetime, 5)
status <- as.integer(lifetime <= 5)
chi_J <- fitcensBayes(x, status, "chi-squared", "jeffreys", seed = 5)
chi_R <- fitcensBayes(x, status, "chi-squared", "reference", seed = 5)
print(list(Jeffreys = chi_J, Reference = chi_R))

set.seed(106)
lifetime <- (2 / rexp(150))^(1 / 3)
x <- pmin(lifetime, 1.6)
status <- as.integer(lifetime <= 1.6)
fre_J <- fitcensBayes(x, status, "frechet", "jeffreys", seed = 6)
fre_R <- fitcensBayes(x, status, "frechet", "reference", seed = 6)
print(list(Jeffreys = fre_J, Reference = fre_R))

set.seed(107)
lifetime <- 2 * expm1(rexp(150) / 3)
x <- pmin(lifetime, 1)
status <- as.integer(lifetime <= 1)
lom_J <- fitcensBayes(x, status, "lomax", "jeffreys", seed = 7)
print(lom_J)

set.seed(108)
lifetime <- sqrt(rgamma(150, shape = 2, rate = 2 / 3))
x <- pmin(lifetime, 2)
status <- as.integer(lifetime <= 2)
nak_J <- fitcensBayes(x, status, "nakagami", "jeffreys", seed = 8)
nak_R <- fitcensBayes(x, status, "nakagami", "reference", seed = 8)
print(list(Jeffreys = nak_J, Reference = nak_R))

set.seed(109)
theta <- 0.5
rate <- 1
lifetime <- -log((-expm1((1 - runif(150)) * log(theta))) / (1 - theta)) / rate
x <- pmin(lifetime, 1.1)
status <- as.integer(lifetime <= 1.1)
el_J <- fitcensBayes(x, status, "el", "jeffreys", seed = 9)
el_RT <- fitcensBayes(x, status, "el", "reference-theta", seed = 9)
el_RR <- fitcensBayes(x, status, "el", "reference-rate", seed = 9)
el_MDI <- fitcensBayes(x, status, "el", "mdi", seed = 9)
print(list(Jeffreys = el_J, Reference_theta = el_RT, Reference_rate = el_RR, MDI = el_MDI))

set.seed(110)
lifetime <- sqrt(rnorm(150, 2, 0.7)^2 + rnorm(150, 0, 0.7)^2)
x <- pmin(lifetime, 2.6)
status <- as.integer(lifetime <= 2.6)
ric_J <- fitcensBayes(x, status, "rician", "jeffreys", seed = 10)
print(ric_J)

set.seed(111)
lambda <- 1.5
phi <- 2
component <- rbinom(150, size = 1, prob = phi / (lambda + phi))
lifetime <- rgamma(150, shape = phi + component, rate = lambda)
x <- pmin(lifetime, 2)
status <- as.integer(lifetime <= 2)
wl_J <- fitcensBayes(x, status, "weighted lindley", "jeffreys", seed = 11)
wl_R <- fitcensBayes(x, status, "weighted lindley", "reference", seed = 11)
wl_FR <- fitcensBayes(x, status, "weighted lindley", "first-rule", seed = 11)
wl_IJ <- fitcensBayes(x, status, "weighted lindley", "independence-jeffreys", seed = 11)
wl_RL <- fitcensBayes(x, status, "weighted lindley", "reference-lambda", seed = 11)
wl_RP <- fitcensBayes(x, status, "weighted lindley", "reference-phi", seed = 11)
print(list(Jeffreys = wl_J, Reference = wl_R, First_rule = wl_FR,
           Independence_Jeffreys = wl_IJ, Reference_lambda = wl_RL, Reference_phi = wl_RP))

set.seed(112)
measurement <- rbeta(150, shape1 = 2, shape2 = 5)
x <- pmin(measurement, 0.4)
status <- as.integer(measurement <= 0.4)
beta_J <- fitcensBayes(x, status, "beta", "jeffreys", seed = 12)
beta_R <- fitcensBayes(x, status, "beta", "reference", seed = 12)
print(list(Jeffreys = beta_J, Reference = beta_R))

set.seed(113)
measurement <- rnorm(150, mean = 1, sd = 0.7)
x <- pmin(measurement, 1.4)
status <- as.integer(measurement <= 1.4)
norm_J <- fitcensBayes(x, status, "normal", "jeffreys", seed = 13)
norm_R <- fitcensBayes(x, status, "normal", "reference", seed = 13)
norm_MDI <- fitcensBayes(x, status, "normal", "mdi", seed = 13)
print(list(Jeffreys = norm_J, Reference = norm_R, MDI = norm_MDI))

set.seed(114)
measurement <- rlogis(150, location = 1, scale = 0.7)
x <- pmin(measurement, 1.6)
status <- as.integer(measurement <= 1.6)
logis_J <- fitcensBayes(x, status, "logistic", "jeffreys", seed = 14)
logis_R <- fitcensBayes(x, status, "logistic", "reference", seed = 14)
logis_MDI <- fitcensBayes(x, status, "logistic", "mdi", seed = 14)
print(list(Jeffreys = logis_J, Reference = logis_R, MDI = logis_MDI))

set.seed(115)
measurement <- 1 - 0.7 * log(-log(runif(150)))
x <- pmin(measurement, 2)
status <- as.integer(measurement <= 2)
gum_J <- fitcensBayes(x, status, "gumbel", "jeffreys", seed = 15)
gum_R <- fitcensBayes(x, status, "gumbel", "reference", seed = 15)
gum_MDI <- fitcensBayes(x, status, "gumbel", "mdi", seed = 15)
print(list(Jeffreys = gum_J, Reference = gum_R, MDI = gum_MDI))

set.seed(116)
measurement <- rcauchy(150, location = 1, scale = 0.7)
x <- pmin(measurement, 1.6)
status <- as.integer(measurement <= 1.6)
cau_J <- fitcensBayes(x, status, "cauchy", "jeffreys", seed = 16)
cau_R <- fitcensBayes(x, status, "cauchy", "reference", seed = 16)
cau_MDI <- fitcensBayes(x, status, "cauchy", "mdi", seed = 16)
print(list(Jeffreys = cau_J, Reference = cau_R, MDI = cau_MDI))

set.seed(117)
measurement <- 1 + 0.7 * rt(150, df = 5)
x <- pmin(measurement, 1.6)
status <- as.integer(measurement <= 1.6)
t_J <- fitcensBayes(x, status, "t", "jeffreys", fixed = list(df = 5), seed = 17)
t_R <- fitcensBayes(x, status, "t", "reference", fixed = list(df = 5), seed = 17)
t_MDI <- fitcensBayes(x, status, "t", "mdi", fixed = list(df = 5), seed = 17)
t_IJ <- fitcensBayes(x, status, "t", "independence-jeffreys", seed = 17)
print(list(Jeffreys_fixed_df = t_J, Reference_fixed_df = t_R,
           MDI_fixed_df = t_MDI, Independence_Jeffreys_unknown_df = t_IJ))

set.seed(118)
measurement <- rgeom(150, prob = 0.3)
x <- pmin(measurement, 3)
status <- as.integer(measurement <= 3)
geo_J <- fitcensBayes(x, status, "geometric", "jeffreys", seed = 18)
geo_R <- fitcensBayes(x, status, "geometric", "reference", seed = 18)
geo_MDI <- fitcensBayes(x, status, "geometric", "mdi", seed = 18)
print(list(Jeffreys = geo_J, Reference = geo_R, MDI = geo_MDI))

set.seed(119)
measurement <- rpois(150, lambda = 5)
x <- pmin(measurement, 6)
status <- as.integer(measurement <= 6)
poi_J <- fitcensBayes(x, status, "Poisson", "jeffreys", seed = 19)
poi_R <- fitcensBayes(x, status, "Poisson", "reference", seed = 19)
poi_MDI <- fitcensBayes(x, status, "Poisson", "mdi", seed = 19)
print(list(Jeffreys = poi_J, Reference = poi_R, MDI = poi_MDI))

set.seed(120)
measurement <- rnbinom(150, size = 2, mu = 5)
x <- pmin(measurement, 6)
status <- as.integer(measurement <= 6)
nb_J <- fitcensBayes(x, status, "negative binomial", "jeffreys", fixed = list(size = 2), seed = 20)
nb_R <- fitcensBayes(x, status, "negative binomial", "reference", fixed = list(size = 2), seed = 20)
nb_MDI <- fitcensBayes(x, status, "negative binomial", "mdi", fixed = list(size = 2), seed = 20)
print(list(Jeffreys = nb_J, Reference = nb_R, MDI = nb_MDI))

fitcensBayes_models()
