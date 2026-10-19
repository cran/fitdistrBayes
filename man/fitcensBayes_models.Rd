\name{fitcensBayes_models}
\alias{fitcensBayes_models}
\title{Catalogue of Registered Censored Model-Prior Routes}
\description{
Lists distributions and complete-data objective priors supported by the
separate right-censoring interface, including sufficient propriety
conditions and sampling options.
}
\usage{fitcensBayes_models()}
\value{
A data frame with 59 rows (one per registered model-prior route) for
21 distributions. Its character columns are:
\describe{
  \item{\code{model}}{Canonical distribution name.}
  \item{\code{prior}}{Registered prior name accepted by the fitting function.}
  \item{\code{parameters}}{Names of parameters that are estimated.}
  \item{\code{required_fixed}}{Required fixed parameters, if any.}
  \item{\code{complete_data_engine}}{Sampler registered for complete data;
  not a claim that the same algorithm is used under censoring.}
  \item{\code{exact_subset_sufficient_condition}}{Sufficient condition
  applied to the exact events alone. Failure means not certified, not
  necessarily improper. Dhillon instead lists its exact observed-data
  condition. See also the analytic Exponential exception below.}
  \item{\code{censored_methods}}{Available observed-likelihood and
  data-augmentation methods.}
  \item{\code{censoring}}{Right-censoring status convention.}
  \item{\code{prior_origin}}{The complete-data sampling model.}
}
}
\details{
Prior formulas, parameterizations, and aliases are inherited from
\code{\link{fitdistrBayes}} and \code{\link{fitdistrBayes_routes}}.
The catalogue is not an assertion that objective priors are invariant to
the censoring design. Conditions are sufficient and conservative except
where an exact model-specific condition is explicitly documented.
Exponential MDI additionally permits zero exact events when the total
observed time is positive and finite; the resulting Gamma posterior
has shape two. The function documentation describes this analytic case.
Dhillon requires an event and either distinct exact times or a censoring
limit above the common event time; see \link{Dhillon}.
}
\seealso{\code{\link{fitcensBayes}}, \code{\link{fitdistrBayes_routes}}}
\examples{
routes <- fitcensBayes_models()
routes[routes$model == "gamma", c("model", "prior", "parameters")]
length(unique(routes$model))
nrow(routes)
}
\keyword{models}
