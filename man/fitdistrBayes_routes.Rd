\name{fitdistrBayes_routes}
\alias{fitdistrBayes_routes}
\title{List the Built-In Objective Bayesian Fitting Routes}
\usage{
fitdistrBayes_routes(model = NULL)
}
\arguments{
  \item{model}{Optional recognized distribution name. If supplied, only the
  routes for that distribution are returned.}
}
\value{
A data frame with one row for every enabled model--prior combination and
the following character columns:
\describe{
  \item{\code{model}}{Canonical distribution name.}
  \item{\code{prior}}{Accepted objective-prior label.}
  \item{\code{parameters}}{Comma-separated parameters estimated by the route.}
  \item{\code{fixed}}{Parameters that the user must supply in \code{fixed}, or
  an empty string when none are required.}
  \item{\code{engine}}{Posterior simulation or numerical-integration strategy
  used by the route.}
  \item{\code{condition}}{Concise sample condition under which the posterior
  is proper. The fitting function performs the authoritative check on the
  supplied observations.}
}
}
\description{
Returns the machine-readable catalogue used to document and test the available
model--prior combinations. The condition column is a concise summary; the
fitting function remains the authoritative executable check for the observed
sample.
}
\examples{
fitdistrBayes_routes()
fitdistrBayes_routes("gamma")
fitdistrBayes_routes("t")
}
