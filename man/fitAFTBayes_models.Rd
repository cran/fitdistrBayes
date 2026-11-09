\name{fitAFTBayes_models}
\alias{fitAFTBayes_models}
\title{Catalogue of Objective Bayesian AFT Models}
\description{Lists the lifetime families, registered prior names, estimated
parameters and fixed model arguments supported by \code{fitAFTBayes()}.}
\usage{fitAFTBayes_models()}
\value{A data frame with one row per family and character columns
\code{distribution}, \code{priors}, \code{estimated}, and
\code{fixed_arguments}. These describe accepted model names, available
objective-prior routes, posterior parameter columns, and required/optional
fixed parameters. The catalogue does not certify propriety for arbitrary
data; the fitting function performs the model-specific checks.}
\details{There are nine families. Student degrees of freedom are fixed.
Gamma shape is estimated by default and may optionally be fixed. Exponential
and Rayleigh have fixed scale, so their two registered prior names represent
the same flat measure on regression coefficients.}
\seealso{\code{\link{fitAFTBayes}}}
\examples{fitAFTBayes_models()}
