% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gradients-likelihoods.R
\name{gr_GP_mod_common_hp}
\alias{gr_GP_mod_common_hp}
\title{Gradient of the modified logLikelihood with common HPs for GPs in Magma}
\usage{
gr_GP_mod_common_hp(hp, db, mean, kern, post_cov, pen_diag)
}
\arguments{
\item{hp}{A tibble or data frame containing hyper-parameters for all
individuals.}

\item{db}{A tibble containing the values we want to compute the logL on.
Required columns: ID, Input, Output. Additional covariate columns are
allowed.}

\item{mean}{A vector, specifying the mean of the GPs at the reference inputs.}

\item{kern}{A kernel function.}

\item{post_cov}{A matrix, covariance parameter of the hyper-posterior.
Used to compute the correction term.}

\item{pen_diag}{A jitter term that is added to the covariance matrix to avoid
numerical issues when inverting, in cases of nearly singular matrices.}
}
\value{
A named vector, corresponding to the value of the hyper-parameters'
gradients for the modified Gaussian log-Likelihood involved in Magma with
the 'common HP' setting.
}
\description{
Gradient of the modified logLikelihood with common HPs for GPs in Magma
}
\examples{
TRUE
}
\keyword{internal}
