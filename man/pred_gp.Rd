% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction.R
\name{pred_gp}
\alias{pred_gp}
\title{Gaussian Process prediction}
\usage{
pred_gp(
  data = NULL,
  grid_inputs = NULL,
  mean = NULL,
  hp = NULL,
  kern = "SE",
  get_full_cov = FALSE,
  plot = TRUE,
  pen_diag = 1e-10
)
}
\arguments{
\item{data}{A tibble or data frame. Required columns: 'Input',
'Output'. Additional columns for covariates can be specified.
The 'Input' column should define the variable that is used as
reference for the observations (e.g. time for longitudinal data). The
'Output' column specifies the observed values (the response
variable). The data frame can also provide as many covariates as desired,
with no constraints on the column names. These covariates are additional
inputs (explanatory variables) of the models that are also observed at
each reference 'Input'. If NULL, the prior GP is returned.}

\item{grid_inputs}{The grid of inputs (reference Input and covariates) values
on which the GP should be evaluated. Ideally, this argument should be a
tibble or a data frame, providing the same columns as \code{data}, except
'Output'. Nonetheless, in cases where \code{data} provides only one
'Input' column, the \code{grid_inputs} argument can be NULL (default) or a
vector. This vector would be used as reference input for prediction and if
NULL, a vector of length 500 is defined, ranging between the min and max
Input values of \code{data}.}

\item{mean}{Mean parameter of the GP. This argument can be specified under
various formats, such as:
\itemize{
\item NULL (default). The mean would be set to 0 everywhere.
\item A number. The mean would be a constant function.
\item A tibble or data frame. Required columns: Input, Output. The Input
values should include at least the same values as in the \code{data}
argument.
}}

\item{hp}{A named vector, tibble or data frame of hyper-parameters
associated with \code{kern}. The columns/elements should be named
according to the hyper-parameters that are used in \code{kern}. If NULL
(default), the function \code{\link{train_gp}} is called with random
initial values for learning maximum-likelihood estimators of the
hyper-parameters associated with \code{kern}.}

\item{kern}{A kernel function, defining the covariance structure of the GP.
Several popular kernels
(see \href{https://www.cs.toronto.edu/~duvenaud/cookbook/}{The Kernel
Cookbook}) are already implemented and can be selected within the
following list:
\itemize{
\item "SE": (default value) the Squared Exponential Kernel (also called
Radial Basis Function or Gaussian kernel),
\item "LIN": the Linear kernel,
\item "PERIO": the Periodic kernel,
\item "RQ": the Rational Quadratic kernel.
Compound kernels can be created as sums or products of the above kernels.
For combining kernels, simply provide a formula as a character string
where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
elements are treated sequentially from the left to the right, the product
operator '*' shall always be used before the '+' operators (e.g.
'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
}}

\item{get_full_cov}{A logical value, indicating whether the full posterior
covariance matrix should be returned.}

\item{plot}{A logical value, indicating whether a plot of the results is
automatically displayed.}

\item{pen_diag}{A number. A jitter term, added on the diagonal to prevent
numerical issues when inverting nearly singular matrices.}
}
\value{
A tibble, representing the GP predictions as two column 'Mean' and
'Var', evaluated on the \code{grid_inputs}. The column 'Input' and
additional covariates columns are associated to each predicted values.
If the \code{get_full_cov} argument is TRUE, the function returns a list,
in which the tibble described above is defined as 'pred' and the full
posterior covariance matrix is defined as 'cov'.
}
\description{
Compute the posterior distribution of a standard GP, using the formalism of
Magma. By providing observed data, the prior mean and covariance
matrix (by defining a kernel and its associated hyper-parameters), the mean
and covariance parameters of the posterior distribution are computed on the
grid of inputs that has been specified. This predictive distribution can be
evaluated on any arbitrary inputs since a GP is an infinite-dimensional
object.
}
\examples{
TRUE
}
