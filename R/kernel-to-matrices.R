#' Create covariance matrix from a kernel
#'
#' \code{kern_to_cov()} creates a covariance matrix between input values (that
#'    could be either scalars or vectors) evaluated within a kernel function,
#'    which is characterised by specified hyper-parameters. This matrix is
#'    a finite-dimensional evaluation of the infinite-dimensional covariance
#'    structure of a GP, defined thanks to this kernel.
#'
#' @param input A vector, matrix, data frame or tibble containing all inputs for
#'    one individual. If a vector, the elements are used as reference, otherwise,
#'    one column should be named 'Input' to indicate that it represents the
#'    reference (e.g. 'Input' would contain the timestamps in time-series
#'    applications). The other columns are considered as being covariates. If
#'    no column is named 'Input', the first one is used by default.
#' @param kern A kernel function. Several popular kernels
#'    (see \href{https://www.cs.toronto.edu/~duvenaud/cookbook/}{The Kernel
#'    Cookbook}) are already implemented and can be selected within the
#'    following list:
#'    - "SE": (default value) the Squared Exponential Kernel (also called
#'        Radial Basis Function or Gaussian kernel),
#'    - "PERIO": the Periodic kernel,
#'    - "RQ": the Rational Quadratic kernel.
#' @param hp A list, data frame or tibble containing the hyper-parameters used
#'   in the kernel. The name of the elements (or columns) should correspond
#'   exactly to those used in the kernel definition.
#'
#' @export
#'
#' @return A covariance matrix, where elements are evaluations of the associated
#' kernel for each pair of reference inputs.
#'
#' @examples
#' \dontrun{kern_to_cov(
#'    rbind(c(1, 0, 1), c(2, 1, 2), c(1, 2, 3)),
#'    "SE",
#'    tibble::tibble(sigma = 1, lengthscale = 0.5)
#'    )}
kern_to_cov <- function(input, kern = "SE", hp) {

  kern = dplyr::case_when(
    kern == "SE" ~ se_kernel,
    kern == "PERIO" ~ perio_kernel,
    kern == "RQ" ~ rq_kernel,
    TRUE ~ kern
  )

  # Transform the batches of input into lists
  list_input <- split(t(input), rep(1:ncol(input), each = nrow(input)))

  outer(list_input, list_input, Vectorize(function(x, y) kern(x, y, hp))) %>%
    return()
}


#' Create inverse of a covariance matrix from a kernel
#'
#'\code{kern_to_inv()} creates the inverse of a covariance matrix between
#'    input values (that could be either scalars or vectors) evaluated within
#'    a kernel function, which is characterised by specified hyper-parameters.
#'    This matrix is a finite-dimensional evaluation of the
#'    infinite-dimensional covariance structure of a GP, defined thanks to this
#'    kernel.
#'
#' @param input A vector, matrix, data frame or tibble containing all inputs for
#'    one individual. If a vector, the elements are used as reference, otherwise,
#'    one column should be named 'Input' to indicate that it represents the
#'    reference (e.g. 'Input' would contain the timestamps in time-series
#'    applications). The other columns are considered as being covariates. If
#'    no column is named 'Input', the first one is used by default.
#' @param kern A kernel function. Several popular kernels
#'    (see \href{https://www.cs.toronto.edu/~duvenaud/cookbook/}{The Kernel
#'    Cookbook}) are already implemented and can be selected within the
#'    following list:
#'    - "SE": (default value) the Squared Exponential Kernel (also called
#'        Radial Basis Function or Gaussian kernel),
#'    - "PERIO": the Periodic kernel,
#'    - "RQ": the Rational Quadratic kernel.
#' @param hp A list, data frame or tibble containing the hyper-parameters used
#'   in the kernel. The name of the elements (or columns) should correspond
#'   exactly to those used in the kernel definition.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return The inverse of a covariance matrix, which elements are evaluations of
#'    the associated kernel for each pair of reference inputs.
#'
#' @export
#'
#' @examples
#' \dontrun{kern_to_inv(
#'    rbind(c(1, 0, 1), c(2, 1, 2), c(1, 2, 3)),
#'    kernel_sqrd_exp,
#'    tibble::tibble(sigma = 1, lengthscale = 0.5)
#'    )}
kern_to_inv <- function(input, kern, hp, pen_diag = 0) {

  mat_cov <- kern_to_cov(input = input, kern = kern, hp = hp)
  diag <- diag(x = pen_diag, ncol = ncol(mat_cov), nrow = nrow(mat_cov))

  (mat_cov + diag) %>%
    solve() %>%
    return()
}
