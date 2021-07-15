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
#' @param deriv A character, indicating according to which hyper-parameter the
#'  derivative should be computed. If NULL (defaut), the function simply returns
#'  the covariance matrix.
#'
#' @export
#'
#' @return A covariance matrix, where elements are evaluations of the associated
#' kernel for each pair of reference inputs.
#'
#' @examples
#' kern_to_cov(
#'   rbind(c(1, 0, 1), c(2, 1, 2), c(1, 2, 3)),
#'   "SE",
#'   tibble::tibble(variance = 1, lengthscale = 0.5)
#' )
kern_to_cov <- function(input, kern = "SE", hp, deriv = NULL) {
  if (is.character(kern)) {
    if (kern == "SE") {
      kernel <- se_kernel
    }
    else if (kern == "PERIO") {
      kernel <- perio_kernel
    }
    else if (kern == "RQ") {
      kernel <- rq_kernel
    }
  }
  else if (is.function(kern)) {
    kernel <- kern
  }
  else {
    stop("Error in the 'kern' argument: please choose 'SE', 'PERIO', 'RQ', or
      provide a valid custom kernel function")
  }

  # Transform the batches of input into lists
  if (input %>% is.vector()) {
    list_input <- input
    reference <- as.character(input)
  }
  else {
    list_input <- split(t(input), rep(1:nrow(input), each = ncol(input)))
    if ("Input" %in% colnames(input)) {
      reference <- input %>%
        tibble::as_tibble() %>%
        dplyr::pull("Input") %>%
        as.character()
    } else {
      reference <- as.character(input[, 1])
    }
  }

  ## Detect whether the custom kernel provides derivatives
  if ("deriv" %in% methods::formalArgs(kernel)) {
    mat <- outer(
      list_input, list_input,
      Vectorize(function(x, y) kernel(x, y, hp, deriv = deriv))
    )
  }
  else {
    mat <- outer(
      list_input, list_input,
      Vectorize(function(x, y) kernel(x, y, hp))
    )
  }
  mat %>%
    `rownames<-`(reference) %>%
    `colnames<-`(reference) %>%
    return()
}

#' Create inverse of a covariance matrix from a kernel
#'
#' \code{kern_to_inv()} creates the inverse of a covariance matrix between
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
#' @param deriv A character, indicating according to which hyper-parameter the
#'  derivative should be computed. If NULL (defaut), the function simply returns
#'  the inverse covariance matrix.
#'
#' @return The inverse of a covariance matrix, which elements are evaluations of
#'    the associated kernel for each pair of reference inputs.
#'
#' @export
#'
#' @examples
#' kern_to_inv(
#'   rbind(c(1, 0, 1), c(2, 1, 2), c(1, 2, 3)),
#'   "SE",
#'   tibble::tibble(variance = 1, lengthscale = 0.5)
#' )
kern_to_inv <- function(input, kern, hp, pen_diag = 0, deriv = NULL) {
  mat_cov <- kern_to_cov(input = input, kern = kern, hp = hp, deriv = deriv)
  reference <- row.names(mat_cov)
  diag <- diag(x = pen_diag, ncol = ncol(mat_cov), nrow = nrow(mat_cov))



  inv <- tryCatch((mat_cov + diag) %>% chol() %>% chol2inv(),
    error = function(e) {
      MASS::ginv(mat_cov + diag)
    }
  )
  inv %>%
    `rownames<-`(reference) %>%
    `colnames<-`(reference) %>%
    return()
}

#' Compute a covariance matrix for multiple individuals
#'
#' Compute the covariance matrices associated with all individuals in the
#' database, taking into account their specific inputs and hyper-parameters.
#'
#' @param data A tibble or data frame of input data. Required column: 'ID'.
#'   Suggested column: 'Input' (for indicating the reference input).
#' @param kern A kernel function.
#' @param hp A tibble or data frame, containing the hyper-parameters associated
#' with each individual.
#' @param deriv A character, indicating according to which hyper-parameter the
#'  derivative should be computed. If NULL (defaut), the function simply returns
#'  the list of covariance matrices.
#'
#' @return A named list containing all of the inverse covariance matrices.
#'
#' @examples
#' db <- simu_db(M = 3)
#' hp <- tibble::tibble(ID = unique(db$ID), hp())
#' list_kern_to_cov(db, "SE", hp)
list_kern_to_cov <- function(data, kern, hp, deriv = NULL) {
  floop <- function(i) {
    db_i <- data %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(-.data$ID)
    ## To avoid throwing an error if 'Output' has already been removed
    if ("Output" %in% names(db_i)) {
      db_i <- db_i %>% dplyr::select(-.data$Output)
    }

    hp_i <- hp %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(-.data$ID)

    kern_to_cov(db_i, "SE", hp_i, deriv = deriv) %>%
      return()
  }
  sapply(unique(data$ID), floop, simplify = F, USE.NAMES = T) %>%
    return()
}

#' Compute an inverse covariance matrix for multiple individuals
#'
#' Compute the inverse covariance matrices associated with all individuals
#' in the database, taking into account their specific inputs and
#' hyper-parameters.
#'
#' @param db A tibble or data frame of input data. Required column: 'ID'.
#'   Suggested column: 'Input' (for indicating the reference input).
#' @param kern A kernel function.
#' @param hp A tibble or data frame, containing the hyper-parameters associated
#' with each individual.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#' numerical issues when inverting nearly singular matrices.
#' @param deriv A character, indicating according to which hyper-parameter the
#'  derivative should be computed. If NULL (defaut), the function simply returns
#'  the list of covariance matrices.
#'
#' @return A named list containing all of the inverse covariance matrices.
#'
#' @examples
#' db <- simu_db(M = 3)
#' hp <- tibble::tibble(ID = unique(db$ID), hp())
#' list_kern_to_inv(db, "SE", hp, 0)
list_kern_to_inv <- function(db, kern, hp, pen_diag = 0, deriv = NULL) {
  floop <- function(i) {
    db_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(-.data$ID)
    ## To avoid throwing an error if 'Output' has already been removed
    if ("Output" %in% names(db_i)) {
      db_i <- db_i %>% dplyr::select(-.data$Output)
    }

    hp_i <- hp %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(- .data$ID)

    kern_to_inv(db_i, "SE", hp_i, pen_diag, deriv = deriv) %>%
      return()
  }
  sapply(unique(db$ID), floop, simplify = F, USE.NAMES = T) %>%
    return()
}
