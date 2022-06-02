#' Squared Exponential Kernel
#'
#' @param x A vector (or matrix if vectorized = T) of inputs.
#' @param y A vector (or matrix if vectorized = T) of inputs.
#' @param hp A tibble, data frame or named vector, containing the kernel's
#'    hyperparameters. Required columns: 'se_variance', 'se_lengthscale'.
#' @param deriv A character, indicating according to which hyper-parameter the
#'    derivative should be computed. If NULL (default), the function simply
#'    returns the evaluation of the kernel.
#' @param vectorized A logical value, indicating whether the function provides
#'    a vectorized version for speeded-up calculations. If TRUE, the \code{x}
#'    and \code{y} arguments should be the vector or matrix containing all
#'    inputs for which the kernel is evaluated on all pairs of elements.
#'    If FALSE, the \code{x} and \code{y} arguments are simply two inputs.
#'
#' @return A scalar, corresponding to the evaluation of the kernel.
#'
#' @examples
#' TRUE
se_kernel <- function(
  x,
  y,
  hp,
  deriv = NULL,
  vectorized = FALSE) {
  ## Check whether the Rcpp function for speed-up vectorised computation
  if (vectorized) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    ## Compute directly the full covariance matrix
    distance <- cpp_dist(x, y)
  } else {
    ## Compute one element of the covariance matrix
    distance <- sum((x - y)^2)
  }

  top_term <- exp(-hp[["se_lengthscale"]]) * 0.5 * distance

  if (deriv %>% is.null()) {
    (exp(hp[["se_variance"]] - top_term)) %>% return()
  } else if (deriv == "se_variance") {
    (exp(hp[["se_variance"]] - top_term)) %>% return()
  } else if (deriv == "se_lengthscale") {
    (exp(hp[["se_variance"]]) * top_term * exp(-top_term)) %>% return()
  } else {
    stop("Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'")
  }
}

#' Periodic Kernel
#'
#' @param x A vector (or matrix if vectorized = T) of inputs.
#' @param y A vector (or matrix if vectorized = T) of inputs.
#' @param hp A tibble, data frame or named vector, containing the kernel's
#'    hyperparameters. Required columns: 'perio_variance', 'perio_lengthscale',
#'    and 'period'.
#' @param deriv A character, indicating according to which hyper-parameter the
#'  derivative should be computed. If NULL (default), the function simply returns
#'  the evaluation of the kernel.
#' @param vectorized A logical value, indicating whether the function provides
#'    a vectorized version for speeded-up calculations. If TRUE, the \code{x}
#'    and \code{y} arguments should be the vector or matrix containing all
#'    inputs for which the kernel is evaluated on all pairs of elements.
#'    If FALSE, the \code{x} and \code{y} arguments are simply two inputs.
#'
#'
#' @return A scalar, corresponding to the evaluation of the kernel.
#'
#' @examples
#' TRUE
perio_kernel <- function(
  x,
  y,
  hp,
  deriv = NULL,
  vectorized = FALSE) {
  ## Check whether the Rcpp function for speed-up vectorised computation
  if (vectorized) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    ## Compute directly the full covariance matrix
    perio_term <- cpp_perio(x, y, hp[["period"]])
    sum_deriv <- cpp_perio_deriv(x, y, hp[["period"]])
  } else {
    ## Compute one element of the covariance matrix
    angle <- pi * abs(x - y) / exp(hp[["period"]])
    perio_term <- sin(angle)^2 %>% sum()
    sum_deriv <- sum(2 * sin(angle) * cos(angle) * angle)
  }

  if (deriv %>% is.null()) {
    (exp(hp[["perio_variance"]]) *
      exp(-2 * exp(-hp[["perio_lengthscale"]]) * perio_term)) %>%
      return()
  } else if (deriv == "perio_variance") {
    (exp(hp[["perio_variance"]]) *
      exp(-2 * exp(-hp[["perio_lengthscale"]]) * perio_term)) %>%
      return()
  } else if (deriv == "period") {
    (exp(hp[["perio_variance"]]) * exp(-2 * exp(-hp[["perio_lengthscale"]]) *
      perio_term) * 2 * exp(-hp[["perio_lengthscale"]]) * sum_deriv) %>%
      return()
  } else if (deriv == "perio_lengthscale") {
    (exp(hp[["perio_variance"]]) * 2 * perio_term *
      exp(-hp[["perio_lengthscale"]]) *
      exp(-2 * perio_term * exp(-hp[["perio_lengthscale"]]))) %>%
      return()
  } else {
    stop("Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'")
  }
}

#' Rational Quadratic Kernel

#' @param x A vector (or matrix if vectorized = T) of inputs.
#' @param y A vector (or matrix if vectorized = T) of inputs.
#' @param hp A tibble, data frame or named vector, containing the kernel's
#'    hyperparameters. Required columns: 'rq_variance', 'rq_lengthscale', and
#'    'rq_scale'.
#' @param deriv A character, indicating according to which hyper-parameter the
#'  derivative should be computed. If NULL (default), the function simply returns
#'  the evaluation of the kernel.
#' @param vectorized A logical value, indicating whether the function provides
#'    a vectorized version for speeded-up calculations. If TRUE, the \code{x}
#'    and \code{y} arguments should be the vector or matrix containing all
#'    inputs for which the kernel is evaluated on all pairs of elements.
#'    If FALSE, the \code{x} and \code{y} arguments are simply two inputs.
#'
#' @return A scalar, corresponding to the evaluation of the kernel.
#'
#' @examples
#' TRUE
rq_kernel <- function(
  x,
  y,
  hp,
  deriv = NULL,
  vectorized = FALSE) {
  ## Check whether the Rcpp function for speed-up vectorised computation
  if (vectorized) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    ## Compute directly the full covariance matrix
    distance <- cpp_dist(x, y)
  } else {
    ## Compute one element of the covariance matrix
    distance <- sum((x - y)^2)
  }

  term <- (1 + distance * exp(-hp[["rq_lengthscale"]]) / (2 * hp[["rq_scale"]]))

  if (deriv %>% is.null()) {
    (exp(hp[["rq_variance"]]) * term^(-hp[["rq_scale"]])) %>%
      return()
  } else if (deriv == "rq_variance") {
    (exp(hp[["rq_variance"]]) * term^(-hp[["rq_scale"]])) %>%
      return()
  } else if (deriv == "rq_scale") {
    (exp(hp[["rq_variance"]]) * term^(-hp[["rq_scale"]]) *
      (distance * exp(-hp[["rq_lengthscale"]]) / (2 * hp[["rq_scale"]] * term)
        - log(term))) %>%
      return()
  } else if (deriv == "rq_lengthscale") {
    (exp(hp[["rq_variance"]]) * distance * 0.5 *
      exp(-hp[["rq_lengthscale"]]) * term^(-hp[["rq_scale"]] - 1)) %>%
      return()
  } else {
    stop("Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'")
  }
}

#' Linear Kernel
#'
#' @param x A vector (or matrix if vectorized = T) of inputs.
#' @param y A vector (or matrix if vectorized = T) of inputs.
#' @param hp A tibble, data frame or named vector, containing the kernel's
#'    hyperparameters. Required columns: 'lin_slope' and 'lin_offset'.
#' @param deriv A character, indicating according to which hyper-parameter the
#'    derivative should be computed. If NULL (default), the function simply
#'    returns the evaluation of the kernel.
#' @param vectorized A logical value, indicating whether the function provides
#'    a vectorized version for speeded-up calculations. If TRUE, the \code{x}
#'    and \code{y} arguments should be the vector or matrix containing all
#'    inputs for which the kernel is evaluated on all pairs of elements.
#'    If FALSE, the \code{x} and \code{y} arguments are simply two inputs.
#'
#' @return A scalar, corresponding to the evaluation of the kernel.
#'
#' @examples
#' TRUE
lin_kernel <- function(
  x,
  y,
  hp,
  deriv = NULL,
  vectorized = FALSE) {
  ## Check whether the Rcpp function for speed-up vectorised computation
  if (vectorized) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    ## Compute directly the full covariance matrix
    prod <- cpp_prod(x, y)
    ## Create a dummy matrix of one for the 'lin_offset' derivative matrix
    mat_one <- matrix(1, ncol = ncol(prod), nrow = nrow(prod))
  } else {
    ## Compute one element of the covariance matrix
    prod <- x %*% y
    mat_one <- 1
  }

  if (deriv %>% is.null()) {
    (exp(hp[["lin_slope"]]) * prod + exp(hp[["lin_offset"]])) %>%
      return()
  } else if (deriv == "lin_offset") {
    exp(hp[["lin_offset"]] * mat_one) %>% return()
  } else if (deriv == "lin_slope") {
    (exp(hp[["lin_slope"]]) * prod) %>% return()
  } else {
    stop("Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'")
  }
}

#' Generate random hyper-parameters
#'
#' Generate a set of random hyper-parameters, specific to the chosen type of
#' kernel, under the format that is used in Magma.
#'
#' @param kern A function, or a character string indicating the chosen type of
#'    kernel among:
#'    - "SE": the Squared Exponential kernel,
#'    - "LIN": the Linear kernel,
#'    - "PERIO": the Periodic kernel,
#'    - "RQ": the Rational Quadratic kernel.
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#'
#'    In case of a custom kernel function, the argument \code{list_hp} has to be
#'    provided as well, for designing a tibble with the correct names of
#'    hyper-parameters.
#' @param list_ID A vector, associating an \code{ID} value with each individual
#'    for whom hyper-parameters are generated. If NULL (default) only one set of
#'    hyper-parameters is return without the \code{ID} column.
#' @param list_hp A vector of characters, providing the name of each
#'    hyper-parameter, in case where \code{kern} is a custom kernel function.
#' @param noise A logical value, indicating whether a 'noise' hyper-parameter
#'    should be included.
#' @param common_hp A logical value, indicating whether the set of
#'    hyper-parameters is assumed to be common to all indiviuals.
#'
#' @return A tibble, providing a set of random hyper-parameters associated with
#'   the kernel specified through the argument \code{kern}.
#'
#' @export
#'
#' @examples
#' TRUE
hp <- function(
  kern = "SE",
  list_ID = NULL,
  list_hp = NULL,
  noise = FALSE,
  common_hp = FALSE) {
  ## Initiate interval boundaries
  min <- 0
  max <- 3
  min_noise <- -2
  max_noise <- 0

  len <- 1

  ## Split the 'kern' argument into successive kernels
  if (kern %>% is.character()) {
    str_kern <- strsplit(kern, " +")[[1]]
  } else {
    str_kern <- 1
  }

  if (is.null(list_ID)) {
    hp <- tibble::tibble(.rows = 1)
  } else {
    hp <- tibble::tibble(ID = as.character(list_ID))
  }

  for (i in str_kern)
  {
    if (is.null(list_ID)) {
      if (is.function(kern) | is.null(kern)) {
        temp_hp <- tibble::tibble(.rows = 1)
        for (j in list_hp) {
          temp_hp[j] <- stats::runif(1, min, max)
        }
        hp <- hp %>% dplyr::bind_cols(temp_hp)
      } else if (i == "SE") {
        hp <- hp %>% dplyr::bind_cols(
          tibble::tibble(
            se_variance = stats::runif(1, min, max),
            se_lengthscale = stats::runif(1, min, max)
          )
        )
      } else if (i == "PERIO") {
        hp <- hp %>% dplyr::bind_cols(
          tibble::tibble(
            perio_variance = stats::runif(1, min, max),
            perio_lengthscale = stats::runif(1, min, max),
            period = stats::runif(1, 0, 2 * pi)
          )
        )
      } else if (i == "RQ") {
        hp <- hp %>% dplyr::bind_cols(
          tibble::tibble(
            rq_variance = stats::runif(1, min, max),
            rq_lengthscale = stats::runif(1, min, max),
            rq_scale = stats::runif(1, min, max)
          )
        )
      } else if (i == "LIN") {
        hp <- hp %>% dplyr::bind_cols(
          tibble::tibble(
            lin_slope = stats::runif(1, min, max),
            lin_offset = stats::runif(1, min, max),
          )
        )
      }
    } else {
      if (common_hp) len <- 1 else len <- length(list_ID)
      if (is.function(kern) | is.null(kern)) {
        temp_hp <- tibble::tibble(ID = as.character(list_ID))
        for (j in list_hp) {
          temp_hp[j] <- stats::runif(len, min, max)
        }
        hp <- hp %>% dplyr::left_join(temp_hp, by = "ID")
      } else if (i == "SE") {
        hp <- hp %>% dplyr::left_join(
          tibble::tibble(
            ID = as.character(list_ID),
            se_variance = stats::runif(len, min, max),
            se_lengthscale = stats::runif(len, min, max)
          ),
          by = "ID"
        )
      } else if (i == "PERIO") {
        hp <- hp %>% dplyr::left_join(
          tibble::tibble(
            ID = as.character(list_ID),
            perio_variance = stats::runif(len, min, max),
            perio_lengthscale = stats::runif(len, min, max),
            period = stats::runif(len, 0, 2 * pi)
          ),
          by = "ID"
        )
      } else if (i == "RQ") {
        hp <- hp %>% dplyr::left_join(
          tibble::tibble(
            ID = as.character(list_ID),
            rq_variance = stats::runif(len, min, max),
            rq_lengthscale = stats::runif(len, min, max),
            rq_scale = stats::runif(len, min, max)
          ),
          by = "ID"
        )
      } else if (i == "LIN") {
        hp <- hp %>% dplyr::left_join(
          tibble::tibble(
            ID = as.character(list_ID),
            lin_slope = stats::runif(len, min, max),
            lin_offset = stats::runif(len, min, max)
          ),
          by = "ID"
        )
      }
    }
  }
  ## Add noise if necessary
  if (noise) {
    hp <- hp %>%
      dplyr::mutate("noise" = stats::runif(len, min_noise, max_noise))
  }
  return(hp)
}
