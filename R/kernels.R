#' Squared Exponential Kernel
#'
#' @param x A vector of inputs.
#' @param y A vector of inputs.
#' @param hp A tibble, data frame or named vector, containing the kernel's
#'    hyperparameters. Required columns: 'variance', 'lengthscale', and 'scale'.
#' @param deriv A character, indicating according to which hyper-parameter the
#'  derivative should be computed. If NULL (defaut), the function simply returns
#'  the evaluation of the kernel.
#'
#' @return A scalar, conresponding to the evaluation of the kernel.
#'
#' @examples
#' se_kernel(
#'   c(1, 0), c(0, 1),
#'   tibble::tibble(variance = 1, lengthscale = 0.5)
#' )
se_kernel <- function(x, y, hp, deriv = NULL) {
  top_term <- exp(-hp[['lengthscale']]) * 0.5 * sum((x - y)^2)

  if(deriv %>% is.null()){
    (exp(hp[['variance']] - top_term)) %>% return()
  }
  else if(deriv == 'variance'){
    (exp(hp[['variance']] - top_term)) %>% return()
  }
  else if(deriv == 'lengthscale'){
    (exp(hp[['variance']]) * top_term * exp(-top_term)) %>% return()
  }
  else{
    stop("Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'")
  }
}

#' Periodic Kernel
#'
#' @param x A vector of inputs.
#' @param y A vector of inputs.
#' @param hp A tibble, data frame or named vector, containing the kernel's
#'    hyperparameters. Required columns: 'variance', 'lengthscale', and 'scale'.
#' @param deriv A character, indicating according to which hyper-parameter the
#'  derivative should be computed. If NULL (defaut), the function simply returns
#'  the evaluation of the kernel.
#'
#' @return A scalar, conresponding to the evaluation of the kernel.
#'
#' @examples
#' perio_kernel(
#'   c(1, 0), c(0, 1),
#'   tibble::tibble(variance = 1, lengthscale = 0.5, period = 2)
#' )
perio_kernel <- function(x, y, hp, deriv = NULL) {
  distance <- abs(x - y) %>% sum()
  angle <- pi * distance / hp[['period']]

  if(deriv %>% is.null()){
    (exp(hp[['variance']])*exp(-2 * sin(angle)^2*exp(-hp[['lengthscale']]))) %>%
      return()
  }
  else if(deriv == 'variance'){
    (exp(hp[['variance']])*exp(-2 * sin(angle)^2*exp(-hp[['lengthscale']]))) %>%
      return()
  }
  else if(deriv == 'period'){
    (exp(hp[['variance']]) * (pi * distance / hp[['period']]^2) *
      cos(angle) * (2 * exp(-hp[['lengthscale']]) * 2 * sin(angle)) *
      exp(-2 * sin(angle)^2 * exp(-hp[['lengthscale']]))) %>%
      return()
  }
  else if(deriv == 'lengthscale'){
    (exp(hp[['variance']]) * 2 * sin(angle)^2 * exp(-hp[['lengthscale']]) *
      exp(-2 * sin(angle)^2 * exp(-hp[['lengthscale']]))) %>%
      return()
  }
  else{
    stop("Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'")
  }
}

#' Rational Quadratic Kernel

#' @param x A vector of inputs.
#' @param y A vector of inputs.
#' @param hp A tibble, data frame or named vector, containing the kernel's
#'    hyperparameters. Required columns: 'variance', 'lengthscale', and 'scale'.
#' @param deriv A character, indicating according to which hyper-parameter the
#'  derivative should be computed. If NULL (defaut), the function simply returns
#'  the evaluation of the kernel.
#'
#' @return A scalar, conresponding to the evaluation of the kernel.
#'
#' @examples
#' rq_kernel(
#'   c(1, 0), c(0, 1),
#'   tibble::tibble(variance = 1, lengthscale = 0.5, scale = 3)
#' )
rq_kernel <- function(x, y, hp, deriv = NULL) {
  distance <- sum((x - y)^2)
  term <- (1 + distance * exp(-hp[['lengthscale']]) / (2 * hp[['scale']]))

  if(deriv %>% is.null()){
    (exp(hp[['variance']]) * term^(-hp[['scale']])) %>%
      return()
  }
  else if(deriv == 'variance'){
    (exp(hp[['variance']]) * term^(-hp[['scale']])) %>%
      return()
  }
  else if(deriv == 'scale'){
    (exp(hp[['variance']]) * term^(-hp[['scale']]) *
       (distance * exp(-hp[['lengthscale']]) / (2 * hp[['scale']] * term)
        - log(term))) %>%
      return()
  }
  else if(deriv == 'lengthscale'){
    (exp(hp[['variance']]) * distance * 0.5 *
       exp(-hp[['lengthscale']]) * term^(-hp[['scale']] - 1)) %>%
      return()
  }
  else{
    stop("Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'")
  }
}

#' Generate random hyper-parameters
#'
#' Generate a set of random hyper-parameters, specific to the chosen type of
#'  kernel.
#'
#' @param kern A character string indicating the type of kernel among:
#'  - "SE": the Squared Exponential kernel,
#'  - "PERIO": the Periodic kernel,
#'  - "RQ": the Rational Quadratic kernel.
#' @param list_ID A vector, associating an \code{ID} value with each individual
#'    for whom hyper-parameters are generated. If NULL (defaut) only one set of
#'    hyper-parameters is return without the \code{ID} column.
#'
#' @return A tibble, gathering a set of hyper-parameters.
#'
#' @examples
#' hp("PERIO")
hp <- function(kern = "SE", list_ID = NULL) {
  if (is.null(list_ID)) {
    if (kern == "SE") {
      hp <- tibble::tibble(
        variance = stats::runif(1, 1, 5),
        lengthscale = stats::runif(1, 1, 5)
      )
    }
    else if (kern == "PERIO") {
      hp <- tibble::tibble(
        variance = stats::runif(1, 1, 5),
        lengthscale = stats::runif(1, 1, 5),
        period = stats::runif(1, 0, 2 * pi)
      )
    }
    else if (kern == "RQ") {
      hp <- tibble::tibble(
        variance = stats::runif(1, 1, 5),
        lengthscale = stats::runif(1, 1, 5),
        scale = stats::runif(1, 1, 5)
      )
    }
  }
  else {
    len <- length(list_ID)
    if (kern == "SE") {
      hp <- tibble::tibble(
        ID = list_ID,
        variance = stats::runif(len, 1, 5),
        lengthscale = stats::runif(len, 1, 5)
      )
    }
    else if (kern == "PERIO") {
      hp <- tibble::tibble(
        ID = list_ID,
        variance = stats::runif(len, 1, 5),
        lengthscale = stats::runif(len, 1, 5),
        period = stats::runif(len, 0, 2 * pi)
      )
    }
    else if (kern == "RQ") {
      hp <- tibble::tibble(
        ID = list_ID,
        variance = stats::runif(len, 1, 5),
        lengthscale = stats::runif(len, 1, 5),
        scale = stats::runif(len, 1, 5)
      )
    }
  }
  return(hp)
}
