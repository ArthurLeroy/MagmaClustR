#' Squared Exponential Kernel
#'
#' @param x A vector of inputs.
#' @param y A vector of inputs.
#' @param hp A tibble containing the kernel's hyperparameters.
#' Required columns: 'variance', 'lengthscale', and 'scale'.
#'
#' @return The evaluation of the kernel.
#'
#' @examples
#' se_kernel(
#'   c(1, 0), c(0, 1),
#'   tibble::tibble(variance = 1, lengthscale = 0.5)
#' )
se_kernel <- function(x, y, hp) {
  top_term <- exp(-hp$lengthscale) * 0.5 * sum((x - y)^2)
  kern <- exp(hp$variance - top_term)

  attr(kern, "variance") <- exp(hp$variance - top_term)

  attr(kern, "lengthscale") <- exp(hp$variance) * top_term * exp(-top_term)

  return(kern)
}


#' Periodic Kernel
#'
#' @param x A vector of inputs.
#' @param y A vector of inputs.
#' @param hp A tibble containing the kernel's hyperparameters.
#' Required columns: 'variance', 'lengthscale', and 'scale'.
#'
#' @return The evaluation of the kernel.
#'
#' @examples
#' kernel_period(
#'   c(1, 0), c(0, 1),
#'   tibble::tibble(variance = 1, lengthscale = 0.5, period = 2)
#' )
perio_kernel <- function(x, y, hp) {
  distance <- abs(x - y) %>% sum()
  angle <- pi * distance / hp$period

  kern <- exp(hp$variance) * exp(-2 * sin(angle)^2 * exp(-hp$lengthscale))

  attr(kern, "variance") <- exp(hp$variance) *
    exp(-2 * sin(angle)^2 * exp(-hp$lengthscale))

  attr(kern, "period") <- exp(hp$variance) * (pi * distance / hp$period^2) *
    cos(angle) * (2 * exp(-hp$lengthscale) * 2 * sin(angle)) *
    exp(-2 * sin(angle)^2 * exp(-hp$lengthscale))

  attr(kern, "lengthscale") <- exp(hp$variance) * 2 * sin(angle)^2 *
    exp(-hp$lengthscale) * exp(-2 * sin(angle)^2 * exp(-hp$lengthscale))

  return(kern)
}

#' Rational Quadratic Kernel

#' @param x A vector of inputs.
#' @param y A vector of inputs.
#' @param hp A tibble containing the kernel's hyperparameters.
#' Required columns: 'variance', 'lengthscale', and 'scale'.
#'
#' @return The evaluation of the kernel.
#'
#' @examples
#' kernel_quad(
#'   c(1, 0), c(0, 1),
#'   tibble::tibble(variance = 1, lengthscale = 0.5, scale = 3)
#' )
rq_kernel <- function(x, y, hp) {
  distance <- sum((x - y)^2)
  term <- (1 + distance * exp(-hp$lengthscale) / (2 * hp$scale))

  kern <- exp(hp$variance) * term^(-hp$scale)

  attr(kern, "variance") <- exp(hp$variance) * term^(-hp$scale)

  attr(kern, "scale") <- exp(hp$variance) * term^(-hp$scale) *
    (distance * exp(-hp$lengthscale) / (2 * hp$scale * term) - log(term))

  attr(kern, "lengthscale") <- exp(hp$variance) * distance * 0.5 *
    exp(-hp$lengthscale) * term^(-hp$scale - 1)

  return(kern)
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
#'
#' @return A tibble, gathering a set of hyper-parameters.
#'
#' @examples
#' hp("PERIO")
hp <- function(kern = "SE", list_ID = NULL) {
  if (is.null(list_ID)) {
    if (kern == "SE") {
      hp <- tibble::tibble(
        variance = runif(1, 1, 5),
        lengthscale = runif(1, 1, 5)
      )
    }
    else if (kern == "PERIO") {
      hp <- tibble::tibble(
        variance = runif(1, 1, 5),
        lengthscale = runif(1, 1, 5),
        period = runif(1, 0, 2 * pi)
      )
    }
    else if (kern == "RQ") {
      hp <- tibble::tibble(
        variance = runif(1, 1, 5),
        lengthscale = runif(1, 1, 5),
        scale = runif(1, 1, 5)
      )
    }
  }
  else {
    len <- length(list_ID)
    if (kern == "SE") {
      hp <- tibble::tibble(
        ID = list_ID,
        variance = runif(len, 1, 5),
        lengthscale = runif(len, 1, 5)
      )
    }
    else if (kern == "PERIO") {
      hp <- tibble::tibble(
        ID = list_ID,
        variance = runif(len, 1, 5),
        lengthscale = runif(len, 1, 5),
        period = runif(len, 0, 2 * pi)
      )
    }
    else if (kern == "RQ") {
      hp <- tibble::tibble(
        ID = list_ID,
        variance = runif(len, 1, 5),
        lengthscale = runif(len, 1, 5),
        scale = runif(len, 1, 5)
      )
    }
  }
  return(hp)
}
