#' Squared Exponential Kernel
#'
#' @param x a vector/matrix
#' @param y a vector/matrix
#' @param hp a tibble of the hyperameters
#'
#' @return The value of dot production <f(t1),f(t2)> computed from the kernel
#'
#' @examples
#' se_kernel(c(1, 0), c(0, 1), tibble::tibble(variance = 1, lengthscale = 0.5))
se_kernel <- function(x, y, hp) {
  distance <- sum((x - y)^2)
  kern <- exp(hp$variance - exp(-hp$lengthscale) * 0.5 * distance)

  attr(kern, "derivative_variance") <- function(x_1 = x, y_1 = y, hp_1 = hp) {
    exp(hp_1$variance - exp(-hp_1$lengthscale / 2 * as.double(stats::dist(rbind(x_1, y_1)))))
  }

  attr(kern, "derivative_lengthscale") <- function(x_1 = x, y_1 = y, hp_1 = hp) {
    distance <- as.double(stats::dist(rbind(x_1, y_1)))

    grad <- 0.5 * exp(-hp_1$lengthscale) * distance
    derivative_2 <- (exp(hp_1$variance) * grad * exp(-grad))

    derivative_2 %>% return()
  }

  return(kern)
}


#' Periodic Kernel
#'
#' @param x a vector/matrix
#' @param hp a tibble of the hyperparameters
#' @param y a vector/matrix
#'
#' @return The value of dot production <f(t1),f(t2)> computed from the kernel
#'
#' @examples
#' kernel_period(c(1, 0), c(0, 1), tibble::tibble(variance = 1, lengthscale = 0.5, period = 2))
perio_kernel <- function(x, y, hp) {
  distance <- abs(x - y) %>% sum()
  kern <- exp(hp$variance - exp(-hp$lengthscale) * 2 * sin(pi * distance / hp$period)^2)

  attr(kern, "variance") <- function(x_1 = x, y_1 = y, hp_1 = hp) {
    distance <- as.double(stats::dist(rbind(x_1, y_1)))
    grad <- pi * distance / hp_1$period

    (2 * hp_1$variance * exp(-2 * sin(grad)^2 * exp(-hp_1$lengthscale))) %>% return()
  }


  attr(kern, "period") <- function(x_1 = x, y_1 = y, hp_1 = hp) {
    distance <- as.double(stats::dist(rbind(x_1, y_1)))
    grad <- pi * distance / hp_1$period

    (4 * pi * distance * sin(grad) * cos(grad) * exp(hp_1$variance - 2 *
      exp(-hp_1$lengthscale) * sin(grad)^2 - hp_1$period * hp_1$lengthscale)) %>% return()
  }

  attr(kern, "lengthscale") <- function(x_1 = x, y_1 = y, hp_1 = hp) {
    distance <- as.double(stats::dist(rbind(x_1, y_1)))
    grad <- pi * distance / hp_1$period

    ((4 * sin(grad)^2 * exp(hp_1$variance - 2 * sin(grad)^2 * exp(-hp_1$lengthscale))) / hp_1$lengthscale^3) %>% return()
  }

  kern %>% return()
}

#' Rational Quadratic Kernel


#' @param x a vector/matrix
#' @param y a vector/matrix
#' @param hp a tibble of the hyperparameters
#'
#' @return The value of dot production <f(t1),f(t2)> computed from the kernel
#' @export
#'
#' @examples
#' kernel_quad(c(1, 0), c(0, 1), tibble::tibble(variance = 1, lengthscale = 0.5, scale = 3))
rq_kernel <- function(x, y, hp) {
  distance <-  sum((x - y)^2)
  kern <- (exp(hp$variance) * (1 + exp(-hp$lengthscale) / (2 * hp$scale) * distance)^(-hp$scale))

  attr(kern, "variance") <- function(x_1 = x, y_1 = y, hp_1 = hp) {
    distance <- as.double(stats::dist(rbind(x_1, y_1)))
    grad <- (1 + distance * exp(-hp_1$lengthscale) / (2 * hp_1$scale))
    (2 * hp_1$variance * grad^(-hp_1$scale)) %>% return()
  }

  attr(kern, "scale") <- function(x_1 = x, y_1 = y, hp_1 = hp) {
    distance <- as.double(stats::dist(rbind(x_1, y_1)))
    grad <- (1 + distance * exp(-hp_1$lengthscale) / (2 * hp_1$scale))
    (exp(hp_1$variance) * grad^(-hp_1$scale) * (distance * exp(-hp_1$lengthscale) / (2 * hp_1$scale * grad) - log(grad))) %>% return()
  }

  attr(kern, "lengthscale") <- function(x_1 = x, y_1 = y, hp_1 = hp) {
    distance <- as.double(stats::dist(rbind(x_1, y_1)))
    grad <- (1 + distance * exp(-hp_1$lengthscale) / (2 * hp_1$scale))
    (exp(hp_1$variance) * distance * grad^(-hp_1$scale - 1) / (hp_1$lengthscale^3)) %>% return()
  }

  kern %>% return()
}
