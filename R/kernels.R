#' Squared Exponential Kernel
#'
#' @param x a vector/matrix
#' @param y a vector/matrix
#' @param hp a tibble of the hyperameters
#'
#' @return The value of dot production <f(t1),f(t2)> computed from the kernel
#'
#' @examples
#' se_kernel(c(1,0),c(0,1), tibble::tibble(sigma = 1, lengthscale = 0.5))
se_kernel <- function(x,y, hp)
{
  distance =  as.double(stats::dist(rbind(x, y)))
  kern <- exp(hp$sigma - exp(-hp$lengthscale )/2 * distance)

  attr(kern,"derivative_sigma") <- function(x_1=x, y_1=y, hp_1=hp){
    exp(hp_1$sigma - exp(-hp_1$lengthscale/2 * as.double(stats::dist(rbind(x_1, y_1)))))
  }


  attr(kern, "derivative_lengthscale") <- function(x_1=x, y_1=y, hp_1=hp){
    distance =  as.double(stats::dist(rbind(x_1, y_1)))

    grad = 0.5 * exp(- hp_1$lengthscale) * distance
    derivative_2 = (exp(hp_1$sigma) * grad * exp(- grad) )

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
#' kernel_period(c(1,0),c(0,1),tibble::tibble(sigma = 1, lengthscale = 0.5, period = 2))
perio_kernel <-function(x,y,hp)
{
  distance =  as.double(stats::dist(rbind(x, y)))
  kern = exp(hp$sigma - exp(- hp$lengthscale) * 2 * sin(pi * distance / hp$period)^2)


  attr(kern, "derivative_sigma") <- function(x_1=x, y_1=y, hp_1=hp){
    distance =  as.double(stats::dist(rbind(x_1, y_1)))
    grad = pi * distance / hp_1$period

    (2*hp_1$sigma * exp(-2 * sin(grad)^2 * exp(- hp_1$lengthscale))) %>% return()}


  attr(kern, "derivative_period") <- function(x_1=x, y_1=y, hp_1=hp){
    distance =  as.double(stats::dist(rbind(x_1, y_1)))
    grad = pi * distance / hp_1$period

    (4*pi*distance* sin(grad) * cos(grad)* exp(hp_1$sigma - 2*
                                                 exp(- hp_1$lengthscale)*sin(grad)^2 - hp_1$period* hp_1$lengthscale) ) %>% return()}


  attr(kern, "derivative_lengthscale") <- function(x_1=x, y_1=y, hp_1=hp){
    distance =  as.double(stats::dist(rbind(x_1, y_1)))
    grad = pi * distance / hp_1$period

    ((4*sin(grad)^2 * exp(hp_1$sigma - 2* sin(grad)^2 * exp(- hp_1$lengthscale)) ) / hp_1$lengthscale^3) %>% return()}


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
#' kernel_quad(c(1,0),c(0,1),tibble::tibble(sigma = 1, lengthscale = 0.5, alpha = 3))

rq_kernel <- function(x,y,hp)
{
  distance =  as.double(stats::dist(rbind(x, y)))
  kern = (exp(hp$sigma)*(1+ exp(-hp$lengthscale)/(2*hp$alpha) * distance)^(-hp$alpha))

  attr(kern, "derivative_sigma") <- function(x_1=x, y_1=y, hp_1=hp){
    distance =  as.double(stats::dist(rbind(x_1, y_1)))
    grad = (1 + distance * exp(-hp_1$lengthscale)/(2* hp_1$alpha) )
    (2 * hp_1$sigma * grad^(-hp_1$alpha)) %>% return()
  }

  attr(kern, "derivative_alpha") <- function(x_1=x, y_1=y, hp_1=hp){
    distance =  as.double(stats::dist(rbind(x_1, y_1)))
    grad = (1 + distance * exp(-hp_1$lengthscale)/(2* hp_1$alpha) )
    (exp(hp_1$sigma) * grad^(-hp_1$alpha) * (distance * exp(-hp_1$lengthscale) /(2*hp_1$alpha * grad) - log(grad))) %>% return() }

  attr(kern, "derivative_lengthscale") <- function(x_1=x, y_1=y, hp_1=hp){
    distance =  as.double(stats::dist(rbind(x_1, y_1)))
    grad = (1 + distance * exp(-hp_1$lengthscale)/(2* hp_1$alpha) )
    (exp(hp_1$sigma) * distance * grad ^(-hp_1$alpha - 1) / (hp_1$lengthscale ^3)) %>% return()}


  kern %>% return()
}
