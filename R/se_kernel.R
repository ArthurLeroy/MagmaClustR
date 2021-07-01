#' Squared Exponential Kernel
#'
#' @param x a vector/matrix
#' @param y a vector/matrix
#' @param hp a tibble of the hyperameters
#'
#' @return The value of dot production <f(t1),f(t2)> computed from the kernel
#' @export
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
