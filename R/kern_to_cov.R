#' Kernel to Covariance
#'
#' @param input tibble of all inputs for each individuals.
#' @param kern indicates which kernel function to use to compute the covariance matrix
#' @param hp the hyperparameters used in your kernel
#'
#' @return list of covariance matrices (1 by individual) of the input vectors according to kernel() function
#' @export
#'
#' @examples
#' x = rbind(c(1,0,1),c(2,1,2),c(1,2,3))
#' kern_to_cov(x,kernel_sqrd_exp, tibble::tibble(sigma = 1, lengthscale = 0.5))
#'
kern_to_cov <- function(input, kern, hp)
{
  x = t(input)

  #transform the batches of input into lists

  list_input = split(x, rep(1:ncol(x), each = nrow(x)) )
  outer(list_input,list_input, Vectorize(function(x,y) kern(x,y, hp )) ) %>% return()
}
