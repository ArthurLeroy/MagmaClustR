#' Gradient Gaussian Process
#'
#' @param hp Set of hyper-parameter.
#' @param db Full database with all individuals. Columns required : ID, Input, Output
#' @param mean Mean of your Gaussian Process.
#' @param kern Kernel used to compute the covariance matrix.
#' @param new_cov Posterior covariance matrix of the mean GP (mu_0).
#'
#' @return
#' @export
#'
#' @examples
gr_GP = function(hp, db, mean, kern, new_cov)
{
  list_hp = names(hp)
  output = db$Output
  ## Extract the reference Input
  input = db$Input
  ## Extract the input variables (reference Input + Covariates)
  inputs = db %>% dplyr::select(- .data$Output)

  cov = kern_to_cov(inputs, kern, hp) + new_cov
  inv = tryCatch(solve(cov), error = function(e){MASS::ginv(cov)})
  prod_inv = inv %*% (output - mean)
  common_term = prod_inv %*% t(prod_inv) - inv

  floop = function(deriv){
    (- 1/2 * (common_term %*% kern_to_cov(inputs, kern, hp, deriv))) %>%
      diag() %>%
      sum() %>%
      return()
  }
  sapply(list_hp, floop) %>%
    return()
}

#' Gradient Gaussian Process modif
#'
#' @param hp set of hyper-parameter
#' @param db full database with all individuals. Columns required : ID, input, Output
#' @param mean mean of your Gaussian Process
#' @param kern kernel used to compute the covariance matrix
#' @param new_cov posterior covariance matrix of the mean GP (mu_0).
#' @param pen_diag value of the penalization of the diagonal
#'
#' @return Gradient of the Gaussian Process modified
#' @export
#'
#' @examples
#' TRUE
gr_GP_mod = function(hp, db, mean, kern, new_cov, pen_diag = NULL)
{
  list_hp = names(hp)
  output = db$Output
  ## Extract the reference Input
  input = db$Input
  ## Extract the input variables (reference Input + Covariates)
  inputs = db %>% dplyr::select(- .data$Output)

  inv <- kern_to_inv(inputs, kern, hp, pen_diag)
  prod_inv = inv %*% (output - mean)
  common_term = prod_inv %*% t(prod_inv) +
    inv %*% ( new_cov %*% inv - diag(1, length(input)) )

  floop = function(deriv){
    (- 1/2 * (common_term %*% kern_to_cov(inputs, kern, hp, deriv))) %>%
      diag() %>%
      sum() %>%
      return()
  }
  sapply(list_hp, floop) %>%
    return()
}

#' Gradient common Gaussian Process
#'
#' @param hp set of hyper-parameter
#' @param db full database with all individuals. Columns required : ID, input, Output
#' @param mean mean of the GP at union of observed input
#' @param kern kernel used to compute the covariance matrix at corresponding inputs
#' @param new_cov posterior covariance matrix of the mean GP (mu_0).
#' @param pen_diag A jitter
#'
#' @return Gradient of the common Gaussian Process
#' @export
#'
#' @examples
#' TRUE
gr_GP_mod_common_hp = function(hp, db, mean, kern, new_cov, pen_diag)
{
  list_hp = names(hp)
  ## Initialise the the sum of the gradients and inputs
  gr_hp = rep(0, length(hp)) %>% `names<-`(hp %>% names())
  input_i_old = NULL

  ## Loop over individuals to compute the sum of log-Likelihoods
  for(i in unique(db$ID)){
    ## Extract the i-th specific reference inputs
    input_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Input)
    ## Extract the Output values associated with the i-th specific inputs
    output_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Output)
    ## Extract the mean values associated with the i-th specific inputs
    mean_i = mean %>%
      dplyr::filter(.data$Input %in% input_i) %>%
      dplyr::pull(.data$Output)
    ## Extract the covariance values associated with the i-th specific inputs
    new_cov_i = new_cov[as.character(input_i), as.character(input_i)]

    ## Update the inverse cov matrix only if necessary (if different inputs)
    if (!identical(input_i, input_i_old)) {
      ## Extract the i-th specific inputs (reference + covariates)
      inputs_i = db %>%
        dplyr::filter(.data$ID == i) %>%
        dplyr::select(- c(.data$ID, .data$Output))
      inv <- kern_to_inv(inputs_i, kern, hp, pen_diag)
    }

    ## Compute the terms common to all hyper-parameters
    prod_inv = inv %*% (output_i - mean_i)
    common_term = prod_inv %*% t(prod_inv) +
      inv %*% (new_cov_i  %*% inv - diag(1, length(input_i)) )
    ## Loop over the derivatives of hyper-parameters for computing the gradient
    floop = function(deriv){
      (- 1/2 * (common_term %*% kern_to_cov(inputs_i, kern, hp, deriv))) %>%
        diag() %>%
        sum() %>%
        return()
    }
    ## Add a new term to the sum of gradients
    gr_hp = gr_hp + sapply(list_hp, floop)
    ## Keep track of the reference inputs
    input_i_old <- input_i
  }
  return(gr_hp)
}
