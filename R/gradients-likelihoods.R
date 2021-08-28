#' Gradient of the logLikelihood of a Gaussian Process
#'
#' @param hp A tibble, data frame or named vector containing hyper-parameters.
#' @param db A tibble containing the values we want to compute the logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mean A vector, specifying the mean of the GP at the reference inputs.
#' @param kern A kernel function.
#' @param post_cov (optional) A matrix, corresponding to covariance parameter of
#'    the hyper-posterior. Used to compute the hyper-prior distribution of a new
#'    individual in Magma.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A named vector, corresponding to the value of the hyper-parameters
#'    gradients for the Gaussian log-Likelihood (where the covariance can be the
#'    sum of the individual and the hyper-posterior's mean process covariances).
#'
#' @examples
#' db <- tibble::tibble(Input = 1:5, Output = 2:6)
#' mean <- rep(0, 5)
#' hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
#' post_cov = kern_to_cov(1:5, 'SE', hp)
#' MagmaClustR:::gr_GP(hp, db, mean, "SE", post_cov, 0.001)
gr_GP <- function(hp, db, mean, kern, post_cov, pen_diag) {
  list_hp <- names(hp)
  output <- db$Output
  ## Extract the reference Input
  input <- db$Input
  ## Extract the input variables (reference Input + Covariates)
  inputs <- db %>% dplyr::select(-.data$Output)

  cov <- kern_to_cov(inputs, kern, hp) + post_cov
  diag <- diag(x = pen_diag, ncol = ncol(cov), nrow = nrow(cov))
  inv <- tryCatch((cov + diag) %>% chol() %>% chol2inv(),
    error = function(e) {
      MASS::ginv(cov + diag)
    }
  )

  ## Compute the term common to all partial derivatives
  prod_inv <- inv %*% (output - mean)
  common_term <- prod_inv %*% t(prod_inv) - inv

  ## Loop over the derivatives of hyper-parameters for computing the gradient
  floop <- function(deriv) {
    (-0.5 * (common_term %*% kern_to_cov(inputs, kern, hp, deriv))) %>%
      diag() %>%
      sum() %>%
      return()
  }
  sapply(list_hp, floop) %>%
    return()
}


#' Gradient of the modified logLikelihood for GPs in Magma
#'
#' @param hp A tibble, data frame or named vector containing hyper-parameters.
#' @param db A tibble containing the values we want to compute the logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mean A vector, specifying the mean of the GPs at the reference inputs.
#' @param kern A kernel function.
#' @param post_cov A matrix, covariance parameter of the hyper-posterior.
#'    Used to compute the correction term.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A named vector, corresponding to the value of the hyper-parameters
#'    gradients for the modified Gaussian log-Likelihood involved in Magma.
#'
#' @examples
#' db <- tibble::tibble(Input = 1:5, Output = 2:6)
#' mean <- rep(0, 5)
#' hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
#' post_cov = kern_to_cov(1:5, 'SE', hp)
#' MagmaClustR:::gr_GP_mod(hp, db, mean, "SE", post_cov, 0.001)
gr_GP_mod <- function(hp, db, mean, kern, post_cov, pen_diag) {
  list_hp <- names(hp)
  output <- db$Output
  ## Extract the reference Input
  input <- db$Input
  ## Extract the input variables (reference Input + Covariates)
  inputs <- db %>% dplyr::select(-.data$Output)

  inv <- kern_to_inv(inputs, kern, hp, pen_diag)
  ## Compute the term common to all partial derivatives
  prod_inv <- inv %*% (output - mean)
  common_term <- prod_inv %*% t(prod_inv) +
    inv %*% (post_cov %*% inv - diag(1, length(input)))

  ## Loop over the derivatives of hyper-parameters for computing the gradient
  floop <- function(deriv) {
    (-1 / 2 * (common_term %*% kern_to_cov(inputs, kern, hp, deriv))) %>%
      diag() %>%
      sum() %>%
      return()
  }
  sapply(list_hp, floop) %>%
    return()
}

#' Gradient of the modified logLikelihood with common HPs for GPs in Magma
#'
#' @param hp A tibble or data frame containing hyper-parameters for all
#'    individuals.
#' @param db A tibble containing the values we want to compute the logL on.
#'    Required columns: ID, Input, Output. Additional covariate columns are
#'    allowed.
#' @param mean A vector, specifying the mean of the GPs at the reference inputs.
#' @param kern A kernel function.
#' @param post_cov A matrix, covariance parameter of the hyper-posterior.
#'    Used to compute the correction term.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A named vector, corresponding to the value of the hyper-parameters'
#'    gradients for the modified Gaussian log-Likelihood involved in Magma with
#'    the 'common HP' setting.
#'
#' @examples
#' db <- simu_db(N = 10, common_input = TRUE)
#' mean <- tibble::tibble(Input = unique(db$Input), Output = 0)
#' hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
#' post_cov <- kern_to_cov(unique(db$Input), 'SE', hp)
#' MagmaClustR:::gr_GP_mod_common_hp(hp, db, mean, "SE", post_cov, 0.001)
gr_GP_mod_common_hp <- function(hp, db, mean, kern, post_cov, pen_diag) {

  list_hp <- names(hp)

  ## Loop over individuals to compute the sum of log-Likelihoods
  funloop <- function(i) {
    ## Extract the i-th specific reference Input
    input_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Input)
    ## Extract the i-th specific inputs (reference + covariates)
    inputs_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(- c(.data$ID, .data$Output))
    ## Extract the i-th specific Inputs and Output
    output_i = db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Output)
    ## Extract the mean values associated with the i-th specific inputs
    mean_i = mean %>%
      dplyr::filter(.data$Input %in% input_i) %>%
      dplyr::pull(.data$Output)
    ## Extract the covariance values associated with the i-th specific inputs
    post_cov_i = post_cov[as.character(input_i), as.character(input_i)]

    ## Compute the term common to all partial derivatives
    inv <- kern_to_inv(inputs_i, kern, hp, pen_diag)
    prod_inv <- inv %*% (output_i - mean_i)
    common_term <- prod_inv %*% t(prod_inv) +
      inv %*% (post_cov_i %*% inv - diag(1, length(input_i)))
    ## Loop over the derivatives of hyper-parameters for computing the gradient
    floop <- function(deriv) {
      (- 0.5 * (common_term %*% kern_to_cov(inputs_i, kern, hp, deriv))) %>%
        diag() %>%
        sum() %>%
        return()
    }
    sapply(list_hp, floop) %>%
      return()
  }
  sapply(unique(db$ID), funloop) %>%
    rowSums() %>%
    return()
}
