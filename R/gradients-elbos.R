#' Gradient multi-Gaussian Process
#'
#' @param hp A tibble, data frame or named vector containing hyper-parameters.
#' @param db A tibble containing the values we want to compute the elbo on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mu_k_param List of parameters for the K mean Gaussian processes.
#' @param kern A kernel function.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return The gradient of multi gaussian processes for clustering
#'
#' @examples
#' TRUE
gr_clust_multi_GP = function(hp, db, mu_k_param, kern, pen_diag)
{
  #browser()
  list_hp = names(hp)

  names_k = mu_k_param$mean %>% names()
  t_i = db$Input
  y_i = db$Output
  #inputs = db %>% dplyr::select(-.data$Output)
  i = unique(db$ID)


  corr1 = 0
  corr2 = 0

  for(k in (names_k) )
  {
    hp_mixture = mu_k_param$hp_mixture[k][i,]
    mean_mu_k = mu_k_param$mean[[k]] %>% dplyr::filter(.data$Input %in% t_i) %>%
      dplyr::pull(.data$Output)
    corr1 = corr1 + as.double(hp_mixture) * mean_mu_k
    corr2 = corr2 + as.double(hp_mixture) *
      ( mean_mu_k %*% t(mean_mu_k) +
          mu_k_param$cov[[k]][as.character(t_i), as.character(t_i)] )
  }

  inv = kern_to_inv(t_i, kern, hp, pen_diag)
  prod_inv = inv %*% y_i
  common_term = (prod_inv - 2 * inv %*% corr1) %*% t(prod_inv)  +
    inv %*% ( corr2 %*% inv - diag(1, length(t_i)) )

  ## Loop over the derivatives of hyper-parameters for computing the gradient
  floop = function(deriv){
    (- 1/2 * (common_term %*% kern_to_cov(t_i, kern, hp, deriv))) %>%
      diag() %>%
      sum() %>%
      return()
  }
  #browser()
  sapply(list_hp, floop) %>%
    return()
}


#' Modified common Gaussian Process for each cluster
#'
#' @param hp A tibble, data frame or named vector containing hyper-parameters..
#' @param db A tibble containing the values we want to compute the elbo on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param kern A kernel function
#' @param mean A list of the k means of the GP at union of observed timestamps.
#' @param post_cov list of the k posterior covariance of the mean GP (mu_k).
#' Used to compute correction term (cor_term)
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return The value of the modified Gaussian log-likelihood for
#' the sum of the k mean GPs with same HPs.
#'
#' @examples
#' TRUE
gr_GP_mod_common_hp_k = function(hp, db, mean, kern, post_cov, pen_diag = NULL)
{
  #browser()
  if('ID' %in% names(hp)){
    hp <- hp %>% dplyr::select(- .data$ID)
  }

  list_ID_k = names(db)
  list_hp <- names(hp)

  ## Extract the i-th specific reference Input
  input_k <- db[[1]] %>%
   dplyr::pull(.data$Input)
  ## Extract the i-th specific inputs (reference + covariates)
  inputs_k <- db[[1]] %>%
    dplyr::select(- .data$Output)
  if('ID' %in% names(db[[1]])){
    inputs_k <- inputs_k %>% dplyr::select(- .data$ID)
  }



  inv =  kern_to_inv(inputs_k, kern, hp, pen_diag)

  ## Loop over individuals to compute the sum of log-Likelihoods
  funloop <- function(k) {
    ## Extract the i-th specific Inputs and Output
    output_k = db[[k]] %>%
      dplyr::pull(.data$Output)
    ## Extract the mean values associated with the i-th specific inputs
    mean_k = mean[[k]]
    ## Extract the covariance values associated with the i-th specific inputs
    post_cov_k = post_cov[[k]]

    prod_inv = inv %*% (output_k - mean_k)
    common_term = prod_inv %*% t(prod_inv) +
      inv %*% ( post_cov_k %*% inv - diag(1, length(input_k)) )

    floop = function(deriv){
      (- 1/2 * (common_term %*% kern_to_cov(inputs_k, kern, hp, deriv))) %>%
        diag() %>%
        sum() %>%
        return()
    }

     sapply(list_hp, floop) %>%
       return()
  }

  sapply(list_ID_k, funloop) %>% rowSums() %>%
    return()

}


#' Modified common Gaussian Process for each variations.
#'
#' @param hp A tibble, data frame or name vector of hyper-parameters.
#' @param db A tibble containing values we want to compute elbo on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param kern A kernel function used to compute the covariance matrix at
#' corresponding timestamps.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#' @param mu_k_param List of parameters for the K mean Gaussian processes.
#'
#' @return The value of the modified Gaussian log-likelihood for
#' one GP as it appears in the model.
#'
#' @examples
#' TRUE
gr_clust_multi_GP_common_hp_i = function(hp, db, mu_k_param, kern, pen_diag = NULL)
{
  list_hp <- names(hp)
  names_k = mu_k_param$mean %>% names()

  ## Loop over individuals to compute the sum of log-Likelihoods
  funloop <- function(i){
    #browser()
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

    corr1 = 0
    corr2 = 0

    for(k in (names_k) )
    {
      #browser()
      ## Extract the covariance values associated with the i-th specific inputs
      post_cov_i = mu_k_param$cov[[k]][as.character(input_i), as.character(input_i)]

      hp_mixture = mu_k_param$hp_mixture[k][i,]
      mean_mu_k = mu_k_param$mean[[k]] %>% dplyr::filter(.data$Input %in% input_i) %>%
        dplyr::pull(.data$Output)
      corr1 = corr1 + as.double(hp_mixture) * mean_mu_k
      corr2 = corr2 + as.double(hp_mixture) *
        ( mean_mu_k %*% t(mean_mu_k) + post_cov_i )
    }

    inv = kern_to_inv(inputs_i, kern, hp, pen_diag)

    prod_inv = inv %*% output_i
    common_term = (prod_inv - 2 * inv %*% corr1) %*% t(prod_inv)  +
      inv %*% ( corr2 %*% inv - diag(1, length(input_i)) )

    floop = function(deriv){
      (- 1/2 * (common_term %*% kern_to_cov(inputs_i, kern, hp, deriv))) %>%
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