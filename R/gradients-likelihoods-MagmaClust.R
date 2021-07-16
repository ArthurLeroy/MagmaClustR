#' Gradient multi-Gaussian Process
#'
#' @param hp A tibble, data frame or name vector of hyper-parameters.
#' @param db A tibble containing values we want to compute logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mu_k_param List of parameters for the K mean Gaussian processes.
#' @param kern A kernel function used to compute the covariance matrix at corresponding timestamps.
#'
#' @return The gradient of multi gaussian processes for clustering
#' @export
#'
#' @examples
gr_clust_multi_GP = function(hp, db, mu_k_param, kern)
{
  list_hp = names(hp)

  names_k = mu_k_param$mean %>% names()
  t_i = db$Input
  y_i = db$Output
  i = unique(db$ID)


  corr1 = 0
  corr2 = 0
  for(k in seq_len(length(names_k)))
  {
    tau_i_k = mu_k_param$tau_i_k[[k]][[i]]
    mean_mu_k = mu_k_param$mean[[k]] %>% dplyr::filter(.data$Input %in% t_i) %>% dplyr::pull(.data$Output)
    corr1 = corr1 + tau_i_k * mean_mu_k
    corr2 = corr2 + tau_i_k * ( mean_mu_k %*% t(mean_mu_k) + mu_k_param$cov[[k]][paste0('X', t_i), paste0('X', t_i)] )
  }

  inv = kern_to_inv(t_i, kern, hp)
  prod_inv = inv %*% y_i
  common_term = (prod_inv - 2 * inv %*% corr1) %*% t(prod_inv)  + inv %*% ( corr2 %*% inv - diag(1, length(t_i)) )

  floop = function(deriv){
    (- 1/2 * (common_term %*% kern_to_cov(t_i, kern, hp, deriv))) %>%
      diag() %>%
      sum() %>%
      return()
  }
  sapply(list_hp, floop) %>%
    return()
}

#' Modified common Gaussian Process for each variations.
#'
#' @param hp A tibble, data frame or name vector of hyper-parameters.
#' @param db A tibble containing values we want to compute logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mean mean of the GP at union of observed input
#' @param kern A kernel function used to compute the covariance matrix at corresponding timestamps.
#' @param new_cov posterior covariance matrix of the mean GP (mu_0).
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @returnhe gradient of modified gaussian processes for clustering
#' @export
#'
#' @examples
gr_GP_mod_common_hp_k = function(hp, db, mean, kern, new_cov, pen_diag = NULL)
{
  # list_ID_k = names(db)
  # #t_k = db[[1]] %>% pull(Timestamp)
  # t_k = dplyr::select(-.data$ID) %>% dplyr::slice(1) %>% dplyr::pull(Input)
  # inv =  kern_to_inv(t_k, kern, hp, pen_diag)
  #
  #
  #   y_k = db[[k]] %>% pull(Output)
  #
  #   prod_inv = inv %*% (y_k - mean[[k]])
  #   cste_term = prod_inv %*% t(prod_inv) + inv %*% ( new_cov[[k]] %*% inv - diag(1, length(t_k)) )


  ##############################################
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
