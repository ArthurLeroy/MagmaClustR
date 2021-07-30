#' Gradient multi-Gaussian Process
#'
#' @param hp A tibble, data frame or name vector of hyper-parameters.
#' @param db A tibble containing values we want to compute logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mu_k_param List of parameters for the K mean Gaussian processes.
#' @param kern A kernel function used to compute the covariance matrix at corresponding timestamps.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return The gradient of multi gaussian processes for clustering
#' @export
#'
#' @examples
gr_clust_multi_GP = function(hp, db, mu_k_param, kern, pen_diag)
{
  #browser()
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
    corr2 = corr2 + tau_i_k * ( mean_mu_k %*% t(mean_mu_k) + mu_k_param$cov[[k]][as.character(t_i), as.character(t_i)] )
  }

  inv = kern_to_inv(t_i, kern, hp, pen_diag)
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


#' Modified common Gaussian Process for each cluster
#'
#' @param hp A tibble, data frame or name vector of hyper-parameters.
#' @param db A tibble containing values we want to compute logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param kern A kernel function used to compute the covariance matrix at corresponding timestamps.
#' @param mean A vector, specifying the mean of the GP at the reference inputs.
#' @param new_cov posterior covariance matrix of the mean GP (mu_0). Used to compute correction term (cor_term)#'
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return gradient of modified gaussian processes for clustering
#' @export
#'
#' @examples
gr_GP_mod_common_hp_k = function(hp, db, mean, kern, new_cov, pen_diag = NULL)
{
  list_ID_k = names(db)
  t_k = db[[1]] %>% dplyr::pull(.data$Input)
  inv =  kern_to_inv(t_k, kern, hp, pen_diag)

  g_1 = 0
  g_2 = 0
  g_3 = 0

  for(k in list_ID_k)
  {
    y_k = db[[k]] %>% dplyr::pull(.data$Output)

    prod_inv = inv %*% (y_k - mean[[k]])
    common_term = prod_inv %*% t(prod_inv) + inv %*% ( new_cov[[k]] %*% inv - diag(1, length(t_k)) )

    floop = function(deriv){
      (- 1/2 * (common_term %*% kern_to_cov(t_k, kern, hp, deriv))) %>%
        diag() %>%
        sum() %>%
        return()
    }

     sapply(list_ID_k, floop) %>%
       return()
  }

  #sapply(list_ID_k, floop) %>%
  #  return()

}


#' Modified common Gaussian Process for each variations.
#'
#' @param hp A tibble, data frame or name vector of hyper-parameters.
#' @param db A tibble containing values we want to compute logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param kern A kernel function used to compute the covariance matrix at corresponding timestamps.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#' @param mu_k_param List of parameters for the K mean Gaussian processes.
#'
#' @return gradient of modified gaussian processes for clustering
#' @export
#'
#' @examples
gr_clust_multi_GP_common_hp_i = function(hp, db, mu_k_param, kern, pen_diag = NULL)
{
  names_k = mu_k_param$mean %>% names()
  t_i_old = NULL

  for(i in unique(db$ID))
  {
    t_i = db %>% dplyr::filter(.data$ID == i) %>% dplyr::pull(.data$Input)
    input_i = as.character(t_i)
    y_i = db %>% dplyr::filter(.data$ID == i) %>% dplyr::pull(.data$Output)

    corr1 = 0
    corr2 = 0

    for(k in seq_len(length(names_k)))
    {
      tau_i_k = mu_k_param$tau_i_k[[k]][[i]]
      mean_mu_k = mu_k_param$mean[[k]] %>% dplyr::filter(.data$Input %in% t_i) %>% dplyr::pull(.data$Output)
      corr1 = corr1 + tau_i_k * mean_mu_k
      corr2 = corr2 + tau_i_k * ( mean_mu_k %*% t(mean_mu_k) + mu_k_param$cov[[k]][input_i, input_i] )
    }

    if( !identical(t_i, t_i_old) )
    { ## We update the inverse cov matrix only if necessary (if different timestamps)
      inv = kern_to_inv(t_i, kern, hp, pen_diag)
    }

    prod_inv = inv %*% y_i
    common_term = (prod_inv - 2 * inv %*% corr1) %*% t(prod_inv)  +
      inv %*% ( corr2 %*% inv - diag(1, length(t_i)) )

    floop = function(deriv){
      (- 1/2 * (common_term %*% kern_to_cov(t_i, kern, hp, deriv))) %>%
        diag() %>%
        sum() %>%
        return()
    }

    sapply(names_k, floop) %>%
      return()
  }

}
