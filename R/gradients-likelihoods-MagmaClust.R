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
#' @param kern A kernel function used to compute the covariance matrix at corresponding timestamps.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#' @param mu_k_param List of parameters for the K mean Gaussian processes.
#'
#' @returnhe gradient of modified gaussian processes for clustering
#' @export
#'
#' @examples
gr_clust_multi_GP_common_hp_i = function(hp, db, mu_k_param, kern, pen_diag = NULL)
{

}
