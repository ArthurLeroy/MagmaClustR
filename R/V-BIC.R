#' Bayesian Information Criterion
#'
#' @param hp_0 A tibble, data frame or name vector of hyper-parameters for your initial hyperparameters.
#' @param hp_k A tibble, data frame or name vector of hyper-parameters at corresponding variations.
#' @param hp_i A tibble, data frame or name vector of hyper-parameters at corresponding inputs.
#' @param db A tibble containing values we want to compute logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param kern_i Kernel used to compute the covariance matrix of individuals GP at corresponding inputs (Psi_i).
#' @param kern_0 Kernel used to compute the covariance matrix of the mean GP at corresponding inputs (K_0).
#' @param mu_k_param parameters of the variational distributions of mean GPs (mu_k).
#' @param m_k prior value of the mean parameter of the mean GPs (mu_k). Length = 1 or nrow(db).
#' @param clust A Boolean. If True [...] if False [...].
#'
#' @return Variational BIC
#' @export
#'
#' @examples
BIC = function(hp_0, hp_i, hp_k, db, kern_0, kern_i, m_k, mu_k_param, clust = T)
{
  list_clust = mu_k_param$mean %>% names
  nb_indiv = db %>% dplyr::pull(.data$ID) %>% dplyr::n_distinct

  nb_hp_i = hp_i %>% unlist %>% dplyr::n_distinct
  nb_hp_k = hp_k %>% unlist %>% dplyr::n_distinct
  nb_clust = ifelse(clust, hp_k$pi %>% dplyr::n_distinct, 1)
  size_mu = db %>% dplyr::pull(.data$Input) %>% dplyr::n_distinct
  ## Penalty terms for the BIC
  pen =  0.5 * (nb_hp_i + nb_hp_k + nb_clust - 1)*log(nb_indiv)

  floop = function(i)
  {
    t_i = db %>% dplyr::filter(.data$ID == i) %>% dplyr::pull(.data$Input)
    input_t_i = paste0('X', t_i)
    y_i = db %>% dplyr::filter(.data$ID == i) %>% dplyr::pull(.data$Output)
    inv_i = kern_to_inv(t_i, kern_i, hp_i)

    if(clust)
    {
      sum_LL = 0
      LL_Z = 0
      for(k in list_clust)
      {
        mu_k_ti = mu_k_param$mean[[k]] %>% dplyr::filter(.data$Input %in% t_i) %>% dplyr::pull(.data$Output)
        pi_k = hp_k$pi[[k]]
        tau_i_k = mu_k_param$tau_i_k[[k]][[i]]
        cov_k_hat = mu_k_param$cov[[k]][input_t_i, input_t_i]

        sum_LL = sum_LL + tau_i_k * (dmnorm(y_i, mu_k_ti, inv_i, log = T) - 0.5 * sum(inv_i * cov_k_hat))
        LL_Z = LL_Z + tau_i_k * ifelse((tau_i_k == 0|(pi_k == 0)), 0, log(pi_k/ tau_i_k)) # To avoid log(0)
      }
    }
    else
    {
      mu_0 = mu_k_param$mean %>% dplyr::filter(.data$Timestamp %in% t_i) %>% dplyr::pull(.data$Output)
      cov_k_hat = mu_k_param$cov[input_t_i, input_t_i]

      sum_LL = dmnorm(y_i, mu_0, inv_i, log = T) - 0.5 * sum(inv_i * cov_k_hat)
      LL_Z = 0
    }
    #paste0('LL = ', sum_LL, '|| LL_Z = ', LL_Z) %>% print()
    return( sum_LL + LL_Z )
  }
  LL_i = sapply(unique(db$ID), floop) %>% sum

  t = db %>% dplyr::pull(.data$Input) %>% unique
  input_t = paste0('X', t)
  if(clust)
  {
    LL_k = 0
    for(k in list_clust)
    {
      mu_k_t  = mu_k_param$mean[[k]] %>% dplyr::filter(.data$Input %in% t) %>% dplyr::pull(.data$Output)
      mean_k = rep(m_k[[1]], length(t))
      inv_k = kern_to_inv(t  , kern_0, hp_k)
      cov_k_hat = mu_k_param$cov[[k]][input_t, input_t] #+ diag(0.01, size_mu) # jitter term
      diag = diag(cov_k_hat) %>% mean %>% `*`(0) %>% diag(size_mu) # jitter term for numerical problems

      LL_k = LL_k + dmnorm(mu_k_t, mean_k, inv_k, log = T) - 0.5 * sum(inv_k * cov_k_hat) +
        0.5 * (size_mu +  size_mu * log(2*pi) + maotai::pdeterminant(cov_k_hat + diag)$modulus )
    }
  }
  else
  {
    mu_0 = mu_k_param$mean %>% dplyr::filter(.data$Input %in% t) %>% dplyr::pull(.data$Output)
    mean_0 = rep(m_k[[1]], length(t))
    pen_diag = sapply(hp$theta_i, function(x) x[[3]]) %>% mean
    inv_0 = kern_to_inv(t , kern_0, hp_0, pen_diag)
    cov_hat = mu_k_param$cov[input_t, input_t] #+ diag(0.01, size_mu)
    diag = diag(cov_hat) %>% mean %>% `*`(0) %>% diag(size_mu) # jitter term for numerical problems

    LL_k = dmnorm(mu_0, mean_0, inv_0, log = T) - 0.5 * sum(inv_0 * cov_hat) +
      0.5 * (size_mu +  size_mu * log(2*pi) + maotai::pdeterminant(cov_hat + diag)$modulus )
  }
  print(c(LL_i, LL_k, - pen))
  return(LL_i + LL_k - pen)
}

