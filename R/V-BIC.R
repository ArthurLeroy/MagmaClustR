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


#' Selection of the model
#'
#' @param db A tibble containing values we want to compute logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param k_grid A vector, indicating the grid of additional reference
#'    inputs on which the mean process' hyper-posterior should be evaluated.
#' @param ini_hp_k named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_0}, the mean process' kernel. The
#'    columns/elements should be named according to the hyper-parameters
#'    that are used in \code{kern_0}.
#' @param ini_hp_i A tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}, the individual processes' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each individual/task. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_i}.
#' @param kern_0 A kernel function, associated with the mean GP.
#' @param kern_i A kernel function, associated with the individual GPs.
#' @param ini_tau_i_k initial values of probabiliy to belong to each cluster for each individuals.
#' @param common_hp_k boolean indicating whether hp are common among mean GPs (for each mu_k)
#' @param common_hp_i boolean indicating whether hp are common among individual GPs (for each y_i)
#' @param plot boolean indicating whether you want to plot or not
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#'
#' @return TRUE
#' @export
#'
#' @examples
#' TRUE
model_selection = function(db, k_grid = 1:5, ini_hp_k, ini_hp_i,
                           kern_0, kern_i, ini_tau_i_k = NULL, common_hp_k = T, common_hp_i = T,
                           plot = T, pen_diag)
{
  # floop = function(K)
  # {
  #   print(paste0('K = ', K))
  #   prior_mean_k = rep(0, K) %>% stats::setNames(paste0('K', seq_len(K))) %>% as.list
  #
  #   if(K == 1)
  #   {
  #     model = train_magma(db, 0, ini_hp_k, ini_hp_i, kern_0, kern_i, common_hp_i, pen_diag)
  #   }
  #   else
  #   {
  #     model = train_magma_VEM(db, prior_mean_k, ini_hp_k, ini_hp_i, kern_0, kern_i, ini_tau_i_k, common_hp_k, common_hp_i, pen_diag)
  #   }
  #   model[['BIC']] = BIC(model$hp_0, model$hp_i, hp_k ,db, kern_0, kern_i, prior_mean_k, model$param, K != 1)
  #   ## Comment hp_k ?
  #   return(model)
  # }
  # res = k_grid %>% sapply(floop, simplify = FALSE, USE.NAMES = TRUE) %>% stats::setNames(paste0('K = ', k_grid))
  #
  # db_plot = tibble::tibble('K' = k_grid, 'BIC' = res %>% purrr::map_dbl('BIC'))
  #
  # res$K_max_BIC = db_plot %>% dplyr::filter(BIC == max(BIC)) %>% dplyr::pull(.data$K)
  #
  # # if(plot)
  # # {
  # #   res$plot = ggplot2::ggplot(db_plot, ggplot2::aes(x = K, y = BIC)) + ggplot2::geom_point() +
  # #     ggplot2::geom_line() + ggplot2::theme_classic()
  # # }
  #
  # return(res)
}
