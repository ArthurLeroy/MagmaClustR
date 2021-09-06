#' Log-Likelihood function for clusturing mutli-Gaussian Process
#'
#' @param hp A tibble, data frame or name vector containing hyper-parameters.
#' @param db A tibble containing the values we want to compute the logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mu_k_param List of parameters for the K mean Gaussian processes.
#' @param kern A kernel function used to compute the covariance matrix at corresponding timestamps.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return The value of the modified Gaussian log-likelihood for one GP as it appears in the model.
#' @export
#'
#' @examples
logL_clust_multi_GP = function(hp, db, mu_k_param, kern, pen_diag)
{
  #browser()
  names_k = mu_k_param$mean %>% names()
  t_i = db$Input
  y_i = db$Output
  i = unique(db$ID)

  inv =  kern_to_inv(t_i, kern, hp, pen_diag)

  LL_norm = - dmnorm(y_i, rep(0, length(y_i)), inv, log = T) ## classic gaussian centered loglikelihood

  corr1 = 0
  corr2 = 0

  #floop = function(k){
  for(k in (names_k))
  {
    tau_i_k = mu_k_param$tau_i_k[k][i,]
    mean_mu_k = mu_k_param$mean[[k]] %>% dplyr::filter(.data$Input %in% t_i) %>% dplyr::pull(.data$Output)
    corr1 = corr1 + as.double(tau_i_k) * mean_mu_k
    corr2 = corr2 + as.double(tau_i_k) * ( mean_mu_k %*% t(mean_mu_k) + mu_k_param$cov[[k]][as.character(t_i), as.character(t_i)] )
  }

  #sapply(seq_len(length(names_k)), floop) %>%
  #  return()

  ( LL_norm - y_i %*% inv %*% corr1 + 0.5 * sum(inv * corr2) ) %>% return()
}

#' Modified Gaussian log-likelihood for the sum of the k mean GPs with same HPs
#'
#' @param hp A tibble, data frame or name vector containing hyper-parameters.
#' @param db  A tibble containing values we want to compute logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mean list of the k means of the GP at union of observed timestamps
#' @param kern A kernel function used to compute the covariance matrix at corresponding timestamps.
#' @param post_cov A List of the k posterior covariance of the mean GP (mu_k). Used to compute correction term (cor_term).
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#'
#' @return value of the modified Gaussian log-likelihood for the sum of the k mean GPs with same HPs
#' @export
#'
#' @examples
logL_GP_mod_common_hp_k = function(hp, db, mean, kern, post_cov, pen_diag = NULL)
{
  ## To avoid pathological behaviour of the opm optimization function in rare cases
  #if((hp %>% abs %>% sum) > 20){return(10^10)}
  #browser()

  list_ID_k = names(db)
  #t_k = db[[1]] %>% dplyr::pull(.data$Input)
  t_k = db[[1]] %>%
    dplyr::pull(.data$Input)
  inv =  kern_to_inv(t_k, kern, hp, pen_diag)

  LL_norm = 0
  cor_term = 0

  for(k in list_ID_k)
  {
    y_k = db[[k]] %>% dplyr::pull(.data$Output)

    LL_norm = LL_norm - dmnorm(y_k, rep(mean[[k]], length(t_k)), inv, log = T)
    cor_term = cor_term + 0.5 * (inv * post_cov[[k]]) %>% sum()  ##(0.5 * Trace(inv %*% post_cov))
  }
  return(LL_norm + cor_term)
}

#' Modified Gaussian log-likelihood for for the sum of all indiv with same HPs
#'
#' @param hp A tibble, data frame or name vector containing hyper-parameters.
#' @param db A tibble containing values we want to compute logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mu_k_param List of parameters for the K mean Gaussian processes.
#' @param kern A kernel function used to compute the covariance matrix at corresponding timestamps.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return The value of the modified Gaussian log-likelihood for for the sum of all indiv with same HPs.
#' @export
#'
#' @examples
logL_clust_multi_GP_common_hp_i = function(hp, db, mu_k_param, kern, pen_diag)
{
  ## To avoid pathological behaviour of the opm optimization function in rare cases
  #if((hp %>% abs %>% sum) > 20){return(10^10)}

  names_k = mu_k_param$mean %>% names()
  t = unique(db$Input)

  sum_i = 0
  t_i_old = NULL

  for(i in unique(db$ID))
  {
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

    #t_i = db %>% dplyr::filter(.data$ID == i) %>% dplyr::pull(.data$Input)
    #input_i = as.character(t_i)
    #y_i = db %>% dplyr::filter(.data$ID == i) %>% dplyr::pull(.data$Output)

    corr1 = 0
    corr2 = 0

    for(k in (names_k) )
    {
      ## Extract the covariance values associated with the i-th specific inputs
      post_cov_i = mu_k_param$cov[[k]][as.character(input_i), as.character(input_i)]

      tau_i_k = mu_k_param$tau_i_k[k][i,]
      mean_mu_k = mu_k_param$mean[[k]] %>% dplyr::filter(.data$Input %in% input_i) %>% dplyr::pull(.data$Output)
      corr1 = corr1 + as.double(tau_i_k) * mean_mu_k
      corr2 = corr2 + as.double(tau_i_k) * ( mean_mu_k %*% t(mean_mu_k) + post_cov_i)
    }

    inv = kern_to_inv(inputs_i, kern, hp, pen_diag)

    LL_norm = - dmnorm(output_i, rep(0, length(output_i)), inv, log = T) ## classic gaussian centered loglikelihood

    sum_i = sum_i + LL_norm - output_i %*% inv %*% corr1 + 0.5 * sum(inv * corr2)

  }
  return(sum_i)
}

#' Expectation of joint log-likelihood of the model
#'
#' @param hp_k A tibble, data frame or name vector of hyper-parameters at corresponding clusters.
#' @param hp_i A tibble, data frame or name vector of hyper-parameters at corresponding individuals.
#' @param db A tibble containing values we want to compute logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param kern_i Kernel used to compute the covariance matrix of individuals GP at corresponding inputs (Psi_i).
#' @param kern_0 Kernel used to compute the covariance matrix of the mean GP at corresponding inputs (K_0).
#' @param mu_k_param parameters of the variational distributions of mean GPs (mu_k).
#' @param m_k prior value of the mean parameter of the mean GPs (mu_k). Length = 1 or nrow(db).
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.

#'
#' @return Value of expectation of joint log-likelihood of the model. The function to be maximised in step M.
#' @export
#'
#' @examples
logL_monitoring_VEM = function(hp_k, hp_i, db, kern_i, kern_0, mu_k_param, m_k, pen_diag)
{
  floop = function(k)
  {
    logL_GP_mod(hp_k[hp_k$ID == k,], db = mu_k_param$mean[[k]], mean = m_k[[k]] , kern_0, mu_k_param$cov[[k]], pen_diag) %>%
      return()
  }
  sum_ll_k = sapply(names(m_k), floop) %>% sum()

  floop2 = function(i)
  {
    t_i = db %>% dplyr::filter(.data$ID == i) %>% dplyr::pull(.data$Input)
    logL_clust_multi_GP(hp_i[hp_i$ID == i,], db %>% dplyr::filter(.data$ID == i), mu_k_param, kern_0, pen_diag) %>% return()
  }
  sum_ll_i = sapply(unique(db$ID), floop2) %>% sum()
  return(-sum_ll_k - sum_ll_i)
}

