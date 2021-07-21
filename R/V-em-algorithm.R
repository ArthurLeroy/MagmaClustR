#' Variation of the E-Step of the EM algorithm
#'
#' Expectation step of the variation of the EM algorithm to compute the parameters of the
#' hyper-posterior Gaussian distribution of the mean process in Magma.
#'
#' @param db A tibble or data frame. Columns required: ID, Input, Output.
#'    Additional columns for covariates can be specified.
#' @param kern_0 A kernel function, associated with the mean GP.
#' @param kern_i A kernel function, associated with the individual GPs.
#' @param hp_k A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_0}.
#' @param hp_i A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#' @param m_k prior means of the mu_k processes.
#' @param old_tau_i_k values au tau_i_k from previous iterations. List(list(tau)_i)_k
#'
#' @return A named list, containing the elements \code{mean}, a tibble
#' containing the Input and associated Output of the hyper-posterior mean
#' parameter, and \code{cov}, the hyper-posterior covariance matrix.
#'
#' @return A named list, containing the elements \code{mean}, a tibble
#' containing the Input and associated Output of the hyper-posterior mean
#' parameter, and \code{cov}, the hyper-posterior covariance matrix.
#' @export
#'
#' @examples
e_step_VEM = function(db, m_k, kern_0, kern_i, hp_k, hp_i, old_tau_i_k, pen_diag = NULL)
{
  pi_k = hp_k$pi
  all_t = unique(db$Input) %>% sort()
  t_clust = tibble::tibble('ID' = rep(names(m_k), each = length(all_t)),
                           'Timestamp' = rep(all_t, length(m_k)))

  inv_k = kern_to_inv(t_clust, kern_0, hp_k)
  list_inv_i = kern_to_inv(db, kern_i, hp_i)
  value_i = base::split(db$Output, list(db$ID))

  ## Update each mu_k parameters
  floop = function(k)
  {
    new_inv = inv_k[[k]]; tau_i_k = old_tau_i_k[[k]]
    for(x in list_inv_i %>% names())
    {
      inv_i = list_inv_i[[x]]
      common_times = intersect(row.names(inv_i), row.names(new_inv))
      new_inv[common_times, common_times] = new_inv[common_times, common_times] +
        tau_i_k[[x]] * inv_i[common_times, common_times]
    }


    tryCatch(solve(new_inv), error = function(e){MASS::ginv(new_inv)}) %>%
      return()
  }
  cov_k = sapply(names(m_k), floop, simplify = FALSE, USE.NAMES = TRUE)

  floop2 = function(k)
  {
    prior_mean = m_k[[k]]; prior_inv = inv_k[[k]]; tau_i_k = old_tau_i_k[[k]]

    if(length(prior_mean) == 1){prior_mean = rep(prior_mean, ncol(prior_inv))}
    weighted_mean = prior_inv %*% prior_mean
    #row.names(weithed_mean) = row.names(prior_inv)

    for(i in list_inv_i %>% names())
    {
      weighted_i = tau_i_k[[i]] * list_inv_i[[i]] %*% value_i[[i]]
      #row.names(weithed_i) = row.names(list_inv_i[[j]])

      common_times = intersect(row.names(weighted_i), row.names(weighted_mean))
      weighted_mean[common_times,] = weighted_mean[common_times,] + weighted_i[common_times,]
    }


    new_mean = cov_k[[k]] %*% weighted_mean %>% as.vector()
    tibble::tibble('Timestamp' = all_t, 'Output' = new_mean) %>% return()
  }
  mean_k = sapply(names(m_k), floop2, simplify = FALSE, USE.NAMES = TRUE)

  ## Update tau_i_k
  c_i = 0
  c_k = 0
  mat_logL = matrix(NA, nrow = length(names(m_k)), ncol = length(unique(db$ID)) )

  for(i in unique(db$ID))
  { c_i = c_i + 1
  t_i = db %>% dplyr::filter(.data$ID == i) %>% dplyr::pull(.data$Input)
  for(k in names(m_k))
  { c_k = c_k + 1

  mat_logL[c_k,c_i] = - logL_GP_mod(hp_i[[i]], db %>%
                                      dplyr::filter(.data$ID == i), mean_k[[k]] %>%
                                      dplyr::filter(.data$Input %in% t_i) %>%
                                      dplyr::pull(.data$Output), kern_i,
                                    cov_k[[k]][paste0('X',t_i),paste0('X',t_i)])
  if(is.na(mat_logL[c_k,c_i])){print(i)}
  }
  c_k = 0
  }

  ## We need to use the 'log-sum-exp' trick: exp(x - max(x)) / sum exp(x - max(x)) to remain numerically stable
  mat_L = mat_logL %>% apply(2,function(x) exp(x - max(x)))

  tau_i_k = (pi_k * mat_L) %>% apply(2,function(x) x / sum(x)) %>%
    `rownames<-`(names(m_k)) %>%
    `colnames<-`(unique(db$ID)) %>%
    apply(1, as.list)


  list('mean' = mean_k, 'cov' = cov_k, 'tau_i_k' = tau_i_k) %>% return()

}


#' Variation of the M-Step of the EM algorithm
#'
#' Maximization step of the variation of the EM algorithm to compute hyper-parameters of all the
#' kernels involved in Magma.
#'
#' @param db A tibble or data frame. Columns required: ID, Input, Output.
#'    Additional columns for covariates can be specified.
#' @param list_mu_param List of parameters of the K mean GPs. Format list('mean', 'cov', 'tau_i_k')
#' @param kern_0 kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0)
#' @param kern_i kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i)
#' @param m_k prior value of the mean parameter of the mean GPs (mu_k). Length = 1 or nrow(unique(db$Input))
#' @param common_hp_k boolean indicating whether hp are common among mean GPs (for each mu_k)
#' @param common_hp_i boolean indicating whether hp are common among individual GPs (for each y_i)
#' @param old_hp_i A tibble or data frame, containing the hyper-parameters
#'    from the previous M-step (or initialisation) associated with the
#'    individual GPs.
#' @param old_hp_k A tibble or data frame, containing the hyper-parameters
#'    from the previous M-step (or initialisation) associated with the variations
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#' numerical issues when inverting nearly singular matrices.
#'
#' @return Set of optimised hyper parameters for the different kernels of the model, and the pi_k
#' @export
#'
#' @examples
m_step_VEM = function(db, old_hp_k, old_hp_i, list_mu_param, kern_0, kern_i, m_k, common_hp_k, common_hp_i, pen_diag)
{
  list_ID_k = names(m_k)
  list_ID_i = unique(db$ID)

  list_hp_i <- old_hp_i %>%
    dplyr::select(-.data$ID) %>%
    names()

  t1 = Sys.time()
  if(common_hp_i)
  {
    param = optimr::opm(old_hp_i %>% dplyr::select(-.data$ID) %>% dplyr::slice(1),
                logL_clust_multi_GP_common_hp_i,
                gr = gr_clust_multi_GP_common_hp_i,
                db = db,
                mu_k_param = list_mu_param,
                kern = kern_i,
                method = "L-BFGS-B",
                control = list(kkt = F)) %>%
      dplyr::select(list_hp_i) %>%
      tibble::as_tibble() %>%
      tidyr::uncount(weights = length(list_ID_i)) %>%
      dplyr::mutate('ID' = list_ID_i, .before = 1)

    new_theta_i = param %>% list() %>% rep(length(list_ID_i))  %>%
      stats::setNames(nm = list_ID_i)
  }
  else
  {
    funloop2 = function(i)
    {
      t_i = db %>% dplyr::filter(.data$ID == i) %>% dplyr::pull(.data$Input)
      optimr::opm(old_hp_i %>% dplyr::select(-.data$ID) %>% dplyr::slice(i),
          logL_clust_multi_GP,
          gr = gr_clust_multi_GP,
          db = db %>% dplyr::filter(.data$ID == i),
          mu_k_param = list_mu_param,
          kern = kern_i,
          method = "L-BFGS-B",
          control = list(kkt = F)) %>%
        dplyr::select(.data$list_hp_i) %>%
        tibble::as_tibble() %>%
        tidyr::uncount(weights = length(list_ID_i)) %>%
        dplyr::mutate('ID' = list_ID_i, .before = 1) %>% return()
    }
    new_theta_i = sapply(list_ID_i, funloop2, simplify = FALSE, USE.NAMES = TRUE)
  }
  t2 = Sys.time()

  if(common_hp_k)
  {
    param = c(optimr::opm(old_hp_k %>% dplyr::select(-.data$ID) %>% dplyr::slice(1),
                  logL_GP_mod_common_hp_k,
                  gr = gr_GP_mod_common_hp_k,
                  db = list_mu_param$mean,
                  mean = m_k,
                  kern = kern_0,
                  new_cov = list_mu_param$cov,
                  pen_diag = pen_diag,
                  method = "L-BFGS-B",
                  control = list(kkt = F))[1,1:2] %>%
                unlist(), pen_diag)
    new_theta_k = param %>% list() %>% rep(length(list_ID_k))  %>%
      stats::setNames(nm = list_ID_k)
  }
  else
  {
    funloop = function(k)
    {
      c(optimr::opm(old_hp_k %>% dplyr::select(-.data$ID) %>% dplyr::slice(k),
            logL_GP_mod,
            gr = gr_GP_mod,
            db = list_mu_param$mean[[k]],
            mean = m_k[[k]],
            kern = kern_0,
            new_cov = list_mu_param$cov[[k]],
            pen_diag = pen_diag,
            method = "L-BFGS-B",
            control = list(kkt = FALSE))[1,1:2] %>%
          unlist(), pen_diag) %>% return()
    }
    new_theta_k = sapply(list_ID_k, funloop, simplify = FALSE, USE.NAMES = TRUE)
  }

  pi_k = sapply( list_mu_param$tau_i_k, function(x) x %>% unlist() %>% mean() )

  t3 = Sys.time()
  print(c('mu_k:',t3 - t2, 'indiv:', t2 - t1))

  list(
    "hp_k" = new_theta_k,
    "hp_i" = new_theta_i,
    'pi_k' = pi_k
  ) %>%
    return()
}
