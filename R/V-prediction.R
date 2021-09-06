#' Posterior mu of the cluster
#'
#' @param db data A tibble or data frame. Required columns: \code{ID}, \code{Input}
#'    , \code{Output}. Additional columns for covariates can be specified.
#'    The \code{ID} column contains the unique names/codes used to identify each
#'    individual/task (or batch of data).
#'    The \code{Input} column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    \code{Output} column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference \code{Input}.
#' @param timestamps timestamps on which we want a prediction
#' @param kern_0 kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0).
#' @param kern_i kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i).
#' @param m_k prior value of the mean parameter of the mean GPs (mu_k). Length = 1 or nrow(db)
#' @param list_hp list of your hyperparameters
#'
#' @return Pamameters of the mean GP at timestamps chosen.
#' @export
#'
#' @examples
posterior_mu_k = function(db, timestamps, m_k, kern_0, kern_i, list_hp)
{
  #browser()
  hp_i = list_hp$hp_i
  hp_k = list_hp$hp_k
  #for(k in names(hp_k)){hp_k[[k]][3] = 0.1}
  tau_i_k = list_hp$param$tau_i_k
  #t_clust = tibble::tibble('ID' = rep(names(hp_k), each = length(timestamps)) , 'Timestamp' = rep(timestamps, length(hp_k)))
  t_clust = tibble::tibble('ID' = rep(hp_k$ID, each = length(timestamps)) , 'Input' = rep(timestamps, length(hp_k$ID)))
  inv_k = list_kern_to_inv(t_clust, kern_0, hp_k)
  list_inv_i = list_kern_to_inv(db, kern_i, hp_i)
  value_i = base::split(db$Output, list(db$ID))

  ## Update each mu_k parameters
  floop = function(k)
  {
    new_inv = inv_k[[k]]
    for(x in list_inv_i %>% names())
    {
      inv_i = list_inv_i[[x]]
      common_times = intersect(row.names(inv_i), row.names(new_inv))
      new_inv[common_times, common_times] = new_inv[common_times, common_times] +
        tau_i_k[[k]][[x]] * inv_i[common_times, common_times]
    }
    tryCatch(solve(new_inv), error = function(e){
      s_inv<- MASS::ginv(new_inv)
      colnames(s_inv) <- colnames(new_inv)
      rownames(s_inv) <- rownames(new_inv)
      s_inv
      }) %>%
      return()
  }
  cov_k = sapply(hp_k$ID, floop, simplify = FALSE, USE.NAMES = TRUE)

  floop2 = function(k)
  {
    prior_mean <- m_k[[k]]
    if(length(prior_mean) == 1){prior_mean = rep(prior_mean, ncol(inv_k[[k]]))}
    weighted_mean = inv_k[[k]] %*% prior_mean
    #row.names(weithed_mean) = row.names(inv_k[[k]])

    for(i in list_inv_i %>% names())
    {
      weighted_i = tau_i_k[[k]][[i]] * list_inv_i[[i]] %*% value_i[[i]]
      #row.names(weithed_i) = row.names(list_inv_i[[j]])

      common_times = intersect(row.names(weighted_i), row.names(weighted_mean))
      weighted_mean[common_times,] = weighted_mean[common_times,] + weighted_i[common_times,]
    }

    new_mean = cov_k[[k]] %*% weighted_mean %>% as.vector()
    #tibble::tibble('Timestamp' = timestamps, 'Output' = new_mean) %>% return()
    tibble::tibble('Input' = timestamps, 'Output' = new_mean) %>% return()
  }
  mean_k = sapply(hp_k$ID, floop2, simplify = FALSE, USE.NAMES = TRUE)

  #names(mean_mu) = paste0('X', t_mu)
  list('mean' = mean_k, 'cov' = cov_k, 'tau_i_k' = tau_i_k) %>% return()
}

#' Prediction Gaussian Process on the clustering
#'
#' @param db A tibble or data frame. Required columns: 'Input',
#'    'Output'. Additional columns for covariates can be specified.
#'    The 'Input' column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    'Output' column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference 'Input'.
#' @param timestamps timestamps on which we want a prediction
#' @param hp A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern}. The columns/elements should be named
#'    according to the hyper-parameters that are used in \code{kern}.
#'    List containing hyper-parameters and tau_k for the new individual
#' @param kern A kernel function, defining the covariance structure of the GP.
#'    Several popular kernels
#'    (see \href{https://www.cs.toronto.edu/~duvenaud/cookbook/}{The Kernel
#'    Cookbook}) are already implemented and can be selected within the
#'    following list:
#'    - "SE": (default value) the Squared Exponential Kernel (also called
#'        Radial Basis Function or Gaussian kernel),
#'    - "LIN": the Linear kernel,
#'    - "PERIO": the Periodic kernel,
#'    - "RQ": the Rational Quadratic kernel.
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#' @param list_mu List containing mean and cov of the K mean GPs.
#' - mean   : value of mean GP at timestamps (obs + pred) (matrix dim: timestamps x 1, with Input rownames)
#' - cov_mu : covariance of mean GP at timestamps (obs + pred) (square matrix, with Input row/colnames)
#'
#' @return pamameters of the gaussian density predicted at timestamps
#' @export
#'
#' @examples
#'
#' k = seq_len(2)
#' m_k <- c("K1" = 0, "K2" = 0)
#'
#' db <- simu_db(N = 2, common_input = FALSE)
#' hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
#' hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))

#' ini_hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))
#' old_tau_i_k = MagmaClustR:::ini_tau_i_k(db = db, k = length(k), nstart = 50)
#'
#' training_test = train_magma_VEM(db, m_k, hp_k, ini_hp_i, "SE", "SE", old_tau_i_k, FALSE, FALSE, 0.1)
#'
#' timestamps = seq(0.01, 10, 0.01)
#' mu_k <- posterior_mu_k(db, timestamps, m_k, "SE", "SE", training_test)
#'
#'
#' list_hp <- train_magma_VEM(db, m_k, hp_k, hp_i, "SE", "SE", old_tau_i_k, FALSE, FALSE, 0.1)
#'
#' new_indiv <- train_new_gp_EM(simu_db(M=1, covariate = FALSE), mu_k, ini_hp_i, "SE", hp_i = list_hp$hp_i)
#'
#' pred_gp_clust(simu_db(M=1, covariate = FALSE), timestamps, mu_k, kern = se_kernel, new_indiv)
pred_gp_clust = function(db, timestamps = NULL, list_mu, kern, hp)
{
  #browser()

  inputs = db %>% dplyr::pull(.data$Input)
  #input_n <- db['Input']
  yn = db %>% dplyr::pull(.data$Output)

  if(is.null(timestamps)){timestamps = seq(min(inputs), max(inputs), length.out = 500)}

  hp_new = hp$theta_new
  tau_k = hp$tau_k

  floop = function(k)
  {
    mean_mu_obs = list_mu$mean[[k]] %>% dplyr::filter(.data$Input %in% inputs) %>% dplyr::pull(.data$Output)
    mean_mu_pred = list_mu$mean[[k]] %>% dplyr::filter(.data$Input %in% timestamps) %>% dplyr::pull(.data$Output)
    cov_mu = list_mu$cov[[k]]

    cov_tn_tn = (kern_to_cov(inputs, kern, hp_new) + cov_mu[as.character(inputs), as.character(inputs)])
    inv_mat = tryCatch(solve(cov_tn_tn), error = function(e){MASS::ginv(cov_tn_tn)})
    cov_tn_t = kern(inputs, timestamps, hp_new) + cov_mu[as.character(inputs), as.character(timestamps)]
    cov_t_t = kern_to_cov(timestamps, kern, hp_new) + cov_mu[as.character(timestamps), as.character(timestamps)]

    tibble::tibble('Input' = timestamps, ##'Timestamp' = timestamps,
           'Mean' = (mean_mu_pred + t(cov_tn_t) %*% inv_mat %*% (yn - mean_mu_obs)) %>% as.vector(),
           'Var' =  (cov_t_t - t(cov_tn_t) %*% inv_mat %*% cov_tn_t) %>% diag %>% as.vector,
           'tau_k' = tau_k[[k]]) %>% return()
  }
  pred = sapply(names(list_mu$mean), floop, simplify = FALSE, USE.NAMES = TRUE)

  return(pred)
}


#' Prediction animate of the Gaussian Process on the clustering
#'
#' @param db A tibble or data frame. Required columns: 'Input',
#'    'Output'. Additional columns for covariates can be specified.
#'    The 'Input' column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    'Output' column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference 'Input'.
#' @param timestamps timestamps on which we want a prediction
#' @param hp A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern}. The columns/elements should be named
#'    according to the hyper-parameters that are used in \code{kern}.
#'    List containing hyper-parameters and tau_k for the new individual.
#' @param kern A kernel function, defining the covariance structure of the GP.
#'    Several popular kernels
#'    (see \href{https://www.cs.toronto.edu/~duvenaud/cookbook/}{The Kernel
#'    Cookbook}) are already implemented and can be selected within the
#'    following list:
#'    - "SE": (default value) the Squared Exponential Kernel (also called
#'        Radial Basis Function or Gaussian kernel),
#'    - "LIN": the Linear kernel,
#'    - "PERIO": the Periodic kernel,
#'    - "RQ": the Rational Quadratic kernel.
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#' @param list_mu List containing mean and cov of the K mean GPs.
#' - mean   : value of mean GP at timestamps (obs + pred) (matrix dim: timestamps x 1, with Input rownames)
#' - cov_mu : covariance of mean GP at timestamps (obs + pred) (square matrix, with Input row/colnames)
#'
#' @return tibble of classic GP predictions but with an inscreasing number of data points considered as 'observed'
#' @export
#'
#' @examples
#'
#' k = seq_len(2)
#' m_k <- c("K1" = 0, "K2" = 0)
#'
#' db <- simu_db(N = 2, common_input = FALSE)
#' hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
#' hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))

#' ini_hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))
#' old_tau_i_k = MagmaClustR:::ini_tau_i_k(db = db, k = length(k), nstart = 50)
#'
#' training_test = train_magma_VEM(db, m_k, hp_k, ini_hp_i, "SE", "SE", old_tau_i_k, FALSE, FALSE, 0.1)
#'
#' timestamps = seq(0.01, 10, 0.01)
#' mu_k <- posterior_mu_k(db, timestamps, m_k, "SE", "SE", training_test)
#'
#'
#' list_hp <- train_magma_VEM(db, m_k, hp_k, hp_i, "SE", "SE", old_tau_i_k, FALSE, FALSE, 0.1)
#' new_indiv <- train_new_gp_EM(db, mu_k, ini_hp_i, "SE", hp_i = list_hp$hp_i)
#'
#' pred_gp_clust_animate(db, timestamps, mu_k, kern = se_kernel, new_indiv)
pred_gp_clust_animate = function(db, timestamps = NULL, list_mu, kern, hp)
{
  browser()
  db <- db %>% dplyr::arrange(.data$Input)
  all_pred = tibble::tibble()

  if(is.null(timestamps)){timestamps = seq(min(db$Input), max(db$Input), length.out = 500)}

  for(j in 1:nrow(db))
  {
    pred_j = pred_gp_clust(db[1:j,], timestamps, list_mu, kern, hp)

    for(k in list_mu$mean %>% names)
    {
      pred_j[[k]] = pred_j[[k]] %>% dplyr::mutate(Nb_data = j, Cluster = k)
    }
    pred_j_all = pred_j %>% dplyr::bind_rows

    all_pred = all_pred %>% rbind(pred_j_all)
  }
  return(all_pred)
}

#' Prediction of the maximum of cluster
#'
#' @param tau_i_k tau_i_k
#'
#' @return Prediction of the maximum of cluster
#' @export
#'
#' @examples
pred_max_cluster = function(tau_i_k)
{
  tau_i_k %>% tibble::as_tibble %>% tidyr::unnest %>% apply(1, which.max) %>% return()
}
