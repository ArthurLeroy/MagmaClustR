#' Evidence Lower Bound for a mixture of GPs
#'
#' @param hp A tibble, data frame or named vector containing hyper-parameters.
#' @param db A tibble containing the values we want to compute the elbo on.
#'    Required columns: Task_ID, Input, Output, Output_ID. Additional covariate
#'    columns are allowed.
#' @param hyperpost List of parameters for the K mean GPs.
#' @param kern A kernel function used to compute the covariance matrix at
#' corresponding timestamps.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#' @param hp_col_names A character vector with the names of the hyper-parameters
#'   (e.g., c("l_t", "S_t")).
#' @param output_ids A character vector with the unique IDs of the outputs for
#'   the current task.
#'
#' @return The value of the penalised Gaussian elbo for a mixture of GPs
#'
#' @keywords internal
#'
#' @examples
#' TRUE
elbo_clust_multi_GP <- function(hp,
                                db,
                                hyperpost,
                                kern,
                                pen_diag,
                                hp_col_names,
                                output_ids) {

  names_k <- hyperpost$mean %>% names()
  t_t <- db$Reference
  y_t <- db$Output
  t <- unique(db$Task_ID)
  output_ids <- unique(db$Output_ID)

  if(!(hp %>% tibble::is_tibble()) && length(output_ids) > 1){
    # Reconstruct the structured HP tibble from the flat vector
    hp_tibble <- reconstruct_hp(
      par_vector = hp,
      hp_names = hp_col_names,
      output_ids = output_ids
    )
  } else if (!(hp %>% tibble::is_tibble()) && length(output_ids) == 1){
    hp_tibble <- hp %>%
      t() %>%
      tibble::as_tibble() %>%
      stats::setNames(hp_col_names)

  } else {
    hp_tibble <- hp
  }

  if("Task_ID" %in% names(db)){
    inputs <- db %>% dplyr::select(-c(Output, -Task_ID))
  } else{
    inputs <- db %>% dplyr::select(-c(Output))
  }

  if(length(output_ids) == 1){
    inputs <- inputs %>% dplyr::select(-Output_ID)
  }

  inv <- kern_to_inv(inputs, kern, hp_tibble, pen_diag)

  ## classic Gaussian centred log likelihood
  LL_norm <- -dmnorm(y_t, rep(0, length(y_t)), inv, log = T)

  corr1 <- 0
  corr2 <- 0

  for (k in (names_k))
  {
    tau_t_k <- tau_t_k <- hyperpost$mixture %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::pull(k)
    mean_mu_k <- hyperpost$mean[[k]] %>%
      dplyr::filter(Reference %in% t_t) %>%
      dplyr::pull(Output)
    corr1 <- corr1 + tau_t_k * mean_mu_k
    corr2 <- corr2 + tau_t_k *
      (mean_mu_k %*% t(mean_mu_k) +
         hyperpost$cov[[k]][as.character(t_t), as.character(t_t)])
  }

  (LL_norm - y_t %*% inv %*% corr1 + 0.5 * sum(inv * corr2)) %>% return()
}

#' #' Penalised elbo for multiple mean GPs with shared HPs
#' #'
#' #' @param hp A tibble, data frame or named vector containing hyper-parameters.
#' #' @param db  A tibble containing values we want to compute elbo on.
#' #'    Required columns: Task_ID, Input, Output, Output_ID. Additional covariate
#' #'    columns are allowed.
#' #' @param mean A list of the K mean GPs at union of observed timestamps.
#' #' @param kern A kernel function used to compute the covariance matrix at
#' #' corresponding timestamps.
#' #' @param post_cov A List of the K posterior covariance of the mean GP (mu_k).
#' #' Used to compute correction term (cor_term).
#' #' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#' #'    numerical issues when inverting, in cases of nearly singular matrices.
#' #' @param hp_col_names A character vector with the names of the hyper-parameters
#' #'   (e.g., c("l_t", "S_t")).
#' #' @param output_ids A character vector with the unique IDs of the outputs for
#' #'   the current task.
#' #'
#' #'
#' #' @return The value of the penalised Gaussian elbo for
#' #'    the sum of the k mean GPs with common HPs.
#' #'
#' #' @keywords internal
#' #'
#' #' @examples
#' #' TRUE
#' elbo_GP_mod_shared_hp_clust <- function(hp,
#'                                     db,
#'                                     mean,
#'                                     kern,
#'                                     post_cov,
#'                                     pen_diag,
#'                                     hp_col_names,
#'                                     output_ids) {
#'
#'   # browser() # -> Fonctionne !
#'   list_ID_k <- names(db)
#'   if(!(hp %>% tibble::is_tibble()) && length(output_ids) > 1){
#'     # Reconstruct the structured HP tibble from the flat vector
#'     hp_tibble <- reconstruct_hp(
#'       par_vector = hp,
#'       hp_names = hp_col_names,
#'       output_ids = output_ids
#'     )
#'   } else if (!(hp %>% tibble::is_tibble()) && length(output_ids) == 1){
#'     hp_tibble <- hp %>%
#'       t() %>%
#'       tibble::as_tibble() %>%
#'       stats::setNames(hp_col_names)
#'
#'   } else {
#'     hp_tibble <- hp
#'   }
#'
#'   if("Task_ID" %in% names(db)){
#'     inputs <- db[[1]] %>% dplyr::select(-c(Output, Task_ID, Output_ID))
#'   } else{
#'     inputs <- db[[1]] %>% dplyr::select(-c(Output))
#'   }
#'
#'   if(length(output_ids) == 1){
#'     inputs <- inputs %>% dplyr::select(-Output_ID)
#'   }
#'
#'   inv <- kern_to_inv(inputs, kern, hp_tibble, pen_diag)
#'
#'   LL_norm <- 0
#'   cor_term <- 0
#'
#'   for (k in list_ID_k)
#'   {
#'     y_k <- db[[k]] %>% dplyr::pull(Output)
#'
#'     LL_norm <- LL_norm - dmnorm(y_k, mean[[k]], inv, log = T)
#'     cor_term <- cor_term + 0.5 * (inv * post_cov[[k]]) %>% sum()
#'   }
#'   return(LL_norm + cor_term)
#' }

#' Penalised elbo for multiple individual GPs with shared HPs
#'
#' @param hp A tibble, data frame or named vector containing hyper-parameters.
#' @param db A tibble containing values we want to compute elbo on.
#'    Required columns: Task_ID, Input, Output, Output_ID. Additional covariate
#'    columns are allowed.
#' @param hyperpost List of parameters for the K mean Gaussian processes.
#' @param kern A kernel function used to compute the covariance matrix at
#'    corresponding timestamps.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#' @param hp_col_names A character vector with the names of the hyper-parameters
#'   (e.g., c("l_t", "S_t")).
#' @param output_ids A character vector with the unique IDs of the outputs for
#'   the current task.
#'
#' @return The value of the penalised Gaussian elbo for
#'    the sum of the M individual GPs with common HPs.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
elbo_clust_multi_GP_shared_hp_tasks <- function(hp,
                                            db,
                                            hyperpost,
                                            kern,
                                            pen_diag,
                                            hp_col_names,
                                            output_ids) {

  # browser() # -> Fonctionne !
  names_k <- hyperpost$mean %>% names()

  if(!(hp %>% tibble::is_tibble()) && length(output_ids) > 1){
    # Reconstruct the structured HP tibble from the flat vector
    hp_tibble <- reconstruct_hp(
      par_vector = hp,
      hp_names = hp_col_names,
      output_ids = output_ids
    )
  } else if (!(hp %>% tibble::is_tibble()) && length(output_ids) == 1){
    hp_tibble <- hp %>%
      t() %>%
      tibble::as_tibble() %>%
      stats::setNames(hp_col_names)

  } else {
    hp_tibble <- hp
  }

  sum_t <- 0
  for (t in unique(db$Task_ID))
  {
    ## Extract the t-th specific reference Input
    input_t <- db %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::pull(Reference)
    ## Extract the t-th specific inputs (reference + covariates)
    inputs_t <- db %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::select(-c(Task_ID, Output))
    ## Extract the t-th specific Inputs and Output
    output_t <- db %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::pull(Output)

    corr1 <- 0
    corr2 <- 0

    for (k in (names_k))
    {
      ## Extract the covariance values associated with the t-th specific inputs
      post_cov_t <- hyperpost$cov[[k]][
        as.character(input_t),
        as.character(input_t)
      ]

      tau_t_k <- hyperpost$mixture %>%
        dplyr::filter(Task_ID == t) %>%
        dplyr::pull(k)
      mean_mu_k <- hyperpost$mean[[k]] %>%
        dplyr::filter(Reference %in% input_t) %>%
        dplyr::pull(Output)
      corr1 <- corr1 + tau_t_k * mean_mu_k
      corr2 <- corr2 + tau_t_k *
        (mean_mu_k %*% t(mean_mu_k) + post_cov_t)
    }

    if(length(output_ids) == 1){
      inputs_t <- inputs %>% dplyr::select(-Output_ID)
    }
    inv <- kern_to_inv(inputs_t, kern, hp_tibble, pen_diag)

    ## Classic Gaussian centred log-likelihood
    LL_norm <- -dmnorm(output_t, rep(0, length(output_t)), inv, log = T)

    sum_t <- sum_t +
      LL_norm - output_t %*% inv %*% corr1 + 0.5 * sum(inv * corr2)
  }
  return(sum_t)
}

#' Evidence Lower Bound maximised in MagmaClust
#'
#' @param hp_k A tibble, data frame or named vector of hyper-parameters
#'    for each cluster.
#' @param hp_t A tibble, data frame or named vector of hyper-parameters
#'    for each task.
#' @param db A tibble containing values we want to compute elbo on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param kern_t Kernel used to compute the covariance matrix of task GPs
#'    at corresponding inputs.
#' @param kern_k Kernel used to compute the covariance matrix of the mean GPs
#'    at corresponding inputs.
#' @param hyperpost A list of parameters for the variational distributions
#'    of the K mean GPs.
#' @param m_k Prior value of the mean parameter of the mean GPs (mu_k).
#'    Length = 1 or nrow(db).
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#' @param hp_col_names A character vector with the names of the hyper-parameters
#'   (e.g., c("l_t", "S_t")).
#' @param output_ids A character vector with the unique IDs of the outputs for
#'   the current task.
#'
#'
#' @return Value of the elbo that is maximised during the VEM algorithm used for
#'    training in MagmaClust.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
elbo_monitoring_VEM <- function(hp_k,
                                hp_t,
                                db,
                                kern_t,
                                kern_k,
                                hyperpost,
                                m_k,
                                pen_diag,
                                hp_col_names,
                                output_ids) {
  floop <- function(k) {
    logL_GP_mod(
      hp_k[hp_k$Cluster_ID == k, ],
      db = hyperpost$mean[[k]],
      mean = m_k[[k]],
      kern_k,
      hyperpost$cov[[k]],
      pen_diag
    ) %>%
      return()
  }
  sum_ll_k <- sapply(names(m_k), floop) %>% sum()

  floop2 <- function(t) {
    t_t <- db %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::pull(Reference)

    elbo_clust_multi_GP(
      hp_t[hp_t$Task_ID == t, ],
      db %>% dplyr::filter(Task_ID == t),
      hyperpost,
      kern_t,
      pen_diag
    ) %>%
      return()
  }
  sum_ll_t <- sapply(unique(db$Task_ID), floop2) %>% sum()

  floop3 <- function(k) {
    sum_tau <- 0
    det <- 0
    ## Extract the proportion in the k-th cluster
    pi_k <- hp_k %>%
      dplyr::filter(Cluster_ID == k) %>%
      dplyr::pull(prop_mixture)

    for (t in unique(db$Task_ID)) {
      ## Extract the probability of the t-th indiv to be in the k-th cluster
      tau_t_k <- hyperpost$mixture %>%
        dplyr::filter(Task_ID == t) %>%
        dplyr::pull(k)
      ## To avoid numerical issues if evaluating log(0/0)
      log_frac <- ifelse((tau_t_k == 0 | (pi_k == 0)), 0, log(pi_k / tau_t_k))

      sum_tau <- sum_tau + tau_t_k * log_frac
    }
    ## Compute the sum of log-determinant terms using Cholesky decomposition
    ## log(det(A)) = 2*sum(log(diag(chol(A))))
    det <- det + hyperpost$cov[[k]] %>%
      chol() %>%
      diag() %>%
      log() %>%
      sum()

    return(sum_tau + det)
  }

  sum_corr_k <- sapply(names(m_k), floop3) %>% sum()

  return(-sum_ll_k - sum_ll_t + sum_corr_k)
}
