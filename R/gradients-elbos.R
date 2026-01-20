#' Gradient of the elbo for a mixture of GPs
#'
#' @param hp A tibble, data frame or named vector containing hyper-parameters.
#' @param db A tibble containing the values we want to compute the elbo on.
#'    Required columns: Task_ID, Input, Output, Output_ID. Additional covariate
#'    columns are allowed.
#' @param hyperpost List of parameters for the K mean Gaussian processes.
#' @param kern A kernel function.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#' @param hp_col_names A character vector with the names of the hyper-parameters.
#' @param output_ids A character vector with the unique IDs of the outputs.
#'
#' @return The gradient of the penalised Gaussian elbo for a mixture of GPs
#'
#' @keywords internal
#'
#' @examples
#' TRUE
gr_clust_multi_GP <- function(hp,
                              db,
                              hyperpost,
                              kern,
                              pen_diag,
                              hp_col_names,
                              output_ids) {
  list_hp <- names(hp)

  names_k <- hyperpost$mean %>% names()
  t_t <- db$Reference
  y_t <- db$Output
  output_ids <- db$Output_ID %>% unique()

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
    inputs <- db %>% dplyr::select(-c(Output, Task_ID))
  } else{
    inputs <- db %>% dplyr::select(-c(Output))
  }

  t <- unique(db$Task_ID)

  corr1 <- 0
  corr2 <- 0

  for (k in (names_k))
  {
    tau_t_k <- hyperpost$mixture %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::pull(k)

    mean_mu_k <- hyperpost$mean[[k]] %>%
      dplyr::filter(Reference %in% t_t) %>%
      dplyr::arrange(match(Reference, t_t)) %>%
      dplyr::pull(Output)

    corr1 <- corr1 + tau_t_k * mean_mu_k
    corr2 <- corr2 + tau_t_k *
      (mean_mu_k %*% t(mean_mu_k) +
         hyperpost$cov[[k]][as.character(t_t), as.character(t_t)])
  }

  if(length(output_ids) == 1){
    inputs <- inputs %>% dplyr::select(-Output_ID)
  }

  inv <- kern_to_inv(inputs, kern, hp_tibble, pen_diag)
  prod_inv <- inv %*% y_t

  common_term <- (prod_inv - 2 * inv %*% corr1) %*% t(prod_inv) +
    inv %*% (corr2 %*% inv - diag(1, length(t_t)))

  ## Loop over the derivatives of hyper-parameters for computing the gradient
  floop <- function(deriv) {
    (-1 / 2 * (common_term %*% kern_to_cov(inputs, kern, hp_tibble, deriv))) %>%
      diag() %>%
      sum() %>%
      return()
  }
  sapply(list_hp, floop) %>%
    return()
}



#' Gradient of the penalised elbo for multiple task GPs with shared HPs
#'
#' @param hp A tibble, data frame or name vector of hyper-parameters.
#' @param db A tibble containing values we want to compute elbo on.
#'    Required columns: Task_ID, Input, Output, Output_ID. Additional covariate
#'    columns are allowed.
#' @param kern A kernel function used to compute the covariance matrix at
#' corresponding timestamps.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#' @param hyperpost List of parameters for the K mean Gaussian processes.
#' @param hp_col_names A character vector with the names of the hyper-parameters.
#' @param output_ids A character vector with the unique IDs of the outputs.
#'
#' @return The gradient of the penalised Gaussian elbo for
#'    the sum of the T task GPs with shared HPs.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
gr_clust_multi_GP_shared_hp_tasks <- function(hp,
                                          db,
                                          hyperpost,
                                          kern,
                                          pen_diag = NULL,
                                          hp_col_names,
                                          output_ids) {
  list_hp <- names(hp)
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

  ## Loop over tasks to compute the sum of log-Likelihoods
  funloop <- function(t) {
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
      ## Extract the covariance values associated with the i-th specific inputs
      post_cov_t <- hyperpost$cov[[k]][
        as.character(input_t),
        as.character(input_t)
      ]

      tau_t_k <- hyperpost$mixture %>%
        dplyr::filter(Task_ID == t) %>%
        dplyr::pull(k)

      mean_mu_k <- hyperpost$mean[[k]] %>%
        dplyr::filter(Reference %in% input_t) %>%
        dplyr::arrange(match(Reference, input_t)) %>% # CORRECTION
        dplyr::pull(Output)
      names(mean_mu_k) <- (hyperpost$mean[[k]] %>%
                             dplyr::filter(Reference %in% input_t) %>%
                             dplyr::arrange(match(Reference, input_t)))$Reference

      corr1 <- corr1 + tau_t_k * mean_mu_k
      corr2 <- corr2 + tau_t_k *
        (mean_mu_k %*% t(mean_mu_k) + post_cov_t)
    }

    if(length(output_ids) == 1){
      inputs_t <- inputs_t %>% dplyr::select(-Output_ID)
    }

    inv <- kern_to_inv(inputs_t, kern, hp_tibble, pen_diag)

    prod_inv <- inv %*% output_t
    common_term <- (prod_inv - 2 * inv %*% corr1) %*% t(prod_inv) +
      inv %*% (corr2 %*% inv - diag(1, length(input_t)))

    floop <- function(deriv) {
      (-1 / 2 * (common_term %*% kern_to_cov(inputs_t, kern, hp_tibble, deriv))) %>%
        diag() %>%
        sum() %>%
        return()
    }

    sapply(list_hp, floop) %>%
      return()
  }

  sapply(unique(db$Task_ID), funloop) %>%
    rowSums() %>%
    return()
}
