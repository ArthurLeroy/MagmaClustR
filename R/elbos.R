#' Evidence Lower Bound for a mixture of GPs
#'
#' @param hp A tibble, data frame or named vector containing hyper-parameters.
#' @param db A tibble containing the values we want to compute the elbo on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param hyperpost List of parameters for the K mean GPs.
#' @param kern A kernel function used to compute the covariance matrix at
#' corresponding timestamps.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return The value of the penalised Gaussian elbo for a mixture of GPs
#'
#' @keywords internal
#'
#' @examples
#' TRUE
elbo_clust_multi_GP <- function(hp, db, hyperpost, kern, pen_diag) {
  names_k <- hyperpost$mean %>% names()
  t_i <- db$Input
  y_i <- db$Output
  i <- unique(db$ID)

  inv <- kern_to_inv(t_i, kern, hp, pen_diag)

  ## classic Gaussian centred log likelihood
  LL_norm <- -dmnorm(y_i, rep(0, length(y_i)), inv, log = T)

  corr1 <- 0
  corr2 <- 0

  for (k in (names_k))
  {
    tau_i_k <- tau_i_k <- hyperpost$mixture %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(k)
    mean_mu_k <- hyperpost$mean[[k]] %>%
      dplyr::filter(.data$Input %in% t_i) %>%
      dplyr::pull(.data$Output)
    corr1 <- corr1 + tau_i_k * mean_mu_k
    corr2 <- corr2 + tau_i_k *
      (mean_mu_k %*% t(mean_mu_k) +
        hyperpost$cov[[k]][as.character(t_i), as.character(t_i)])
  }

  (LL_norm - y_i %*% inv %*% corr1 + 0.5 * sum(inv * corr2)) %>% return()
}

#' Penalised elbo for multiple mean GPs with common HPs
#'
#' @param hp A tibble, data frame or named vector containing hyper-parameters.
#' @param db  A tibble containing values we want to compute elbo on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mean A list of the K mean GPs at union of observed timestamps.
#' @param kern A kernel function used to compute the covariance matrix at
#' corresponding timestamps.
#' @param post_cov A List of the K posterior covariance of the mean GP (mu_k).
#' Used to compute correction term (cor_term).
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#'
#' @return The value of the penalised Gaussian elbo for
#'    the sum of the k mean GPs with common HPs.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
elbo_GP_mod_common_hp_k <- function(
  hp,
  db,
  mean,
  kern,
  post_cov,
  pen_diag
  ) {

  list_ID_k <- names(db)
  # t_k = db[[1]] %>% dplyr::pull(.data$Input)
  t_k <- db[[1]] %>%
    dplyr::pull(.data$Input)

  inv <- kern_to_inv(t_k, kern, hp, pen_diag)

  LL_norm <- 0
  cor_term <- 0

  for (k in list_ID_k)
  {
    y_k <- db[[k]] %>% dplyr::pull(.data$Output)

    LL_norm <- LL_norm - dmnorm(y_k, mean[[k]], inv, log = T)
    cor_term <- cor_term + 0.5 * (inv * post_cov[[k]]) %>% sum()
  }
  return(LL_norm + cor_term)
}

#' Penalised elbo for multiple individual GPs with common HPs
#'
#' @param hp A tibble, data frame or named vector containing hyper-parameters.
#' @param db A tibble containing values we want to compute elbo on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param hyperpost List of parameters for the K mean Gaussian processes.
#' @param kern A kernel function used to compute the covariance matrix at
#'    corresponding timestamps.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return The value of the penalised Gaussian elbo for
#'    the sum of the M individual GPs with common HPs.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
elbo_clust_multi_GP_common_hp_i <- function(hp, db, hyperpost, kern, pen_diag) {
  names_k <- hyperpost$mean %>% names()

  sum_i <- 0
  for (i in unique(db$ID))
  {
    ## Extract the i-th specific reference Input
    input_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Input)
    ## Extract the i-th specific inputs (reference + covariates)
    inputs_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(-c(.data$ID, .data$Output))
    ## Extract the i-th specific Inputs and Output
    output_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Output)

    corr1 <- 0
    corr2 <- 0

    for (k in (names_k))
    {
      ## Extract the covariance values associated with the i-th specific inputs
      post_cov_i <- hyperpost$cov[[k]][
        as.character(input_i),
        as.character(input_i)
      ]

      tau_i_k <- hyperpost$mixture %>%
        dplyr::filter(.data$ID == i) %>%
        dplyr::pull(k)
      mean_mu_k <- hyperpost$mean[[k]] %>%
        dplyr::filter(.data$Input %in% input_i) %>%
        dplyr::pull(.data$Output)
      corr1 <- corr1 + tau_i_k * mean_mu_k
      corr2 <- corr2 + tau_i_k *
        (mean_mu_k %*% t(mean_mu_k) + post_cov_i)
    }

    inv <- kern_to_inv(inputs_i, kern, hp, pen_diag)

    ## Classic Gaussian centred log-likelihood
    LL_norm <- -dmnorm(output_i, rep(0, length(output_i)), inv, log = T)

    sum_i <- sum_i +
      LL_norm - output_i %*% inv %*% corr1 + 0.5 * sum(inv * corr2)
  }
  return(sum_i)
}

#' Evidence Lower Bound maximised in MagmaClust
#'
#' @param hp_k A tibble, data frame or named vector of hyper-parameters
#'    for each clusters.
#' @param hp_i A tibble, data frame or named vector of hyper-parameters
#'    for each individuals.
#' @param db A tibble containing values we want to compute elbo on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param kern_i Kernel used to compute the covariance matrix of individuals GPs
#'    at corresponding inputs.
#' @param kern_k Kernel used to compute the covariance matrix of the mean GPs
#'    at corresponding inputs.
#' @param hyperpost A list of parameters for the variational distributions
#'    of the K mean GPs.
#' @param m_k Prior value of the mean parameter of the mean GPs (mu_k).
#'    Length = 1 or nrow(db).
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return Value of the elbo that is maximised during the VEM algorithm used for
#'    training in MagmaClust.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
elbo_monitoring_VEM <- function(hp_k,
                                hp_i,
                                db,
                                kern_i,
                                kern_k,
                                hyperpost,
                                m_k,
                                pen_diag) {
  floop <- function(k) {
    logL_GP_mod(
      hp_k[hp_k$ID == k, ],
      db = hyperpost$mean[[k]],
      mean = m_k[[k]],
      kern_k,
      hyperpost$cov[[k]],
      pen_diag
    ) %>%
      return()
  }
  sum_ll_k <- sapply(names(m_k), floop) %>% sum()

  floop2 <- function(i) {
    t_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Input)

    elbo_clust_multi_GP(
      hp_i[hp_i$ID == i, ],
      db %>% dplyr::filter(.data$ID == i),
      hyperpost,
      kern_i,
      pen_diag
    ) %>%
      return()
  }
  sum_ll_i <- sapply(unique(db$ID), floop2) %>% sum()

  floop3 <- function(k) {
    sum_tau <- 0
    det <- 0
    ## Extract the proportion in the k-th cluster
    pi_k <- hp_k %>%
      dplyr::filter(.data$ID == k) %>%
      dplyr::pull(.data$prop_mixture)

    for (i in unique(db$ID)) {
      ## Extract the probability of the i-th indiv to be in the k-th cluster
      tau_i_k <- hyperpost$mixture %>%
        dplyr::filter(.data$ID == i) %>%
        dplyr::pull(k)
      ## To avoid numerical issues if evaluating log(0/0)
      log_frac <- ifelse((tau_i_k == 0 | (pi_k == 0)), 0, log(pi_k / tau_i_k))

      sum_tau <- sum_tau + tau_i_k * log_frac
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

  return(-sum_ll_k - sum_ll_i + sum_corr_k)
}
