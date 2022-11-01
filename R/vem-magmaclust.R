#' E-Step of the VEM algorithm
#'
#' Expectation step of the Variational EM algorithm used to compute
#' the parameters of the hyper-posteriors distributions
#' for the mean processes and mixture variables involved in MagmaClust.
#'
#' @param db A tibble or data frame. Columns required: ID, Input, Output.
#'    Additional columns for covariates can be specified.
#' @param m_k A named list of vectors, corresponding to the prior mean
#'    parameters of the K mean GPs.
#' @param kern_k A kernel function, associated with the K mean GPs.
#' @param kern_i A kernel function, associated with the M individual GPs.
#' @param hp_k A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_k}.
#' @param hp_i A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}.
#' @param old_mixture A list of mixture values from the previous iteration.
#' @param iter A number, indicating the current iteration of the VEM algorithm.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#'
#' @return A named list, containing the elements \code{mean}, a tibble
#'    containing the Input and associated Output of the hyper-posterior mean
#'    parameters, \code{cov}, the hyper-posterior covariance matrices,
#'    and \code{mixture}, the probabilities to belong to each cluster for each
#'    individual.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
ve_step <- function(db,
                    m_k,
                    kern_k,
                    kern_i,
                    hp_k,
                    hp_i,
                    old_mixture,
                    iter,
                    pen_diag) {

  ## Extract the union of all reference inputs provided in the training data
  all_inputs <- db %>%
    dplyr::select(-.data$ID,-.data$Output) %>%
    unique() %>%
    dplyr::arrange(.data$Reference)

  all_input <- all_inputs %>% dplyr::pull(.data$Reference)

  ## Sort the database according to Reference
  db <- db %>% dplyr::arrange(.data$Reference, .by_group = TRUE)

  prop_mixture_k <- hp_k %>%
    dplyr::pull(.data$prop_mixture, name = .data$ID)

  ## Format a sequence of inputs for all clusters
  t_clust <- tidyr::expand_grid("ID" = names(m_k),
                                all_inputs
  )

  ## Compute all the inverse covariance matrices
  list_inv_k <- list_kern_to_inv(t_clust, kern_k, hp_k, pen_diag)
  list_inv_i <- list_kern_to_inv(db, kern_i, hp_i, pen_diag)

  ## Create a named list of Output values for all individuals
  list_output_i <- base::split(db$Output, list(db$ID))

  ## Update each mu_k parameters for each cluster ##
  floop <- function(k) {
    post_inv <- list_inv_k[[k]]
    tau_k <- old_mixture %>% dplyr::select(.data$ID, k)
    for (i in list_inv_i %>% names())
    {
      ## Extract the corresponding mixture probability
      tau_i_k <- tau_k %>%
        dplyr::filter(.data$ID == i) %>%
        dplyr::pull(k)

      inv_i <- list_inv_i[[i]]
      ## Collect input's common indices between mean and individual processes
      co_input <- intersect(row.names(inv_i), row.names(post_inv))
      ## Sum the common inverse covariance's terms
      post_inv[co_input, co_input] <- post_inv[co_input, co_input] +
        tau_i_k * inv_i[co_input, co_input]
    }

    post_inv %>%
      chol_inv_jitter(pen_diag = pen_diag) %>%
      `rownames<-`(all_input) %>%
      `colnames<-`(all_input) %>%
      return()
  }
  cov_k <- sapply(tidyselect::all_of(names(m_k)),
                  floop,
                  simplify = FALSE,
                  USE.NAMES = TRUE)

  ## Update the posterior mean for each cluster ##

  floop2 <- function(k) {
    prior_mean <- m_k[[k]]
    prior_inv <- list_inv_k[[k]]
    tau_k <- old_mixture %>% dplyr::select(.data$ID, k)

    weighted_mean <- prior_inv %*% prior_mean

    for (i in list_inv_i %>% names())
    {
      ## Extract the corresponding mixture probability
      tau_i_k <- tau_k %>%
        dplyr::filter(.data$ID == i) %>%
        dplyr::pull(k)
      ## Compute the weighted mean for the i-th individual
      weighted_i <- tau_i_k * list_inv_i[[i]] %*% list_output_i[[i]]
      ## Collect input's common indices between mean and individual processes
      co_input <- intersect(row.names(weighted_i), row.names(weighted_mean))
      ## Sum the common weighted mean's terms
      weighted_mean[co_input, ] <- weighted_mean[co_input, ] +
        weighted_i[co_input, ]
    }

    ## Compute the updated mean parameter
    new_mean <- cov_k[[k]] %*% weighted_mean %>% as.vector()
    tibble::tibble(all_inputs,
                   "Output" = new_mean) %>% return()
  }
  mean_k <- sapply(names(m_k), floop2, simplify = FALSE, USE.NAMES = TRUE)

  ## Update mixture (skip first iteration to avoid bad HP initialisation issues)
  if(iter == 1){
    mixture <- old_mixture
  } else{
    mixture <- update_mixture(
      db,
      mean_k,
      cov_k,
      hp_i,
      kern_i,
      prop_mixture_k,
      pen_diag
    )
  }

  list(
    "mean" = mean_k,
    "cov" = cov_k,
    "mixture" = mixture
  ) %>%
    return()
}


#' V-Step of the VEM algorithm
#'
#' Maximization step of the Variational EM algorithm used to compute
#' hyper-parameters of all the kernels involved in MagmaClust.
#'
#' @param db A tibble or data frame. Columns required: ID, Input, Output.
#'    Additional columns for covariates can be specified.
#' @param list_mu_param List of parameters of the K mean GPs.
#' @param kern_k A kernel used to compute the covariance matrix of the mean GP
#'    at corresponding timestamps.
#' @param kern_i A kernel used to compute the covariance matrix of individuals
#'    GP at corresponding timestamps.
#' @param m_k A named list of prior mean parameters for the K mean GPs.
#'    Length = 1 or nrow(unique(db$Input))
#' @param common_hp_k A boolean indicating whether hp are common among
#'    mean GPs (for each mu_k)
#' @param common_hp_i A boolean indicating whether hp are common among
#'    individual GPs (for each y_i)
#' @param old_hp_i A named vector, tibble or data frame, containing the
#'    hyper-parameters from the previous  M-step (or initialisation) associated
#'    with the individual GPs.
#' @param old_hp_k A named vector, tibble or data frame, containing the
#'    hyper-parameters from the previous M-step (or initialisation) associated
#'    with the mean GPs.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#' numerical issues when inverting nearly singular matrices.
#'
#' @return A named list, containing the elements \code{hp_k}, a tibble
#'    containing the hyper-parameters associated with each cluster,
#'    \code{hp_i}, a tibble containing the hyper-parameters
#'    associated with the individual GPs, and \code{prop_mixture_k},
#'    a tibble containing the hyper-parameters associated with each individual,
#'    indicating the probabilities to belong to each cluster.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
vm_step <- function(db,
                    old_hp_k,
                    old_hp_i,
                    list_mu_param,
                    kern_k,
                    kern_i,
                    m_k,
                    common_hp_k,
                    common_hp_i,
                    pen_diag) {

  list_ID_k <- names(m_k)
  list_ID_i <- unique(db$ID)

  list_hp_i <- old_hp_i %>%
    dplyr::select(-.data$ID) %>%
    names()

  list_hp_k <- old_hp_k %>%
    dplyr::select(-.data$ID) %>%
    dplyr::select(-.data$prop_mixture) %>%
    names()

  ## Detect whether kernel_k provides derivatives for its hyper-parameters
  if (kern_k %>% is.function()) {
    if (!("deriv" %in% methods::formalArgs(kern_k))) {
      gr_GP_mod <- NULL
      gr_GP_mod_common_hp_k <- NULL
    }
  }

  ## Detect whether kernel_i provides derivatives for its hyper-parameters
  if (kern_i %>% is.function()) {
    if (!("deriv" %in% methods::formalArgs(kern_i))) {
      gr_clust_multi_GP_common_hp_i <- NULL
      gr_clust_multi_GP <- NULL
    }
  }

  ## Check whether hyper-parameters are common to all individuals
  if (common_hp_i) {
    ## Extract the hyper-parameters associated with the i-th individual
    par_i <- old_hp_i %>%
      dplyr::select(-.data$ID) %>%
      dplyr::slice(1)

    ## Optimise hyper-parameters of the individual processes
    new_hp_i <- stats::optim(
      par = par_i,
      fn = elbo_clust_multi_GP_common_hp_i,
      gr = gr_clust_multi_GP_common_hp_i,
      db = db,
      hyperpost = list_mu_param,
      kern = kern_i,
      pen_diag = pen_diag,
      method = "L-BFGS-B",
      control = list(factr = 1e13, maxit = 25)
    )$par %>%
      tibble::as_tibble_row() %>%
      tidyr::uncount(weights = length(list_ID_i)) %>%
      dplyr::mutate("ID" = list_ID_i, .before = 1)
  } else {
    loop2 <- function(i) {
      ## Extract the hyper-parameters associated with the i-th individual
      par_i <- old_hp_i %>%
        dplyr::filter(.data$ID == i) %>%
        dplyr::select(-.data$ID)
      ## Extract the data associated with the i-th individual
      db_i <- db %>% dplyr::filter(.data$ID == i)

      ## Optimise hyper-parameters of the individual processes
      stats::optim(
        par = par_i,
        fn = elbo_clust_multi_GP,
        gr = gr_clust_multi_GP,
        db = db_i,
        pen_diag = pen_diag,
        hyperpost = list_mu_param,
        kern = kern_i,
        method = "L-BFGS-B",
        control = list(factr = 1e13, maxit = 25)
      )$par %>%
        tibble::as_tibble_row() %>%
        return()
    }
    new_hp_i <- sapply(list_ID_i, loop2, simplify = FALSE, USE.NAMES = TRUE) %>%
      tibble::enframe(name = "ID") %>%
      tidyr::unnest(cols = .data$value)
  }

  ## Compute the prop mixture of each cluster
  prop_mixture <- list_mu_param$mixture %>%
    dplyr::select(-.data$ID) %>%
    colMeans()

  ## Check whether hyper-parameters are common to all cluster
  if (common_hp_k) {
    ## Extract the hyper-parameters associated with the k-th cluster
    par_k <- old_hp_k %>%
      dplyr::select(-.data$ID) %>%
      dplyr::slice(1) %>%
      dplyr::select(-.data$prop_mixture)

    ## Optimise hyper-parameters of the processes of each cluster
    new_hp_k <- stats::optim(
      par = par_k,
      fn = elbo_GP_mod_common_hp_k,
      gr = gr_GP_mod_common_hp_k,
      db = list_mu_param$mean,
      mean = m_k,
      kern = kern_k,
      post_cov = list_mu_param$cov,
      pen_diag = pen_diag,
      method = "L-BFGS-B",
      control = list(factr = 1e13, maxit = 25)
    )$par %>%
      tibble::as_tibble_row() %>%
      tidyr::uncount(weights = length(list_ID_k)) %>%
      dplyr::mutate("ID" = list_ID_k, .before = 1) %>%
      dplyr::mutate("prop_mixture" = prop_mixture)
  } else {
    loop <- function(k) {
      ## Extract the hyper-parameters associated with the k-th cluster
      par_k <- old_hp_k %>%
        dplyr::filter(.data$ID == k) %>%
        dplyr::select(-.data$ID) %>%
        dplyr::select(-.data$prop_mixture)
      ## Extract the data associated with the k-th cluster
      db_k <- list_mu_param$mean[[k]]
      ## Extract the mean values associated with the k-th specific inputs
      mean_k <- m_k[[k]]
      ## Extract the covariance values associated with the k-th specific inputs
      post_cov_k <- list_mu_param$cov[[k]]

      ## Optimise hyper-parameters of the processes of each cluster
      stats::optim(
        par = par_k,
        logL_GP_mod,
        gr = gr_GP_mod,
        db = db_k,
        mean = mean_k,
        kern = kern_k,
        post_cov = post_cov_k,
        pen_diag = pen_diag,
        method = "L-BFGS-B",
        control = list(factr = 1e13, maxit = 25)
      )$par %>%
        tibble::as_tibble_row() %>%
        return()
    }
    new_hp_k <- sapply(list_ID_k, loop, simplify = FALSE, USE.NAMES = TRUE) %>%
      tibble::enframe(name = "ID") %>%
      tidyr::unnest_wider(.data$value) %>%
      dplyr::mutate("prop_mixture" = prop_mixture)
  }

  list(
    "hp_k" = new_hp_k,
    "hp_i" = new_hp_i
  ) %>%
    return()
}

#' Update the mixture probabilities for each individual and each cluster
#'
#' @param db A tibble or data frame. Columns required: \code{ID},
#'    \code{Input}, \code{Output}. Additional columns for covariates can be
#'    specified.
#' @param mean_k A list of the K hyper-posterior mean parameters.
#' @param cov_k A list of the K hyper-posterior covariance matrices.
#' @param hp A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern}, the individual process' kernel. The
#'    columns/elements should be named according to the hyper-parameters
#'    that are used in \code{kern}.
#' @param kern A kernel function, defining the covariance structure of
#'    the individual GPs.
#' @param prop_mixture A tibble containing the hyper-parameters associated
#'    with each individual, indicating in which cluster it belongs.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#' numerical issues when inverting nearly singular matrices.
#'
#' @return Compute the hyper-posterior multinomial distributions by updating
#'    mixture probabilities.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
update_mixture <- function(db,
                           mean_k,
                           cov_k,
                           hp,
                           kern,
                           prop_mixture,
                           pen_diag) {
  c_i <- 0
  c_k <- 0
  ID_i <- unique(db$ID)
  ID_k <- names(mean_k)
  mat_elbo <- matrix(NA, nrow = length(ID_k), ncol = length(ID_i))
  vec_prop <- c()

  for (i in ID_i)
  {
    c_i <- c_i + 1
    ## Extract the i-th specific Input
    input_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Reference)
    ## Extract the i-th specific hyper-parameters
    hp_i <- hp %>%
      dplyr::filter(.data$ID == i)
    ## Extract the data associated with the i-th individual
    db_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(-.data$ID)

    for (k in ID_k)
    {
      c_k <- c_k + 1

      ## Create a vector of proportion with the clusters in adequate order
      vec_prop[c_k] <- prop_mixture[[k]]
      ## Extract the mean values associated with the i-th specific inputs
      mean_k_i <- mean_k[[k]] %>%
        dplyr::filter(.data$Reference %in% input_i) %>%
        dplyr::pull(.data$Output)
      ## Extract the covariance values associated with the i-th specific inputs
      cov_k_i <- cov_k[[k]][as.character(input_i), as.character(input_i)]

      mat_elbo[c_k, c_i] <- -logL_GP_mod(
        hp_i,
        db_i,
        mean_k_i,
        kern,
        cov_k_i,
        pen_diag
      )
    }
    c_k <- 0
  }

  ## We need to use the 'log-sum-exp' trick: exp(x - max(x))/sum exp(x - max(x))
  ## to remain numerically stable
  mat_L <- mat_elbo %>% apply(2, function(x) exp(x - max(x)))

  (vec_prop * mat_L) %>%
    apply(2, function(x) x / sum(x)) %>%
    `rownames<-`(ID_k) %>%
    t() %>%
    round(5) %>%
    tibble::as_tibble() %>%
    dplyr::mutate("ID" = ID_i, .before = 1) %>%
    return()
}
