#' Compute the Multivariate Gaussian likelihood
#'
#' Modification of the function \code{dmvnorm()} from the package
#' \code{mvtnorm}, providing an implementation of the Multivariate Gaussian
#' likelihood. This version uses inverse of the covariance function as argument
#' instead of the traditional covariance.
#'
#' @param x A vector, containing values the likelihood is evaluated on.
#' @param mu A vector or matrix, specifying the mean parameter.
#' @param inv_Sigma A matrix, specifying the inverse of covariance parameter.
#' @param log A logical value, indicating whether we return the log-likelihood.
#'
#' @return A number, corresponding to the Multivariate Gaussian log-likelihood.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
dmnorm <- function(x, mu, inv_Sigma, log = FALSE) {
  if (is.vector(x)) {
    x <- t(as.matrix(x))
  }
  n <- length(mu)
  if (is.vector(mu)) {
    p <- length(mu)
    if (is.matrix(x)) {
      mu <- matrix(rep(mu, nrow(x)), ncol = p, byrow = TRUE)
    }
  } else {
    p <- ncol(mu)
  }
  if (!all(dim(inv_Sigma) == c(p, p)) || nrow(x) != nrow(mu)) {
    stop("incompatible arguments")
  }

  z <- t(x - mu)
  logdetS <- try(-determinant(inv_Sigma, logarithm = TRUE)$modulus,
    silent = TRUE
  )
  attributes(logdetS) <- NULL

  ssq <- t(z) %*% inv_Sigma %*% z
  loglik <- (-(n * (log(2 * pi)) + logdetS + ssq) / 2) %>% as.vector()
  if (log) {
    return(loglik)
  } else {
    return(exp(loglik))
  }
}

#' Log-Likelihood function of a Gaussian Process
#'
#' @param hp A tibble, data frame or named vector containing hyper-parameters.
#' @param db A tibble containing the values we want to compute the logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mean A vector, specifying the mean of the GP at the reference inputs.
#' @param kern A kernel function.
#' @param post_cov (optional) A matrix, corresponding to covariance parameter of
#'    the hyper-posterior. Used to compute the hyper-prior distribution of a new
#'    individual in Magma.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A number, corresponding to the value of Gaussian
#'    log-Likelihood (where the covariance can be the sum of the individual and
#'    the hyper-posterior's mean process covariances).
#'
#' @keywords internal
#'
#' @examples
#' TRUE
logL_GP <- function(hp,
                    db,
                    mean,
                    kern,
                    post_cov,
                    pen_diag) {

  ## Extract the input variables (reference Input + Covariates)
  inputs <- db %>% dplyr::select(-.data$Output)

  ## Sum the two covariance matrices and inverse the result
  cov <- kern_to_cov(inputs, kern, hp) + post_cov

  inv <- cov %>% chol_inv_jitter(pen_diag = pen_diag)

  (-dmnorm(db$Output, mean, inv, log = T)) %>%
    return()
}


#' Modified log-Likelihood function for GPs in a multi-output context
#'
#' Log-Likelihood function used during the maximization step of the training.
#' The log-Likelihood is defined as a Gaussian likelihood with a correction
#' trace term. This version is adapted for multi-output GPs.
#'
#' @param hp A numeric vector of hyper-parameters, as provided by `stats::optim`.
#' @param db A tibble containing the data to compute the logL on.
#'   Required columns: `Output` and `Output_ID`, plus input coordinates.
#' @param mean A vector specifying the mean of the GP at the reference inputs.
#' @param kern The kernel function (e.g., `convolution_kernel`).
#' @param post_cov A matrix, the covariance parameter of the hyper-posterior,
#'   used for the correction term.
#' @param pen_diag A jitter term added to the covariance matrix for numerical
#'   stability.
#' @param hp_col_names A character vector with the names of the hyper-parameters
#'   (e.g., c("l_t", "S_t")).
#' @param output_ids A character vector with the unique IDs of the outputs for
#'   the current task.
#'
#' @return A number, corresponding to the value of the modified Gaussian
#'   log-Likelihood.
#'
#' @keywords internal
logL_GP_mod <- function(hp,
                        db,
                        mean,
                        kern,
                        post_cov,
                        pen_diag,
                        hp_col_names,
                        output_ids) {
  # browser()
  if(!(hp %>% is_tibble()) && length(output_ids) > 1){
    # 1. Reconstruct the structured HP tibble from the flat vector
    hp_tibble <- reconstruct_hp(
      par_vector = hp,
      hp_names = hp_col_names,
      output_ids = output_ids
    )
  } else if (!(hp %>% is_tibble()) && length(output_ids) == 1){
    # browser()
    hp_tibble <- hp %>%
      t() %>%
      tibble::as_tibble() %>%
      stats::setNames(hp_col_names)

  } else {
    hp_tibble <- hp
  }

  list_output_ID <- db$Output_ID %>% unique()

  # 2. Build and invert the full multi-outputs covariance matrix
  if(length(list_output_ID) > 1 && !(kern %>% is.character())){
    # MO inversion of the TASK covariance
    # It will handle the multi-outputs structure and the noise addition internally.
    # 'kern_t' is expected to be the 'convolution_kernel' function.
    K_task_t <- kern_to_cov(
      input = db %>% dplyr::select(-Output),
      kern = kern,
      hp = hp_tibble
    )

    # Inverse K_task_t
    inv <- K_task_t %>%
      chol_inv_jitter(pen_diag = pen_diag) %>%
      `rownames<-`(row.names(K_task_t)) %>%
      `colnames<-`(colnames(K_task_t))

  } else if (length(list_output_ID) > 1 && (kern %>% is.character())){
    # MO inversion of the MEAN PROCESS covariance
    ## Compute the inverse covariance of the mean process
    inv <- ini_inverse_prior_cov(db, kern, hp, pen_diag)

  } else {
      # Single output case (for both computing the inverse of the task covariance
      # AND the inverse of the post_mean in logL_GP_monitoring())
      # Extract all_inputs to call kern_to_inv() on the single output case
      all_inputs <- db %>%
        dplyr::select(-c(Output, Output_ID)) %>%
        unique() %>%
        dplyr::arrange(Reference)

      if("Task_ID" %in% colnames(all_inputs)){
        all_inputs <- all_inputs %>% select(-Task_ID)
      }

      # Compute the inverse covariance matrix
      inv <- kern_to_inv(
        input = all_inputs,
        kern = kern,
        hp = hp_tibble,
        pen_diag = pen_diag
      )
  }
  # 3. Compute the log-likelihood components
  # Classical Gaussian log-likelihood
  LL_norm <- -dmnorm(db$Output, mean, inv, log = TRUE)
  # Correction trace term (-0.5 * Tr(inv %*% post_cov))
  # cor_term <- 0.5 * sum(inv %*% post_cov)
  cor_term <- 0.5 * sum(diag(inv %*% post_cov))

  return(LL_norm + cor_term)
}


#' Modified log-Likelihood with shared HPs for multi-task GPs
#'
#' Computes the sum of modified log-likelihoods over all tasks, assuming that
#' the hyper-parameters are shared across all tasks.
#'
#' @param hp A numeric vector of hyper-parameters, as provided by `stats::optim`.
#' @param db A tibble containing the data for all tasks.
#'   Required columns: `ID` (task ID), `Output`, `Output_ID`, plus inputs.
#' @param mean A vector, specifying the mean of the GP at the reference inputs.
#' @param kern The kernel function (e.g., `convolution_kernel`).
#' @param post_cov A matrix, covariance parameter of the hyper-posterior.
#' @param pen_diag A jitter term for numerical stability.
#' @param hp_col_names A character vector with the names of the hyper-parameters.
#' @param output_ids A character vector with the unique IDs of the outputs,
#'   assumed to be the same for all tasks.
#'
#' @return A number, the total modified log-Likelihood across all tasks.
#'
#' @keywords internal
logL_GP_mod_shared_tasks <- function(hp,
                                    db,
                                    mean,
                                    kern,
                                    post_cov,
                                    pen_diag,
                                    hp_col_names,
                                    output_ids) {
  # Loop over each task ID to compute and sum its log-likelihood
  funloop <- function(t) {
    # Extract data specific to task 't'
    db_t <- db %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::select(-Task_ID)

    input_t <- db_t %>% dplyr::pull(Reference) # Assumes 'Reference' col exists

    mean_t <- mean %>%
      dplyr::filter(Reference %in% input_t) %>%
      dplyr::pull(Output)

    post_cov_t <- post_cov[as.character(input_t), as.character(input_t)]

    # Call the single-task function, passing all HP-related arguments through
    logL_GP_mod(
      hp = hp,
      db = db_t,
      mean = mean_t,
      kern = kern,
      post_cov = post_cov_t,
      pen_diag = pen_diag,
      hp_col_names = hp_col_names,
      output_ids = output_ids
    )
  }

  sapply(unique(db$Task_ID), funloop) %>%
    sum() %>%
    return()
}


#' Log-Likelihood for monitoring the EM algorithm in Magma
#'
#' @param hp_0 A tibble or data frame, containing the hyper-parameters with the
#'    mean GP.
#' @param hp_t A tibble or data frame, containing the hyper-parameters with the
#'    task GPs.
#' @param db A tibble or data frame. Columns required: ID, Input, Output.
#'    Additional columns for covariates can be specified.
#' @param m_0 A vector, corresponding to the prior mean of the mean GP.
#' @param kern_0 A kernel function, associated with the mean process GP.
#' @param kern_t A kernel function, associated with the task GPs.
#' @param post_mean A tibble, coming out of the E step, containing the Input and
#'    associated Output of the hyper-posterior mean parameter.
#' @param post_cov A matrix, coming out of the E step, being the hyper-posterior
#'    covariance parameter.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A number, expectation of joint log-likelihood of the model. This
#'    quantity is supposed to increase at each step of the EM algorithm, and
#'    thus used for monitoring the procedure.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
logL_monitoring <- function(hp_0,
                            hp_t,
                            db,
                            m_0,
                            kern_0,
                            kern_t,
                            post_mean,
                            post_cov,
                            pen_diag) {

  # browser()
  ## Compute the modified logL for the mean process
  # browser()
  ll_0 <- logL_GP_mod(
    hp = hp_0,
    db = post_mean,
    mean = m_0,
    kern = kern_0,
    post_cov = post_cov,
    pen_diag = pen_diag
  )

  ## Sum over the tasks
  funloop <- function(t) {
    ## Extract the t-th specific reference inputs
    input_t <- db %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::pull(.data$Reference)
    ## Extract the t-th specific hyper-parameters
    hp_t_t <- hp_t %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::select(-Task_ID)
    ## Extract the t-th specific Inputs and Output
    db_t <- db %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::select(-Task_ID)
    ## Extract the mean values associated with the t-th specific inputs
    post_mean_t <- post_mean %>%
      dplyr::filter(Reference %in% input_t) %>%
      dplyr::pull(Output)
    ## Extract the covariance values associated with the t-th specific inputs
    post_cov_t <- post_cov[as.character(input_t),
                           as.character(input_t)
    ]

    ## Compute the modified logL for the individual processes
    logL_GP_mod(hp = hp_t_t,
                db = db_t,
                mean = post_mean_t,
                kern = kern_t,
                post_cov = post_cov_t,
                pen_diag = pen_diag) %>%
      return()
  }
  sum_ll_t <- sapply(unique(db$Task_ID), funloop) %>% sum()

  ## Compute the log-determinant term using Cholesky decomposition
  ## log(det(A)) = 2*sum(log(diag(chol(A))))
  det <- post_cov %>%
    chol() %>%
    diag() %>%
    log() %>%
    sum()

  ## Since the logL_GP_* functions return negative likelihoods for minimisation
  ## in the M-step, we need to x(-1) once more to retrieve the correct logL
  return(-ll_0 - sum_ll_t + det)
}

#' Compute a mixture of Gaussian log-likelihoods
#'
#' During the prediction step of MagmaClust, an EM algorithm is used to compute
#' the maximum likelihood estimator of the hyper-parameters along with
#' mixture probabilities for the new individual/task. This function implements
#' the quantity that is maximised (i.e. a sum of Gaussian log-likelihoods,
#' weighted by their mixture probabilities). It can also be used to monitor the
#' EM algorithm when providing the 'prop_mixture' argument, for proper
#' penalisation of the full log-likelihood.
#'
#' @param hp A tibble, data frame or named vector of hyper-parameters.
#' @param db A tibble containing data we want to evaluate the logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mixture A tibble or data frame, indicating the mixture probabilities
#'    of each cluster for the new individual/task.
#' @param mean A list of hyper-posterior mean parameters for all clusters.
#' @param kern A kernel function.
#' @param post_cov A list of hyper-posterior covariance parameters for all
#'    clusters.
#' @param prop_mixture A tibble or a named vector. Each name of column or
#'    element should refer to a cluster. The value associated with each cluster
#'    is a number between 0 and 1, corresponding to the mixture
#'    proportions.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A number, expectation of mixture of Gaussian log-likelihoods in
#'    the prediction step of MagmaClust. This quantity is supposed to increase
#'    at each step of the EM algorithm, and can be used for monitoring the
#'    procedure.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
sum_logL_GP_clust <- function(hp,
                              db,
                              mixture,
                              mean,
                              kern,
                              post_cov,
                              prop_mixture = NULL,
                              pen_diag) {

  ## Extract the observed (reference) Input
  input_obs <- db %>%
    dplyr::arrange(.data$Reference) %>%
    dplyr::pull(.data$Reference)

  ## Remove 'ID' if present in 'db'
  if ("ID" %in% names(db)) {
    db <- db %>% dplyr::select(-.data$ID)
  }

  ## Loop over the K clusters
  floop <- function(k) {
    tau_k <- mixture[[k]]
    mean_k <- mean[[k]] %>%
      dplyr::filter(.data$Reference %in% input_obs) %>%
      dplyr::pull(.data$Output)

    cov_k <- post_cov[[k]][
      as.character(input_obs),
      as.character(input_obs)
    ]

    sum_LL <- (tau_k * logL_GP(hp, db, mean_k, kern, cov_k, pen_diag))

    ## If prop_mixture is provided, compute full likelihood for monitoring EM
    if (!is.null(prop_mixture)) {
      pi_k <- prop_mixture[[k]]
      ## To avoid numerical issues if evaluating log(0/0)
      log_frac <- ifelse((tau_k == 0 | (pi_k == 0)), 0, log(pi_k / tau_k))

      ## Return -sum_LL because its a quantity that is minimised otherwise
      sum_LL <- -sum_LL + tau_k * log_frac
    }

    return(sum_LL)
  }
  sapply(names(mean), floop) %>%
    sum() %>%
    return()
}


#' Reconstruct the HP tibble from the optimizer's numeric vector
#'
#' This function handles two cases:
#' 1. A simple case where hp_names are like c("l_t", "S_t") and the vector
#'    is reshaped into a table with one row per output.
#' 2. A multi-output "flattened" case where hp_names are like
#'    c("l_t_1", "S_t_1", "l_t_2", "S_t_2", "l_u_t"). The function reshapes
#'    this into a long format tibble.
#'
#' @param par_vector The flat numeric vector handled by `stats::optim`.
#' @param hp_names A character vector containing the HP names.
#' @param output_ids A vector containing the unique output IDs.
#'
#' @return A correctly formatted HP tibble.
reconstruct_hp <- function(par_vector, hp_names, output_ids) {
  browser()
  # --- Détection du format des hyper-paramètres ---
  # On vérifie si au moins un nom contient un suffixe comme "_1", "_2", etc.
  is_flattened_format <- any(grepl("_\\d+$", hp_names))

  if (is_flattened_format) {
    # --- CAS 1 : Format aplati (ex: "l_t_1", "S_t_1", "l_u_t") ---

    shared_params <- hp_names[!grepl("_\\d+$", hp_names)]

    initial_tibble <- tibble::as_tibble_row(par_vector) %>%
      stats::setNames(hp_names)

    hp_tibble <- initial_tibble %>%
      tidyr::pivot_longer(
        cols = -all_of(shared_params),
        names_to = c(".value", "Output_ID"),
        names_pattern = "(.*)_(\\d+)$"
      ) %>%
      dplyr::mutate(Output_ID = as.character(Output_ID))
  } else { # Single output case
    hp_tibble <- par_vector %>%
      t() %>%
      tibble::as_tibble()
  }

  return(hp_tibble)
}
