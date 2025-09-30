#' Gradient of the logLikelihood of a Gaussian Process
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
#' @return A named vector, corresponding to the value of the hyper-parameters
#'    gradients for the Gaussian log-Likelihood (where the covariance can be the
#'    sum of the individual and the hyper-posterior's mean process covariances).
#'
#' @keywords internal
#'
#' @examples
#' TRUE
gr_GP <- function(hp,
                  db,
                  mean,
                  kern,
                  post_cov,
                  pen_diag,
                  hp_col_names) {
  if(!(hp %>% is_tibble()) && length(db$Output_ID %>% unique()) > 1){
    # 1. Reconstruct the structured HP tibble from the flat vector
    hp_tibble <- reconstruct_hp(
      par_vector = hp,
      hp_names = hp_col_names,
      output_ids = db$Output_ID %>% unique()
    )
  } else if (!(hp %>% is_tibble()) && length(db$Output_ID %>% unique()) == 1){
    hp_tibble <- hp %>%
      t() %>%
      tibble::as_tibble() %>%
      stats::setNames(hp_col_names)

  } else {
    hp_tibble <- hp
  }

  list_ID_outputs <- db$Output_ID %>% unique()

  # 2. Build and invert the full multi-output covariance matrix
  #    'inputs' must contain the 'Output_ID' column for the kernel to work
  if(length(list_ID_outputs) > 1 && !(kern %>% is.character())){
    # MO inversion of the TASK covariance
    # Call kern_to_cov directly.
    # It will handle the multi-output structure and the noise addition internally.
    # 'kern_t' is expected to be the 'convolution_kernel' function.
    cov <- kern_to_cov(
      input = db %>%
        dplyr::select(Output_ID, dplyr::starts_with("Input")),
      kern = kern,
      hp = hp_tibble
    ) + post_cov

    # Inverse K_task_t
    inv <- cov %>% chol_inv_jitter(pen_diag = pen_diag)

  } else {
    # Single output case
    # Extract all_inputs to call kern_to_cov() on the single output case
    all_inputs <- db %>%
      dplyr::select(-c(Output, Output_ID)) %>%
      unique() %>%
      dplyr::arrange(Reference)

    # Compute the inverse covariance matrix of the task 't'
    inv <- kern_to_inv(
      input = all_inputs,
      kern = kern,
      hp = hp_tibble,
      pen_diag = pen_diag
    )
  }

  output <- db$Output

  ## Compute the term common to all partial derivatives
  prod_inv <- inv %*% (output - mean)
  common_term <- prod_inv %*% t(prod_inv) - inv

  floop <- function(deriv_name) {
    # Get the derivative of the covariance matrix w.r.t. the current HP
    if(length(list_ID_outputs) > 1){
      dcov <- kern_to_cov(db %>%
                            dplyr::select(Output_ID, dplyr::starts_with("Input")),
                          kern = kern,
                          hp = hp_tibble,
                          deriv = deriv_name)
    } else {
      # Extract all_inputs to call kern_to_cov() on the single output case
      all_inputs <- db %>%
        dplyr::select(-c(Output, Output_ID)) %>%
        unique() %>%
        dplyr::arrange(Reference)

      dcov <- kern_to_cov(input = all_inputs,
                          kern = kern,
                          hp = hp_tibble,
                          deriv = deriv_name)
    }

    # Compute the gradient component for this HP
    gradient_value <- -0.5 * (common_term %*% dcov) %>%
      diag() %>%
      sum() %>%
      return()
  }

  sapply(hp_col_names, floop) %>%
    return()
}


#' Gradient of the modified logLikelihood for GPs in a multi-output context
#'
#' @param hp A numeric vector of hyper-parameters, as provided by `stats::optim`.
#' @param db A tibble containing the data to compute the gradient on.
#'   Required columns: `Output`, `Output_ID`, plus input coordinates.
#' @param mean A vector specifying the mean of the GP at the reference inputs.
#' @param kern The kernel function (e.g., `convolution_kernel`).
#' @param post_cov A matrix, the covariance parameter of the hyper-posterior.
#' @param pen_diag A jitter term for numerical stability.
#' @param hp_col_names A character vector with the names of the hyper-parameters.
#' @param output_ids A character vector with the unique IDs of the outputs.
#' @param priors A list or tibble containing prior information, e.g.,
#'   list(l_t = list(dist = "invgamma", shape = 2, scale = 1), ...).
#'   If NULL, performs ML estimation.
#'
#' @return A named vector of gradients for each hyper-parameter.
#'
#' @keywords internal
gr_GP_mod <- function(hp,
                      db,
                      mean,
                      kern,
                      post_cov,
                      pen_diag,
                      hp_col_names,
                      output_ids,
                      priors,
                      ...) {
  bro
  if(!(hp %>% is_tibble()) && length(output_ids) > 1){
    # 1. Reconstruct the structured HP tibble from the flat vector
    hp_tibble <- reconstruct_hp(
      par_vector = hp,
      hp_names = hp_col_names,
      output_ids = output_ids
    )
  } else if (!(hp %>% is_tibble()) && length(output_ids) == 1){
    hp_tibble <- hp %>%
      t() %>%
      tibble::as_tibble() %>%
      stats::setNames(hp_col_names)

  } else {
    hp_tibble <- hp
  }

  list_ID_outputs <- db$Output_ID %>% unique()
  # 2. Build and invert the full multi-output covariance matrix
  #    'inputs' must contain the 'Output_ID' column for the kernel to work
  if(length(list_ID_outputs) > 1 && !(kern %>% is.character())){
    # MO inversion of the TASK covariance
    # Call kern_to_cov directly.
    # It will handle the multi-output structure and the noise addition internally.
    # 'kern_t' is expected to be the 'convolution_kernel' function.
    K_task_t <- kern_to_cov(
      input = db %>%
                dplyr::select(Output_ID, dplyr::starts_with("Input")),
      kern = kern,
      hp = hp_tibble
    )

    # Inverse K_task_t
    inv <- K_task_t %>% chol_inv_jitter(pen_diag = pen_diag)

  } else if (length(list_ID_outputs) > 1 && kern %>% is.character()){
    # MO inversion of the MEAN PROCESS covariance
    ## Compute the inverse covariance of the mean process
    inv <- ini_inverse_prior_cov(db, kern, hp, pen_diag)

  } else{
    # Single output case
    # Extract all_inputs to call kern_to_cov() on the single output case
    all_inputs <- db %>%
      dplyr::select(-c(Output, Output_ID)) %>%
      unique() %>%
      dplyr::arrange(Reference)

    if("Task_ID" %in% colnames(all_inputs)){
      all_inputs <- all_inputs %>% select(-Task_ID)
    }

    # Compute the inverse covariance matrix of the task 't'
    inv <- kern_to_inv(
      input = all_inputs,
      kern = kern,
      hp = hp_tibble,
      pen_diag = pen_diag
    )
  }

  # 3. Compute the term common to all partial derivatives
  prod_inv <- inv %*% (db$Output - mean)

  common_term <- prod_inv %*% t(prod_inv) +
    inv %*% (post_cov %*% inv - diag(1, nrow(db)))

  # 4. Loop over HPs to compute the gradient for each
  # The derivative names must match what the kernel expects in its 'deriv'
  ## argument
  floop <- function(deriv_name) {
    # Get the derivative of the covariance matrix w.r.t. the current HP
    if(length(list_ID_outputs) > 1){
      dK_dhp <- kern_to_cov(db %>%
                              dplyr::select(Output_ID, dplyr::starts_with("Input")),
                            kern = kern,
                            hp = hp_tibble,
                            deriv = deriv_name)
    } else{
      # Extract all_inputs to call kern_to_cov() on the single output case
      all_inputs <- db %>%
        dplyr::select(-c(Output, Output_ID)) %>%
        unique() %>%
        dplyr::arrange(Reference)

      if("Task_ID" %in% colnames(all_inputs)){
        all_inputs <- all_inputs %>% select(-Task_ID)
      }

      dK_dhp <- kern_to_cov(input = all_inputs,
                            kern = kern,
                            hp = hp_tibble,
                            deriv = deriv_name)
    }

    # Compute the gradient component for this HP
    gradient_value <- -0.5 * (common_term %*% dK_dhp) %>%
      diag() %>%
      sum()
  }

  # Return a named vector of gradients
  neg_grad_Q  <- sapply(hp_col_names, floop, USE.NAMES = TRUE)

  # 5. Add the negative gradient of the log-prior for MAP
  neg_grad_log_prior <- rep(0, length(neg_grad_Q))
  names(neg_grad_log_prior) <- names(neg_grad_Q)

  if (length(priors)!=0) {
    # On itère sur chaque nom de paramètre du vecteur d'optimisation (ex: "l_t_1")
    for (param_name in names(neg_grad_log_prior)) {

      # --- ÉTAPE 1 : Extraire le nom de la colonne (ex: "l_t_1" -> "l_t") ---
      col_name <- gsub("_\\d+$", "", param_name)

      # Vérifier si un prior est défini pour ce type de paramètre
      if (col_name %in% names(priors)) {

        # --- ÉTAPE 2 : Extraire l'index de la ligne (ex: "l_t_1" -> 1) ---
        # On utilise une expression régulière pour trouver le nombre à la fin
        output_index_str <- sub(".*_(\\d+)$", "\\1", param_name)

        # Si un nombre est trouvé, on le convertit. Sinon (cas de "l_u_t"), on prend la 1ère ligne.
        row_index <- if (grepl("\\d+", output_index_str)) {
          as.numeric(output_index_str)
        } else {
          1
        }

        # --- ÉTAPE 3 : Récupérer la valeur unique de l'HP ---
        current_hp_value <- hp_tibble[[col_name]][row_index]

        # --- ÉTAPE 4 : Calculer le gradient comme avant ---
        prior_info <- priors[[col_name]]
        grad_log_p <- 0

        if (prior_info$dist == "invgamma") {
          alpha <- prior_info$shape
          beta <- prior_info$scale
          grad_log_p <- - (alpha + 1) / current_hp_value + beta / (current_hp_value^2)
        } else if (prior_info$dist == "gamma") {
          alpha <- prior_info$shape
          beta <- prior_info$rate
          grad_log_p <- (alpha - 1) / current_hp_value - beta
        } else if (prior_info$dist == "normal") {
          mu <- prior_info$mean
          sigma <- prior_info$sd
          grad_log_p <- - (current_hp_value - mu) / (sigma^2)
        }

        if(is.na(grad_log_p) || is.infinite(grad_log_p)) grad_log_p <- 0

        # --- ÉTAPE 5 : Assigner le gradient à la bonne place ---
        neg_grad_log_prior[param_name] <- -grad_log_p
      }
    }
  }

  final_gradient <- neg_grad_Q + neg_grad_log_prior
  return(final_gradient)
}



#' Gradient of the logLikelihood with shared HPs for multi-task GPs
#'
#' Computes the sum of gradients over all tasks, assuming that the
#' hyper-parameters are shared across all tasks.
#'
#' @param hp A numeric vector of hyper-parameters, as provided by `stats::optim`.
#' @param db A tibble containing the data for all tasks.
#'   Required columns: `ID` (task ID), `Output`, `Output_ID`, plus inputs.
#' @param mean A vector, specifying the mean of the GP at the reference inputs.
#' @param kern The kernel function (e.g., `convolution_kernel`).
#' @param post_cov A matrix, covariance parameter of the hyper-posterior.
#' @param pen_diag A jitter term for numerical stability.
#' @param hp_col_names A character vector with the names of the hyper-parameters.
#' @param output_ids A character vector with the unique IDs of the outputs.
#' @param priors A list or tibble containing prior information, e.g.,
#'   list(l_t = list(dist = "invgamma", shape = 2, scale = 1), ...).
#'   If NULL, performs ML estimation.
#'
#' @return A named vector, the total gradient across all tasks.
#'
#' @keywords internal
gr_GP_mod_shared_tasks <- function(hp,
                                db,
                                mean,
                                kern,
                                post_cov,
                                pen_diag,
                                hp_col_names,
                                output_ids,
                                priors) {

  # Loop over each task ID to compute its gradient vector
  funloop <- function(t) {
    # Extract data specific to task 't'
    db_t <- db %>% dplyr::filter(Task_ID == t) %>%
                   dplyr::select(-Task_ID)

    input_t <- db_t %>% dplyr::pull(Reference)

    mean_t <- mean %>%
      dplyr::filter(Reference %in% input_t) %>%
      dplyr::pull(Output)

    post_cov_t <- post_cov[as.character(input_t), as.character(input_t)]
    # browser()
    # Call the single-task gradient function, passing all arguments through
    gr_GP_mod(
      hp = hp,
      db = db_t,
      mean = mean_t,
      kern = kern,
      post_cov = post_cov_t,
      pen_diag = pen_diag,
      hp_col_names = hp_col_names,
      output_ids = output_ids,
      priors = priors
    )

  }

  # Sum the gradient vectors from all tasks element-wise
  sapply(unique(db$Task_ID), funloop) %>%
    rowSums() %>%
    return()
}


#'  Gradient of the mixture of Gaussian likelihoods
#'
#' Compute the gradient of a sum of Gaussian log-likelihoods, weighted by their
#' mixture probabilities.
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
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A named vector, corresponding to the value of the hyper-parameters'
#'    gradients for the mixture of Gaussian log-likelihoods involved in the
#'    prediction step of MagmaClust.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
gr_sum_logL_GP_clust <- function(hp,
                                 db,
                                 mixture,
                                 mean,
                                 kern,
                                 post_cov,
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
    (tau_k * gr_GP(hp, db, mean_k, kern, cov_k, pen_diag)) %>%
      return()
  }
  sapply(names(mean), floop) %>%
    rowSums() %>%
    return()
}
