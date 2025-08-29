#' Compute the Covariance Matrix for a Multi-Output GP via Convolution
#'
#' @param x The input data. For the vectorized multi-output case, this must be
#'   a tibble/data.frame containing the coordinates and an 'Output_ID' column.
#' @param y The second input data (must have the same format as x).
#' @param hp For the vectorized multi-output case, this must be a tibble
#'   containing the hyperparameters 'l_t', 'S_t', 'l_u_t' and an 'Output_ID'
#'   column.
#' @param vectorized If TRUE, enables the calculation of the full MO covariance
#'   matrix.
#' @param deriv A character string specifying the partial derivative to compute.
#'   Can be one of "l_t_1", "l_t_2", "S_t_1", "S_t_2", "l_u_t".
#' @return The covariance matrix or its partial derivative.
#'
#'
convolution_kernel <- function(x,
                               y,
                               hp,
                               vectorized = FALSE,
                               deriv = NULL) {
  # browser()
  # Le bloc non-vectorisé est principalement pour les calculs point par point
  if (!vectorized) {
    l_i <- exp(hp[["lengthscale_output1"]])
    l_j <- exp(hp[["lengthscale_output2"]])
    l_u <- exp(hp[["lengthscale_u"]])
    S_i <- exp(hp[["variance_output1"]])
    S_j <- exp(hp[["variance_output2"]])

    distance_sq     <- sum((x - y)^2)
    denominator_sum <- l_i + l_j + l_u
    S_prod          <- S_i * S_j

    pre_factor <- S_prod / sqrt(2 * pi * denominator_sum)
    exp_term <- exp(-0.5 * distance_sq / denominator_sum)
    K <- pre_factor * exp_term

  } else {
    # --- Bloc vectorisé pour construire la matrice de covariance complète ---
    if (!"Output_ID" %in% names(x) || !"Output_ID" %in% names(hp)) {
      stop("'input' and 'hp' must contain an 'Output_ID' column for vectorized mode.")
    }

    # Préparation des données et des hyperparamètres
    coord_cols <- names(x)[sapply(x, is.numeric) & names(x) != "Output_ID"]
    x_coords <- as.matrix(x[, coord_cols, drop = FALSE])
    y_coords <- as.matrix(y[, coord_cols, drop = FALSE])

    idx_x <- x$Output_ID
    idx_y <- y$Output_ID

    hp_ordered <- hp %>% dplyr::arrange(Output_ID)
    l_outputs <- hp_ordered$l_t
    S_outputs <- hp_ordered$S_t
    l_u_val <- exp(hp_ordered$l_u_t[1])

    l_vec_1 <- exp(l_outputs[idx_x])
    l_vec_2 <- exp(l_outputs[idx_y])
    S_vec_1 <- exp(S_outputs[idx_x])
    S_vec_2 <- exp(S_outputs[idx_y])

    # Calcul des composantes sous forme de matrices
    distance_sq     <- cpp_dist(x_coords, y_coords)
    denominator_sum <- outer(l_vec_1, l_vec_2, FUN = "+") + l_u_val
    S_prod          <- outer(S_vec_1, S_vec_2, FUN = "*")

    pre_factor <- S_prod / sqrt(2 * pi * denominator_sum)
    exp_term   <- exp(-0.5 * distance_sq / denominator_sum)
    K          <- pre_factor * exp_term
  }

  if (is.null(deriv)) {
    return(K)
  }

  # --- CALCUL DES DÉRIVÉES PARTIELLES (CORRIGÉ) ---

  hp_id_str <- stringr::str_extract(deriv, "\\d+$")
  if (is.na(hp_id_str) && deriv != "l_u_t") {
    stop("Le nom de la dérivée doit se terminer par un ID numérique, ex: 'l_t_1'.")
  }
  hp_id <- as.integer(hp_id_str)

  if (startsWith(deriv, "l_t_")) {
    # Formule de base (inchangée)
    common_deriv_denom <- K * ((-0.5 / denominator_sum) + (0.5 * distance_sq / (denominator_sum^2)))

    # Règle de dérivation en chaîne (inchangée)
    chain_rule_factor <- exp(l_outputs[hp_id])

    # --- CORRECTION : Création d'une matrice de facteurs ---
    N <- nrow(x)
    M <- nrow(y)
    mask_i_is_k <- outer(idx_x == hp_id, rep(TRUE, M))
    mask_j_is_k <- outer(rep(TRUE, N), idx_y == hp_id)

    # On additionne les masques booléens (TRUE=1, FALSE=0).
    # Sur un bloc diagonal (k,k), mask_i_is_k et mask_j_is_k sont tous deux TRUE,
    # donc leur somme vaut 2. Ailleurs, elle vaut 1 ou 0.
    factor_matrix <- mask_i_is_k + mask_j_is_k

    return(common_deriv_denom * chain_rule_factor * factor_matrix)
  } else if (startsWith(deriv, "S_t_")) {
    # La dérivée par rapport à log(S_k) est K * (I(i=k) + I(j=k))
    # où I est la fonction indicatrice.
    N <- nrow(x)
    M <- nrow(y)
    mask_i_is_k <- outer(idx_x == hp_id, rep(TRUE, M))
    mask_j_is_k <- outer(rep(TRUE, N), idx_y == hp_id)

    # L'addition des masques booléens (convertis en 0/1) implémente (I(i=k) + I(j=k))
    pd <- K * (mask_i_is_k + mask_j_is_k)
    return(pd)

  } else if (deriv == "l_u_t") {
    # Formule corrigée pour la dépendance au dénominateur (avec D^2)
    common_deriv_denom <- K * ((-0.5 / denominator_sum) + (0.5 * distance_sq / (denominator_sum^2)))

    # Règle de dérivation en chaîne : on multiplie par l_u
    chain_rule_factor <- l_u_val

    # Pas de masque car l_u est partagé et affecte tous les blocs
    return(common_deriv_denom * chain_rule_factor)

  } else {
    stop("Invalid 'deriv' argument.")
  }
}


#' Compute the Covariance Matrix for a Multi-Output GP via Convolution
#'
#' @param x The input data. For the vectorized multi-output case, this must be
#'   a tibble/data.frame containing the coordinates and an 'Output_ID' column.
#' @param y The second input data (must have the same format as x).
#' @param hp For the vectorized multi-output case, this must be a tibble
#'   containing the hyperparameters 'l_t', 'S_t', 'l_u_t' and an 'Output_ID'
#'   column.
#' @param vectorized If TRUE, enables the calculation of the full MO covariance
#'   matrix.
#' @return The covariance matrix.

# convolution_kernel <- function(x,
#                                y,
#                                hp,
#                                vectorized = FALSE) {
#
#   # ----- NON-VECTORIZED CASE -----
#   # Handles the calculation between two single points.
#   if (!vectorized) {
#     l_1 <- exp(hp[["lengthscale_output1"]])
#     l_2 <- exp(hp[["lengthscale_output2"]])
#     l_u <- exp(hp[["lengthscale_u"]])
#     S_1 <- exp(hp[["variance_output1"]])
#     S_2 <- exp(hp[["variance_output2"]])
#     distance <- sum((x - y)^2)
#     top_term <- - (l_1 + l_2 + l_u)^(-1) * 0.5 * distance
#     return(((S_1 * S_2 / ((2*pi)^(0.5)*(l_1 + l_2 + l_u)^(0.5))) * exp(top_term)))
#   }
#
#   # ----- VECTORIZED CASE -----
#
#   # 1. Input checks
#   if (!"Output_ID" %in% names(x) || !"Output_ID" %in% names(hp)) {
#     stop("For the vectorized convolution kernel, 'input' and 'hp' must contain",
#          " an 'Output_ID' column.")
#   }
#
#   # 2. Separate coordinates from output indices
#   #    We exclude non-numeric columns to get the coordinates.
#   coord_cols <- names(x)[sapply(x, is.numeric) & names(x) != "Output_ID"]
#   x_coords <- as.matrix(x[, coord_cols, drop = FALSE])
#   y_coords <- as.matrix(y[, coord_cols, drop = FALSE])
#   idx_x <- x$Output_ID
#   idx_y <- y$Output_ID
#
#   # 3. Prepare hyperparameter vectors
#   #    Ensure HPs are correctly sorted by Output_ID to guarantee alignment.
#   hp_ordered <- hp %>% dplyr::arrange(Output_ID)
#
#   l_outputs <- hp_ordered$l_t
#   S_outputs <- hp_ordered$S_t
#   l_u <- hp_ordered$l_u_t[1] # Assume l_u is shared across outputs for a given task
#
#   # 4. Compute the covariance matrix (meta-kernel logic)
#   distance_matrix_sq <- cpp_dist(x_coords, y_coords)
#
#   l_vec_x <- exp(l_outputs[idx_x])
#   l_vec_y <- exp(l_outputs[idx_y])
#   S_vec_x <- exp(S_outputs[idx_x])
#   S_vec_y <- exp(S_outputs[idx_y])
#
#   l_sum_mat <- outer(l_vec_x, l_vec_y, FUN = "+")
#   S_prod_mat <- outer(S_vec_x, S_vec_y, FUN = "*")
#
#   denominator_sum <- l_sum_mat + exp(l_u)
#   top_term_mat <- -0.5 * distance_matrix_sq / denominator_sum
#   pre_factor_mat <- S_prod_mat / sqrt(2 * pi * denominator_sum)
#
#   return(pre_factor_mat * exp(top_term_mat))
# }


#' Squared Exponential Kernel
#'
#' @param x A vector (or matrix if vectorized = T) of inputs.
#' @param y A vector (or matrix if vectorized = T) of inputs.
#' @param hp A tibble, data frame or named vector, containing the kernel's
#'    hyperparameters. Required columns: 'se_variance', 'se_lengthscale'.
#' @param deriv A character, indicating according to which hyper-parameter the
#'    derivative should be computed. If NULL (default), the function simply
#'    returns the evaluation of the kernel.
#' @param vectorized A logical value, indicating whether the function provides
#'    a vectorized version for speeded-up calculations. If TRUE, the \code{x}
#'    and \code{y} arguments should be the vector or matrix containing all
#'    inputs for which the kernel is evaluated on all pairs of elements.
#'    If FALSE, the \code{x} and \code{y} arguments are simply two inputs.
#'
#' @return A scalar, corresponding to the evaluation of the kernel.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
se_kernel <- function(
  x,
  y,
  hp,
  deriv = NULL,
  vectorized = FALSE) {
  ## Check whether the Rcpp function for speed-up vectorised computation
  if (vectorized) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    ## Compute directly the full covariance matrix
    distance <- cpp_dist(x, y)
  } else {
    ## Compute one element of the covariance matrix
    distance <- sum((x - y)^2)
  }

  top_term <- exp(-hp[["se_lengthscale"]]) * 0.5 * distance

  if (deriv %>% is.null()) {
    (exp(hp[["se_variance"]] - top_term)) %>% return()
  } else if (deriv == "se_variance") {
    (exp(hp[["se_variance"]] - top_term)) %>% return()
  } else if (deriv == "se_lengthscale") {
    (exp(hp[["se_variance"]]) * top_term * exp(-top_term)) %>% return()
  } else {
    stop("Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'")
  }
}

#' Periodic Kernel
#'
#' @param x A vector (or matrix if vectorized = T) of inputs.
#' @param y A vector (or matrix if vectorized = T) of inputs.
#' @param hp A tibble, data frame or named vector, containing the kernel's
#'    hyperparameters. Required columns: 'perio_variance', 'perio_lengthscale',
#'    and 'period'.
#' @param deriv A character, indicating according to which hyper-parameter the
#'  derivative should be computed. If NULL (default), the function simply returns
#'  the evaluation of the kernel.
#' @param vectorized A logical value, indicating whether the function provides
#'    a vectorized version for speeded-up calculations. If TRUE, the \code{x}
#'    and \code{y} arguments should be the vector or matrix containing all
#'    inputs for which the kernel is evaluated on all pairs of elements.
#'    If FALSE, the \code{x} and \code{y} arguments are simply two inputs.
#'
#'
#' @return A scalar, corresponding to the evaluation of the kernel.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
perio_kernel <- function(
  x,
  y,
  hp,
  deriv = NULL,
  vectorized = FALSE) {
  ## Check whether the Rcpp function for speed-up vectorised computation
  if (vectorized) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    ## Compute directly the full covariance matrix
    perio_term <- cpp_perio(x, y, hp[["period"]])
    sum_deriv <- cpp_perio_deriv(x, y, hp[["period"]])
  } else {
    ## Compute one element of the covariance matrix
    angle <- pi * abs(x - y) / exp(hp[["period"]])
    perio_term <- sin(angle)^2 %>% sum()
    sum_deriv <- sum(2 * sin(angle) * cos(angle) * angle)
  }

  if (deriv %>% is.null()) {
    (exp(hp[["perio_variance"]]) *
      exp(-2 * exp(-hp[["perio_lengthscale"]]) * perio_term)) %>%
      return()
  } else if (deriv == "perio_variance") {
    (exp(hp[["perio_variance"]]) *
      exp(-2 * exp(-hp[["perio_lengthscale"]]) * perio_term)) %>%
      return()
  } else if (deriv == "period") {
    (exp(hp[["perio_variance"]]) * exp(-2 * exp(-hp[["perio_lengthscale"]]) *
      perio_term) * 2 * exp(-hp[["perio_lengthscale"]]) * sum_deriv) %>%
      return()
  } else if (deriv == "perio_lengthscale") {
    (exp(hp[["perio_variance"]]) * 2 * perio_term *
      exp(-hp[["perio_lengthscale"]]) *
      exp(-2 * perio_term * exp(-hp[["perio_lengthscale"]]))) %>%
      return()
  } else {
    stop("Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'")
  }
}

#' Rational Quadratic Kernel

#' @param x A vector (or matrix if vectorized = T) of inputs.
#' @param y A vector (or matrix if vectorized = T) of inputs.
#' @param hp A tibble, data frame or named vector, containing the kernel's
#'    hyperparameters. Required columns: 'rq_variance', 'rq_lengthscale', and
#'    'rq_scale'.
#' @param deriv A character, indicating according to which hyper-parameter the
#'  derivative should be computed. If NULL (default), the function simply returns
#'  the evaluation of the kernel.
#' @param vectorized A logical value, indicating whether the function provides
#'    a vectorized version for speeded-up calculations. If TRUE, the \code{x}
#'    and \code{y} arguments should be the vector or matrix containing all
#'    inputs for which the kernel is evaluated on all pairs of elements.
#'    If FALSE, the \code{x} and \code{y} arguments are simply two inputs.
#'
#' @return A scalar, corresponding to the evaluation of the kernel.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
rq_kernel <- function(
  x,
  y,
  hp,
  deriv = NULL,
  vectorized = FALSE) {
  ## Check whether the Rcpp function for speed-up vectorised computation
  if (vectorized) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    ## Compute directly the full covariance matrix
    distance <- cpp_dist(x, y)
  } else {
    ## Compute one element of the covariance matrix
    distance <- sum((x - y)^2)
  }

  term <- (1 + distance * exp(-hp[["rq_lengthscale"]]) / (2 * hp[["rq_scale"]]))

  if (deriv %>% is.null()) {
    (exp(hp[["rq_variance"]]) * term^(-hp[["rq_scale"]])) %>%
      return()
  } else if (deriv == "rq_variance") {
    (exp(hp[["rq_variance"]]) * term^(-hp[["rq_scale"]])) %>%
      return()
  } else if (deriv == "rq_scale") {
    (exp(hp[["rq_variance"]]) * term^(-hp[["rq_scale"]]) *
      (distance * exp(-hp[["rq_lengthscale"]]) / (2 * hp[["rq_scale"]] * term)
        - log(term))) %>%
      return()
  } else if (deriv == "rq_lengthscale") {
    (exp(hp[["rq_variance"]]) * distance * 0.5 *
      exp(-hp[["rq_lengthscale"]]) * term^(-hp[["rq_scale"]] - 1)) %>%
      return()
  } else {
    stop("Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'")
  }
}

#' Linear Kernel
#'
#' @param x A vector (or matrix if vectorized = T) of inputs.
#' @param y A vector (or matrix if vectorized = T) of inputs.
#' @param hp A tibble, data frame or named vector, containing the kernel's
#'    hyperparameters. Required columns: 'lin_slope' and 'lin_offset'.
#' @param deriv A character, indicating according to which hyper-parameter the
#'    derivative should be computed. If NULL (default), the function simply
#'    returns the evaluation of the kernel.
#' @param vectorized A logical value, indicating whether the function provides
#'    a vectorized version for speeded-up calculations. If TRUE, the \code{x}
#'    and \code{y} arguments should be the vector or matrix containing all
#'    inputs for which the kernel is evaluated on all pairs of elements.
#'    If FALSE, the \code{x} and \code{y} arguments are simply two inputs.
#'
#' @return A scalar, corresponding to the evaluation of the kernel.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
lin_kernel <- function(
  x,
  y,
  hp,
  deriv = NULL,
  vectorized = FALSE) {
  ## Check whether the Rcpp function for speed-up vectorised computation
  if (vectorized) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    ## Compute directly the full covariance matrix
    prod <- cpp_prod(x, y)
    ## Create a dummy matrix of one for the 'lin_offset' derivative matrix
    mat_one <- matrix(1, ncol = ncol(prod), nrow = nrow(prod))
  } else {
    ## Compute one element of the covariance matrix
    prod <- x %*% y
    mat_one <- 1
  }

  if (deriv %>% is.null()) {
    (exp(hp[["lin_slope"]]) * prod + exp(hp[["lin_offset"]])) %>%
      return()
  } else if (deriv == "lin_offset") {
    exp(hp[["lin_offset"]] * mat_one) %>% return()
  } else if (deriv == "lin_slope") {
    (exp(hp[["lin_slope"]]) * prod) %>% return()
  } else {
    stop("Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'")
  }
}


#' Generate Initial Hyperparameters for GP Kernels
#'
#' A helper function to generate a tibble of random initial hyperparameter (HP)
#' values for common Gaussian Process (GP) kernels. It can generate a single
#' set of HPs, create HPs for different tasks or outputs, or for every
#' combination of task and output.
#'
#' @param kern A character string (e.g., "SE", "SE + LIN")
#'   or a function object (e.g., `convolution_kernel`).
#' @param list_task_ID A vector of task or individual IDs.
#' @param list_output_ID A vector of output IDs.
#' @param shared_hp_tasks If TRUE, all tasks share the same hyperparameter
#'   values.
#' @param shared_hp_outputs If TRUE, all outputs within a task share the same
#'   values for l_t, S_t and noise.
#' @param noise A boolean. If TRUE, an additional 'noise' hyperparameter is
#'   added to the tibble. Default is FALSE.
#'
#' @return A `tibble` containing the generated hyperparameters.

hp <- function(kern = "SE",
               list_task_ID = NULL,
               list_output_ID = NULL,
               shared_hp_tasks = TRUE,
               shared_hp_outputs = TRUE,
               noise = FALSE) {

  # Define boundaries for random sampling
  min_val <- -2; max_val <- 2; min_noise <- -5; max_noise <- -1

  #=============================================================================
  # Case 1: Kernel is provided as a function object (for convolution_kernel)
  #=============================================================================
  if (is.function(kern)) {
    kern_name <- deparse(substitute(kern))
    if (kern_name != "convolution_kernel") {
      stop("Currently, only 'convolution_kernel' is supported as a function input.")
    }
    if (is.null(list_task_ID) || is.null(list_output_ID)) {
      stop("For the convolution_kernel, both 'list_task_ID' and 'list_output_ID' must be provided.")
    }

    num_tasks <- length(list_task_ID)
    num_outputs <- length(list_output_ID)

    # Create the base grid of all Task x Output combinations
    base_ids <- tidyr::crossing(Task_ID = as.character(list_task_ID),
                                Output_ID = as.character(list_output_ID))

    # --- Handle the four sharing scenarios for l_t, S_t, and noise ---

    if (shared_hp_tasks && shared_hp_outputs) {
      # SCENARIO A: 1 set of HPs for EVERYTHING (all tasks and outputs)
      hps_to_add <- tibble::tibble(
        l_t = stats::runif(1, min_val, max_val),
        S_t = stats::runif(1, min_val, max_val)
      )
      if (noise) hps_to_add$noise <- stats::runif(1, min_noise, max_noise)

      final_hp <- tidyr::crossing(base_ids, hps_to_add)

    } else if (shared_hp_tasks && !shared_hp_outputs) {
      # SCENARIO B: HPs are shared across tasks, but different for each output.
      hps_per_output <- tibble::tibble(
        Output_ID = as.character(list_output_ID),
        l_t = stats::runif(num_outputs, min_val, max_val),
        S_t = stats::runif(num_outputs, min_val, max_val)
      )
      if (noise) hps_per_output$noise <- stats::runif(num_outputs, min_noise, max_noise)

      final_hp <- dplyr::left_join(base_ids, hps_per_output, by = "Output_ID")

    } else if (!shared_hp_tasks && shared_hp_outputs) {
      # SCENARIO C: HPs are different for each task, but shared by outputs within a task.
      hps_per_task <- tibble::tibble(
        Task_ID = as.character(list_task_ID),
        l_t = stats::runif(num_tasks, min_val, max_val),
        S_t = stats::runif(num_tasks, min_val, max_val)
      )
      if (noise) hps_per_task$noise <- stats::runif(num_tasks, min_noise, max_noise)

      final_hp <- dplyr::left_join(base_ids, hps_per_task, by = "Task_ID")

    } else { # !shared_hp_tasks && !shared_hp_outputs
      # SCENARIO D: Every Task x Output combination gets its own unique HPs.
      n_draws <- nrow(base_ids)
      hps_to_add <- tibble::tibble(
        l_t = stats::runif(n_draws, min_val, max_val),
        S_t = stats::runif(n_draws, min_val, max_val)
      )
      if (noise) hps_to_add$noise <- stats::runif(n_draws, min_noise, max_noise)

      final_hp <- dplyr::bind_cols(base_ids, hps_to_add)
    }

    # --- Handle l_u_t (latent process lengthscale), which only depends on tasks ---
    if (shared_hp_tasks) {
      # If tasks share HPs, they share one l_u_t
      final_hp$l_u_t <- stats::runif(1, -2, 0.001)
    } else {
      # If tasks have different HPs, they get different l_u_t values
      l_u_t_per_task <- tibble::tibble(
        Task_ID = as.character(list_task_ID),
        l_u_t   = stats::runif(length(list_task_ID), -2, 0.001)
      )
      final_hp <- dplyr::left_join(final_hp, l_u_t_per_task, by = "Task_ID")
    }

    return(final_hp)

  } else {
    #===========================================================================
    # Case 2: Kernel is provided as a string (e.g., "SE", "SE + PERIO")
    # This logic remains as it was in your original function.
    #===========================================================================
    base_ids <- NULL
    n_draws <- 1

    if (!is.null(list_task_ID) && !is.null(list_output_ID)) {
      base_ids <- tidyr::crossing(Task_ID = as.character(list_task_ID),
                                  Output_ID = as.character(list_output_ID))
      if (!shared_hp_tasks && !shared_hp_outputs) {
        n_draws <- nrow(base_ids)
      } else if (!shared_hp_tasks && shared_hp_outputs) {
        n_draws <- length(list_task_ID)
      } else if (shared_hp_tasks && !shared_hp_outputs) {
        n_draws <- length(list_output_ID)
      } # else n_draws remains 1
    } else if (!is.null(list_output_ID)) {
      base_ids <- tibble::tibble(Output_ID = as.character(list_output_ID))
      n_draws <- if (shared_hp_outputs) 1 else nrow(base_ids)
    } else if (!is.null(list_task_ID)) {
      base_ids <- tibble::tibble(Task_ID = as.character(list_task_ID))
      n_draws <- if (shared_hp_tasks) 1 else nrow(base_ids)
    }

    str_kern <- strsplit(kern, " +")[[1]]

    # Generate the required number of unique HP sets
    generated_hps <- tibble::tibble(.rows = n_draws)
    for (i in str_kern) {
      temp_hp <- switch(i,
                        "SE" = tibble::tibble(
                          se_variance = stats::runif(n_draws, min_val, max_val),
                          se_lengthscale = stats::runif(n_draws, min_val, max_val)
                        ),
                        "PERIO" = tibble::tibble(
                          perio_variance = stats::runif(n_draws, min_val, max_val),
                          perio_lengthscale = stats::runif(n_draws, min_val, max_val),
                          period = stats::runif(n_draws, 0, 2 * pi)
                        ),
                        "RQ" = tibble::tibble(
                          rq_variance = stats::runif(n_draws, min_val, max_val),
                          rq_lengthscale = stats::runif(n_draws, min_val, max_val),
                          rq_scale = stats::runif(n_draws, min_val, max_val)
                        ),
                        "LIN" = tibble::tibble(
                          lin_slope = stats::runif(n_draws, min_val, max_val),
                          lin_offset = stats::runif(n_draws, min_val, max_val)
                        ),
                        # Default case for operators like '+' or '*'
                        tibble::tibble()
      )
      generated_hps <- dplyr::bind_cols(generated_hps, temp_hp)
    }

    if (noise) {
      generated_hps <- generated_hps %>%
        dplyr::mutate(noise = stats::runif(n_draws, min_noise, max_noise))
    }

    # Construct the Final Tibble by combining IDs and HPs
    if (is.null(base_ids)) {
      return(generated_hps)
    }

    if (!is.null(list_task_ID) && !is.null(list_output_ID)) {
      if (!shared_hp_tasks && !shared_hp_outputs) {
        final_hp <- dplyr::bind_cols(base_ids, generated_hps)
      } else if (!shared_hp_tasks && shared_hp_outputs) {
        task_hps <- dplyr::bind_cols(Task_ID = as.character(list_task_ID), generated_hps)
        final_hp <- dplyr::left_join(base_ids, task_hps, by = "Task_ID")
      } else if (shared_hp_tasks && !shared_hp_outputs) {
        output_hps <- dplyr::bind_cols(Output_ID = as.character(list_output_ID), generated_hps)
        final_hp <- dplyr::left_join(base_ids, output_hps, by = "Output_ID")
      } else {
        final_hp <- tidyr::crossing(base_ids, generated_hps)
      }
    } else { # Case for only one list provided (task or output)
      if (n_draws == 1) { # Shared HPs
        final_hp <- tidyr::crossing(base_ids, generated_hps)
      } else { # Independent HPs
        final_hp <- dplyr::bind_cols(base_ids, generated_hps)
      }
    }

    return(final_hp)
  }
}
