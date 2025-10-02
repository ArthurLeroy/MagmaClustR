#' @param hp For the vectorized multi-output case, this must be a tibble
#'   containing the hyperparameters 'p_t', 'S_t', 'p_u_t' and an 'Output_ID'
#'   column.
#' @param deriv A character string specifying the partial derivative to compute.
#'   Can be one of "p_t_1", "p_t_2", "S_t_1", "S_t_2", "p_u_t".
#'
convolution_kernel <- function(x,
                               y,
                               hp,
                               vectorized = FALSE,
                               deriv = NULL) {

  if (!vectorized) {
    # Ce bloc est moins critique mais pourrait être mis à jour par cohérence
    # Pour l'instant, nous nous concentrons sur le cas vectorisé qui est utilisé
    # par l'optimisation et la prédiction.
    stop("Le mode non-vectorisé n'est pas mis à jour pour la paramétrisation en précision.")
  }

  # --- Bloc vectorisé pour construire la matrice de covariance complète ---
  if (!"Output_ID" %in% names(x) || !"Output_ID" %in% names(hp)) {
    stop("'input' and 'hp' must contain an 'Output_ID' column for vectorized mode.")
  }

  x_ids_numeric <- as.numeric(as.character(x$Output_ID))
  y_ids_numeric <- as.numeric(as.character(y$Output_ID))
  hp_ids_numeric <- as.numeric(as.character(hp$Output_ID))

  coord_cols <- names(x)[sapply(x, is.numeric) & names(x) != "Output_ID"]
  x_coords <- as.matrix(x[, coord_cols, drop = FALSE])
  y_coords <- as.matrix(y[, coord_cols, drop = FALSE])

  indices_x <- match(x_ids_numeric, hp_ids_numeric)
  indices_y <- match(y_ids_numeric, hp_ids_numeric)

  # NOUVELLE PARAMÉTRISATION : p_t et p_u_t sont les log-précisions
  # longueur_echelle = exp(-p / 2)
  l_vec_1 <- exp(-hp$p_t[indices_x] / 2)
  l_vec_2 <- exp(-hp$p_t[indices_y] / 2)
  S_vec_1 <- exp(hp$S_t[indices_x])
  S_vec_2 <- exp(hp$S_t[indices_y])
  l_u_val <- exp(-hp$p_u_t[1] / 2)

  # Le reste du calcul de la covariance est inchangé
  distance_sq     <- cpp_dist(x_coords, y_coords)
  denominator_sum <- outer(l_vec_1, l_vec_2, FUN = "+") + l_u_val
  S_prod          <- outer(S_vec_1, S_vec_2, FUN = "*")

  pre_factor <- S_prod / sqrt(2 * pi * denominator_sum)
  exp_term   <- exp(-0.5 * distance_sq / denominator_sum)
  K          <- pre_factor * exp_term

  if (is.null(deriv)) {
    return(K)
  }

  # --- Calcul des dérivées partielles ---

  hp_id_str <- stringr::str_extract(deriv, "\\d+$")
  if (is.na(hp_id_str) && deriv != "p_u_t") {
    stop("The name of the derivative should end with a number, ex: 'p_t_1'.")
  }
  hp_id <- as.integer(hp_id_str)

  # Dérivée par rapport à la log-précision p_t_k
  if (startsWith(deriv, "p_t_")) {
    # Terme commun de la dérivée par rapport au dénominateur
    common_deriv_denom <- K * ((-0.5 / denominator_sum) + (0.5 * distance_sq / (denominator_sum^2)))

    # Facteur de la règle de dérivation (chain rule)
    # d(denom)/dp_k = d(l_i)/dp_k + d(l_j)/dp_k
    # l_k = exp(-p_k / 2) => dl_k/dp_k = -0.5 * exp(-p_k/2) = -0.5 * l_k
    current_l_k <- exp(-hp$p_t[hp_ids_numeric == hp_id] / 2)
    chain_rule_factor <- -0.5 * current_l_k

    N <- nrow(x)
    M <- nrow(y)
    mask_i_is_k <- outer(x_ids_numeric == hp_id, rep(TRUE, M))
    mask_j_is_k <- outer(rep(TRUE, N), y_ids_numeric == hp_id)
    factor_matrix <- mask_i_is_k + mask_j_is_k

    return(common_deriv_denom * chain_rule_factor * factor_matrix)

    # La dérivée par rapport à la variance S_t ne change pas
  } else if (startsWith(deriv, "S_t_")) {
    N <- nrow(x)
    M <- nrow(y)
    mask_i_is_k <- outer(x_ids_numeric == hp_id, rep(TRUE, M))
    mask_j_is_k <- outer(rep(TRUE, N), y_ids_numeric == hp_id)

    pd <- K * (mask_i_is_k + mask_j_is_k)
    return(pd)

    # Dérivée par rapport à la log-précision p_u_t
  } else if (deriv == "p_u_t") {
    common_deriv_denom <- K * ((-0.5 / denominator_sum) + (0.5 * distance_sq / (denominator_sum^2)))

    # Facteur de la règle de dérivation pour p_u_t
    current_l_u <- exp(-hp$p_u_t[1] / 2)
    chain_rule_factor <- -0.5 * current_l_u

    return(common_deriv_denom * chain_rule_factor)

  } else {
    stop("Invalid 'deriv' argument.")
  }
}


#' #' Compute the Covariance Matrix for a Multi-Output GP via Convolution
#' #'
#' #' @param x The input data. For the vectorized multi-output case, this must be
#' #'   a tibble/data.frame containing the coordinates and an 'Output_ID' column.
#' #' @param y The second input data (must have the same format as x).
#' #' @param hp For the vectorized multi-output case, this must be a tibble
#' #'   containing the hyperparameters 'l_t', 'S_t', 'l_u_t' and an 'Output_ID'
#' #'   column.
#' #' @param vectorized If TRUE, enables the calculation of the full MO covariance
#' #'   matrix.
#' #' @param deriv A character string specifying the partial derivative to compute.
#' #'   Can be one of "l_t_1", "l_t_2", "S_t_1", "S_t_2", "l_u_t".
#' #' @return The covariance matrix or its partial derivative.
#' #'
#' #'
#' convolution_kernel <- function(x,
#'                                y,
#'                                hp,
#'                                vectorized = FALSE,
#'                                deriv = NULL) {
#'   # browser()
#'   # Le bloc non-vectorisé est principalement pour les calculs point par point
#'   if (!vectorized) {
#'     l_i <- exp(hp[["lengthscale_output1"]])
#'     l_j <- exp(hp[["lengthscale_output2"]])
#'     l_u <- exp(hp[["lengthscale_u"]])
#'     S_i <- exp(hp[["variance_output1"]])
#'     S_j <- exp(hp[["variance_output2"]])
#'
#'     distance_sq     <- sum((x - y)^2)
#'     denominator_sum <- l_i + l_j + l_u
#'     S_prod          <- S_i * S_j
#'
#'     pre_factor <- S_prod / sqrt(2 * pi * denominator_sum)
#'     exp_term <- exp(-0.5 * distance_sq / denominator_sum)
#'     K <- pre_factor * exp_term
#'
#'   } else {
#'     # --- Bloc vectorisé pour construire la matrice de covariance complète ---
#'     if (!"Output_ID" %in% names(x) || !"Output_ID" %in% names(hp)) {
#'       stop("'input' and 'hp' must contain an 'Output_ID' column for vectorized mode.")
#'     }
#'
#'     x_ids_numeric <- as.numeric(as.character(x$Output_ID))
#'     y_ids_numeric <- as.numeric(as.character(y$Output_ID))
#'     hp_ids_numeric <- as.numeric(as.character(hp$Output_ID))
#'
#'     coord_cols <- names(x)[sapply(x, is.numeric) & names(x) != "Output_ID"]
#'     x_coords <- as.matrix(x[, coord_cols, drop = FALSE])
#'     y_coords <- as.matrix(y[, coord_cols, drop = FALSE])
#'
#'     # Recherche robuste avec les IDs maintenant numériques
#'     indices_x <- match(x_ids_numeric, hp_ids_numeric)
#'     indices_y <- match(y_ids_numeric, hp_ids_numeric)
#'
#'     l_vec_1 <- exp(hp$l_t[indices_x])
#'     l_vec_2 <- exp(hp$l_t[indices_y])
#'     S_vec_1 <- exp(hp$S_t[indices_x])
#'     S_vec_2 <- exp(hp$S_t[indices_y])
#'     l_u_val <- exp(hp$l_u_t[1])
#'
#'     # Calcul des composantes sous forme de matrices
#'     distance_sq     <- cpp_dist(x_coords, y_coords)
#'     denominator_sum <- outer(l_vec_1, l_vec_2, FUN = "+") + l_u_val
#'     S_prod          <- outer(S_vec_1, S_vec_2, FUN = "*")
#'
#'     pre_factor <- S_prod / sqrt(2 * pi * denominator_sum)
#'     exp_term   <- exp(-0.5 * distance_sq / denominator_sum)
#'     K          <- pre_factor * exp_term
#'   }
#'
#'   if (is.null(deriv)) {
#'     return(K)
#'   }
#'
#'   # Partial derivative computation
#'
#'   hp_id_str <- stringr::str_extract(deriv, "\\d+$")
#'   if (is.na(hp_id_str) && deriv != "l_u_t") {
#'     stop("The name of the derivative should end with a number, ex: 'l_t_1'.")
#'   }
#'   hp_id <- as.integer(hp_id_str)
#'
#'   if (startsWith(deriv, "l_t_")) {
#'     common_deriv_denom <- K * ((-0.5 / denominator_sum) + (0.5 * distance_sq / (denominator_sum^2)))
#'
#'     chain_rule_factor <- exp(hp$l_t[hp_ids_numeric == hp_id])
#'
#'     N <- nrow(x)
#'     M <- nrow(y)
#'     mask_i_is_k <- outer(x_ids_numeric == hp_id, rep(TRUE, M))
#'     mask_j_is_k <- outer(rep(TRUE, N), y_ids_numeric == hp_id)
#'
#'     factor_matrix <- mask_i_is_k + mask_j_is_k
#'
#'     return(common_deriv_denom * chain_rule_factor * factor_matrix)
#'
#'   } else if (startsWith(deriv, "S_t_")) {
#'     N <- nrow(x)
#'     M <- nrow(y)
#'     mask_i_is_k <- outer(x_ids_numeric == hp_id, rep(TRUE, M))
#'     mask_j_is_k <- outer(rep(TRUE, N), y_ids_numeric == hp_id)
#'
#'     pd <- K * (mask_i_is_k + mask_j_is_k)
#'     return(pd)
#'
#'   } else if (deriv == "l_u_t") {
#'     # Common_term
#'     common_deriv_denom <- K * ((-0.5 / denominator_sum) + (0.5 * distance_sq / (denominator_sum^2)))
#'
#'     chain_rule_factor <- exp(hp$l_u_t[1])
#'
#'     return(common_deriv_denom * chain_rule_factor)
#'
#'   } else {
#'     stop("Invalid 'deriv' argument.")
#'   }
#' }


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
  # browser()
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
#' @param kern A character string (e.g., "SE") or a function (e.g., `convolution_kernel`).
#' @param list_task_ID A vector of task IDs.
#' @param list_output_ID A vector of output IDs.
#' @param shared_hp_tasks If TRUE, HPs are shared across tasks.
#' @param shared_hp_outputs If TRUE, HPs are shared across outputs.
#' @param noise If TRUE, a 'noise' hyperparameter is added.
#' @param hp_config A tibble providing min/max bounds for HP generation, specific
#'   to the convolution_kernel. If NULL, default bounds are used.
#'
#' @return A `tibble` containing the generated hyperparameters.
#'

#' @param hp_config A tibble providing min/max bounds for HP generation, specific
#'   to the convolution_kernel. If NULL, default bounds are used.
#'
hp <- function(kern = "SE",
               list_task_ID = NULL,
               list_output_ID = NULL,
               shared_hp_tasks = TRUE,
               shared_hp_outputs = TRUE,
               noise = FALSE,
               hp_config = NULL) {
  ## Les valeurs initiales ne changent pas
  min_val <- -3
  max_val <- 3
  min_noise <- -5
  max_noise <- -1

  if (is.function(kern)) {
    kern_name <- deparse(substitute(kern))
    if (kern_name != "convolution_kernel") {
      stop("Currently, only 'convolution_kernel' is supported as a function input.")
    }
    if (is.null(list_task_ID) || is.null(list_output_ID)) {
      stop("For the convolution_kernel, both 'list_task_ID' and 'list_output_ID' must be provided.")
    }

    if (is.null(hp_config)) {
      message("hp_config not provided for convolution_kernel, using default HP bounds.")
      # MODIFICATION 1: Renommer les colonnes de 'l' en 'p'
      hp_config <- tibble::tibble(
        output_id   = list_output_ID,
        pt_min      = -2, pt_max      = 2,    # Changé de lt_min/lt_max
        St_min      = -2, St_max      = 2,
        pu_min      = -2, pu_max      = 0,    # Changé de lu_min/lu_max
        noise_min   = -5, noise_max   = -2
      )
    }

    num_tasks <- length(list_task_ID)
    num_outputs <- length(list_output_ID)

    base_ids <- tidyr::crossing(Task_ID = as.character(list_task_ID),
                                Output_ID = as.character(list_output_ID))

    if (shared_hp_tasks && shared_hp_outputs) {
      # MODIFICATION 2: Générer 'p_t' au lieu de 'l_t'
      hps_to_add <- tibble::tibble(
        p_t = stats::runif(1, hp_config$pt_min[1], hp_config$pt_max[1]), # Changé de l_t, lt_min, lt_max
        S_t = stats::runif(1, hp_config$St_min[1], hp_config$St_max[1])
      )
      if (noise) {
        hps_to_add$noise <- stats::runif(1, hp_config$noise_min[1], hp_config$noise_max[1])
      }
      final_hp <- tidyr::crossing(base_ids, hps_to_add)

    } else if (shared_hp_tasks && !shared_hp_outputs) {
      # MODIFICATION 3: Générer 'p_t' par sortie
      hps_per_output <- hp_config %>%
        dplyr::transmute(
          Output_ID = as.character(output_id),
          p_t = purrr::map2_dbl(pt_min, pt_max, ~stats::runif(1, .x, .y)), # Changé de l_t, lt_min, lt_max
          S_t = purrr::map2_dbl(St_min, St_max, ~stats::runif(1, .x, .y))
        )
      if (noise) {
        hps_per_output$noise <- purrr::map2_dbl(hp_config$noise_min,
                                                hp_config$noise_max,
                                                ~stats::runif(1, .x, .y))
      }
      final_hp <- dplyr::left_join(base_ids, hps_per_output, by = "Output_ID")

    } else if (!shared_hp_tasks && shared_hp_outputs) {
      # MODIFICATION 4: Générer 'p_t' par tâche
      hps_per_task <- tibble::tibble(
        Task_ID = as.character(list_task_ID),
        p_t = stats::runif(num_tasks, hp_config$pt_min[1], hp_config$pt_max[1]), # Changé de l_t, lt_min, lt_max
        S_t = stats::runif(num_tasks, hp_config$St_min[1], hp_config$St_max[1])
      )
      if (noise) {
        hps_per_task$noise <- stats::runif(num_tasks,
                                           hp_config$noise_min[1],
                                           hp_config$noise_max[1])
      }
      final_hp <- dplyr::left_join(base_ids, hps_per_task, by = "Task_ID")

    } else { # !shared_hp_tasks && !shared_hp_outputs
      # MODIFICATION 5: Générer 'p_t' unique pour chaque combinaison
      hps_unique <- base_ids %>%
        dplyr::left_join(hp_config %>%
                           dplyr::mutate(Output_ID = as.character(output_id)),
                         by = "Output_ID") %>%
        dplyr::mutate(
          p_t = purrr::map2_dbl(pt_min, pt_max, ~stats::runif(1, .x, .y)), # Changé de l_t, lt_min, lt_max
          S_t = purrr::map2_dbl(St_min, St_max, ~stats::runif(1, .x, .y))
        ) %>%
        dplyr::select(Task_ID, Output_ID, p_t, S_t) # Changé de l_t
      if (noise) {
        # La partie bruit reste inchangée mais on la reconstruit pour la clarté
        hps_noise <- base_ids %>%
          dplyr::left_join(hp_config %>%
                             dplyr::mutate(Output_ID = as.character(output_id)),
                           by = "Output_ID") %>%
          dplyr::mutate(
            noise = purrr::map2_dbl(noise_min, noise_max, ~stats::runif(1, .x, .y))
          ) %>%
          dplyr::select(Task_ID, Output_ID, noise)
        final_hp <- dplyr::left_join(hps_unique, hps_noise, by = c("Task_ID", "Output_ID"))
      } else {
        final_hp <- hps_unique
      }
    }

    # MODIFICATION 6: Gestion de 'p_u_t' au lieu de 'l_u_t'
    if (shared_hp_tasks) {
      final_hp$p_u_t <- stats::runif(1, hp_config$pu_min[1], hp_config$pu_max[1]) # Changé de l_u_t, lu_min, lu_max
    } else {
      p_u_t_per_task <- tibble::tibble(
        Task_ID = as.character(list_task_ID),
        p_u_t   = stats::runif(num_tasks, hp_config$pu_min[1], hp_config$pu_max[1]) # Changé de l_u_t, lu_min, lu_max
      )
      # On doit s'assurer de ne pas écraser une colonne existante si elle a déjà été créée
      if ("p_u_t" %in% names(final_hp)) {
        final_hp <- final_hp %>% dplyr::select(-p_u_t)
      }
      final_hp <- dplyr::left_join(final_hp, p_u_t_per_task, by = "Task_ID")
    }

    return(final_hp)

  } else {
    # Case 2: Kernel is provided as a string (e.g., "SE", "SE + PERIO")
    #     # This logic remains as it was in your original function.
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

# hp <- function(kern = "SE",
#                list_task_ID = NULL,
#                list_output_ID = NULL,
#                shared_hp_tasks = TRUE,
#                shared_hp_outputs = TRUE,
#                noise = FALSE,
#                hp_config = NULL) {
#   ## Initiate interval boundaries
#   min_val <- -3
#   max_val <- 3
#   # min_val_u <- -1
#   # max_val_u <- 3
#   min_noise <- -5
#   max_noise <- -1
#
#   # browser()
#   # Convolution case
#   if (is.function(kern)) {
#     kern_name <- deparse(substitute(kern))
#     if (kern_name != "convolution_kernel") {
#       stop("Currently, only 'convolution_kernel' is supported as a function input.")
#     }
#     if (is.null(list_task_ID) || is.null(list_output_ID)) {
#       stop("For the convolution_kernel, both 'list_task_ID' and 'list_output_ID' must be provided.")
#     }
#
#     if (is.null(hp_config)) {
#       message("hp_config not provided for convolution_kernel, using default HP bounds.")
#       hp_config <- tibble::tibble(
#         output_id   = list_output_ID,
#         lt_min      = -2, lt_max      = 2,
#         St_min      = -2, St_max      = 2,
#         lu_min      = -2, lu_max      = 0,
#         noise_min   = -5, noise_max   = -2
#       )
#     }
#
#     num_tasks <- length(list_task_ID)
#     num_outputs <- length(list_output_ID)
#
#     # All combination task x output
#     base_ids <- tidyr::crossing(Task_ID = as.character(list_task_ID),
#                                 Output_ID = as.character(list_output_ID))
#
#     if (shared_hp_tasks && shared_hp_outputs) {
#       hps_to_add <- tibble::tibble(
#         l_t = stats::runif(1, hp_config$lt_min[1], hp_config$lt_max[1]),
#         S_t = stats::runif(1, hp_config$St_min[1], hp_config$St_max[1])
#       )
#       if (noise) {
#         hps_to_add$noise <- stats::runif(1, hp_config$noise_min[1], hp_config$noise_max[1])
#       }
#       final_hp <- tidyr::crossing(base_ids, hps_to_add)
#
#     } else if (shared_hp_tasks && !shared_hp_outputs) {
#       hps_per_output <- hp_config %>%
#         dplyr::transmute(
#           Output_ID = as.character(output_id),
#           l_t = purrr::map2_dbl(lt_min, lt_max, ~stats::runif(1, .x, .y)),
#           S_t = purrr::map2_dbl(St_min, St_max, ~stats::runif(1, .x, .y))
#         )
#       if (noise) {
#         hps_per_output$noise <- purrr::map2_dbl(hp_config$noise_min,
#                                                 hp_config$noise_max,
#                                                 ~stats::runif(1, .x, .y))
#       }
#       final_hp <- dplyr::left_join(base_ids, hps_per_output, by = "Output_ID")
#
#     } else if (!shared_hp_tasks && shared_hp_outputs) {
#       hps_per_task <- tibble::tibble(
#         Task_ID = as.character(list_task_ID),
#         l_t = stats::runif(num_tasks, hp_config$lt_min[1], hp_config$lt_max[1]),
#         S_t = stats::runif(num_tasks, hp_config$St_min[1], hp_config$St_max[1])
#       )
#       if (noise) {
#         hps_per_task$noise <- stats::runif(num_tasks,
#                                            hp_config$noise_min[1],
#                                            hp_config$noise_max[1])
#       }
#       final_hp <- dplyr::left_join(base_ids, hps_per_task, by = "Task_ID")
#
#     } else { # !shared_hp_tasks && !shared_hp_outputs
#       hps_unique <- base_ids %>%
#         dplyr::left_join(hp_config %>%
#                            dplyr::mutate(Output_ID = as.character(output_id)),
#                          by = "Output_ID") %>%
#         dplyr::mutate(
#           l_t = purrr::map2_dbl(lt_min, lt_max, ~stats::runif(1, .x, .y)),
#           S_t = purrr::map2_dbl(St_min, St_max, ~stats::runif(1, .x, .y))
#         ) %>%
#         dplyr::select(Task_ID, Output_ID, l_t, S_t)
#       if (noise) {
#         hps_unique <- base_ids %>%
#           dplyr::left_join(hp_config %>%
#                              dplyr::mutate(Output_ID = as.character(output_id)),
#                            by = "Output_ID") %>%
#           dplyr::mutate(
#             noise = purrr::map2_dbl(noise_min, noise_max, ~stats::runif(1, .x, .y))
#           ) %>%
#           dplyr::select(Task_ID, Output_ID, noise) %>%
#           dplyr::left_join(hps_unique, ., by = c("Task_ID", "Output_ID"))
#       }
#       final_hp <- hps_unique
#     }
#
#     # --- Gestion de l_u_t (ne dépend que des tâches) ---
#     if (shared_hp_tasks) {
#       final_hp$l_u_t <- stats::runif(1, hp_config$lu_min[1], hp_config$lu_max[1])
#     } else {
#       l_u_t_per_task <- tibble::tibble(
#         Task_ID = as.character(list_task_ID),
#         l_u_t   = stats::runif(num_tasks, hp_config$lu_min[1], hp_config$lu_max[1])
#       )
#       final_hp <- dplyr::left_join(final_hp, l_u_t_per_task, by = "Task_ID")
#     }
#
#     return(final_hp)
#
#   } else {
#     # Case 2: Kernel is provided as a string (e.g., "SE", "SE + PERIO")
#     # This logic remains as it was in your original function.
#     base_ids <- NULL
#     n_draws <- 1
#
#     if (!is.null(list_task_ID) && !is.null(list_output_ID)) {
#       base_ids <- tidyr::crossing(Task_ID = as.character(list_task_ID),
#                                   Output_ID = as.character(list_output_ID))
#       if (!shared_hp_tasks && !shared_hp_outputs) {
#         n_draws <- nrow(base_ids)
#       } else if (!shared_hp_tasks && shared_hp_outputs) {
#         n_draws <- length(list_task_ID)
#       } else if (shared_hp_tasks && !shared_hp_outputs) {
#         n_draws <- length(list_output_ID)
#       } # else n_draws remains 1
#     } else if (!is.null(list_output_ID)) {
#       base_ids <- tibble::tibble(Output_ID = as.character(list_output_ID))
#       n_draws <- if (shared_hp_outputs) 1 else nrow(base_ids)
#     } else if (!is.null(list_task_ID)) {
#       base_ids <- tibble::tibble(Task_ID = as.character(list_task_ID))
#       n_draws <- if (shared_hp_tasks) 1 else nrow(base_ids)
#     }
#
#     str_kern <- strsplit(kern, " +")[[1]]
#
#     # Generate the required number of unique HP sets
#     generated_hps <- tibble::tibble(.rows = n_draws)
#     for (i in str_kern) {
#       temp_hp <- switch(i,
#                         "SE" = tibble::tibble(
#                           se_variance = stats::runif(n_draws, min_val, max_val),
#                           se_lengthscale = stats::runif(n_draws, min_val, max_val)
#                         ),
#                         "PERIO" = tibble::tibble(
#                           perio_variance = stats::runif(n_draws, min_val, max_val),
#                           perio_lengthscale = stats::runif(n_draws, min_val, max_val),
#                           period = stats::runif(n_draws, 0, 2 * pi)
#                         ),
#                         "RQ" = tibble::tibble(
#                           rq_variance = stats::runif(n_draws, min_val, max_val),
#                           rq_lengthscale = stats::runif(n_draws, min_val, max_val),
#                           rq_scale = stats::runif(n_draws, min_val, max_val)
#                         ),
#                         "LIN" = tibble::tibble(
#                           lin_slope = stats::runif(n_draws, min_val, max_val),
#                           lin_offset = stats::runif(n_draws, min_val, max_val)
#                         ),
#                         # Default case for operators like '+' or '*'
#                         tibble::tibble()
#       )
#       generated_hps <- dplyr::bind_cols(generated_hps, temp_hp)
#     }
#
#     if (noise) {
#       generated_hps <- generated_hps %>%
#         dplyr::mutate(noise = stats::runif(n_draws, min_noise, max_noise))
#     }
#
#     # Construct the Final Tibble by combining IDs and HPs
#     if (is.null(base_ids)) {
#       return(generated_hps)
#     }
#
#     if (!is.null(list_task_ID) && !is.null(list_output_ID)) {
#       if (!shared_hp_tasks && !shared_hp_outputs) {
#         final_hp <- dplyr::bind_cols(base_ids, generated_hps)
#       } else if (!shared_hp_tasks && shared_hp_outputs) {
#         task_hps <- dplyr::bind_cols(Task_ID = as.character(list_task_ID), generated_hps)
#         final_hp <- dplyr::left_join(base_ids, task_hps, by = "Task_ID")
#       } else if (shared_hp_tasks && !shared_hp_outputs) {
#         output_hps <- dplyr::bind_cols(Output_ID = as.character(list_output_ID), generated_hps)
#         final_hp <- dplyr::left_join(base_ids, output_hps, by = "Output_ID")
#       } else {
#         final_hp <- tidyr::crossing(base_ids, generated_hps)
#       }
#     } else { # Case for only one list provided (task or output)
#       if (n_draws == 1) { # Shared HPs
#         final_hp <- tidyr::crossing(base_ids, generated_hps)
#       } else { # Independent HPs
#         final_hp <- dplyr::bind_cols(base_ids, generated_hps)
#       }
#     }
#
#     return(final_hp)
#   }
# }


#' @title Generate Initial Hyperparameters for Tasks and Outputs
#' @description Creates an initial `hp_t` tibble based on prior distributions,
#'   handling both shared and specific hyperparameters.
#'
#' @param hp_t_config A named list where each element defines a prior.
#'   Names can be generic (`l_t`) or specific (`l_t_T1`, `l_t_Oa`, `l_t_T1_Oa`)
#'   following the format `baseName_taskID_outputID`.
#' @param task_ids A vector of unique task IDs.
#' @param output_ids A vector of unique output IDs.
#' @param base_hp_names A character vector of the base hyperparameter names
#'   (e.g., `c("l_u_t", "l_t", "S_t", "noise")`).
#' @param shared_tasks A boolean. If `TRUE`, HPs are identical across all tasks
#'   (unless a specific output prior is provided when `shared_hp_outputs = FALSE`).
#' @param shared_hp_outputs A boolean. If `TRUE`, HPs are identical across all outputs.
#'
#' @return A tibble `hp_t` with initial values for each combination
#'   of `Task_ID` and `Output_ID`.

generate_initial_hp_t <- function(hp_t_config,
                                  task_ids,
                                  output_ids,
                                  base_hp_names,
                                  shared_tasks = TRUE,
                                  shared_hp_outputs = TRUE) {

  # Create the base tibble with all combinations
  hp_t <- tidyr::crossing(Task_ID = task_ids, Output_ID = output_ids)

  # --- Case 1: HPs shared across all tasks AND all outputs ---
  if (shared_tasks && shared_hp_outputs) {
    initial_values <- list()
    for (name in base_hp_names) {
      prior <- hp_t_config[[name]]$prior
      if (is.null(prior)) stop(paste("Error: Generic prior is missing for HP:", name))

      value <- switch(prior$dist,
                      "gamma"    = prior$shape / prior$rate,
                      "invgamma" = if (prior$shape > 1) prior$scale / (prior$shape - 1) else NA,
                      "normal"   = prior$mean,
                      stop(paste("Unsupported prior distribution:", prior$dist))
      )
      initial_values[[name]] <- value
    }
    hp_t <- dplyr::bind_cols(hp_t, tibble::as_tibble(initial_values))
  }

  # --- Case 2: HPs specific to each task, but shared across outputs ---
  else if (!shared_tasks && shared_hp_outputs) {
    for (hp_name in base_hp_names) {
      values_for_tasks <- sapply(task_ids, function(current_task) {
        specific_name <- paste0(hp_name, "_", current_task)

        prior_name <- if (specific_name %in% names(hp_t_config)) specific_name else hp_name
        prior <- hp_t_config[[prior_name]]$prior
        if (is.null(prior)) stop(paste("Error: No prior for HP:", hp_name, "for task:", current_task))

        switch(prior$dist,
               "gamma"    = prior$shape / prior$rate,
               "invgamma" = if (prior$shape > 1) prior$scale / (prior$shape - 1) else NA,
               "normal"   = prior$mean,
               stop(paste("Unsupported prior:", prior$dist)))
      })
      values_tibble <- tibble(Task_ID = task_ids, !!hp_name := values_for_tasks)
      hp_t <- dplyr::left_join(hp_t, values_tibble, by = "Task_ID")
    }
  }

  # --- Case 3: HPs shared across tasks, but specific to each output ---
  else if (shared_tasks && !shared_hp_outputs) {
    for (hp_name in base_hp_names) {
      values_for_outputs <- sapply(output_ids, function(current_output) {
        specific_name <- paste0(hp_name, "_", current_output)

        prior_name <- if (specific_name %in% names(hp_t_config)) specific_name else hp_name
        prior <- hp_t_config[[prior_name]]$prior
        if (is.null(prior)) stop(paste("Error: No prior for HP:", hp_name, "for output:", current_output))

        switch(prior$dist,
               "gamma"    = prior$shape / prior$rate,
               "invgamma" = if (prior$shape > 1) prior$scale / (prior$shape - 1) else NA,
               "normal"   = prior$mean,
               stop(paste("Unsupported prior:", prior$dist)))
      })
      values_tibble <- tibble(Output_ID = output_ids, !!hp_name := values_for_outputs)
      hp_t <- dplyr::left_join(hp_t, values_tibble, by = "Output_ID")
    }
  }

  # --- Case 4: HPs specific to each combination of task AND output ---
  else if (!shared_tasks && !shared_hp_outputs) {
    # Iterate over each hyperparameter name
    for (hp_name in base_hp_names) {
      # Use purrr::map2_dbl to apply a function over each pair of Task_ID and Output_ID
      hp_values <- purrr::map2_dbl(hp_t$Task_ID, hp_t$Output_ID, function(current_task, current_output) {

        # Hierarchical search logic for the prior
        name_task_output <- paste0(hp_name, "_", current_task, "_", current_output)
        name_task        <- paste0(hp_name, "_", current_task)
        name_output      <- paste0(hp_name, "_", current_output)

        prior_name <- if (name_task_output %in% names(hp_t_config)) {
          name_task_output
        } else if (name_task %in% names(hp_t_config)) {
          name_task
        } else if (name_output %in% names(hp_t_config)) {
          name_output
        } else {
          hp_name
        }

        prior <- hp_t_config[[prior_name]]$prior
        if (is.null(prior)) stop(paste("Error: No prior for HP:", hp_name, "for combo:", current_task, current_output))

        switch(prior$dist,
               "gamma"    = prior$shape / prior$rate,
               "invgamma" = if (prior$shape > 1) prior$scale / (prior$shape - 1) else NA,
               "normal"   = prior$mean,
               stop(paste("Unsupported prior:", prior$dist)))
      })

      # Assign the computed values as a new column in the tibble
      hp_t[[hp_name]] <- hp_values
    }
  }

  return(hp_t)
}
