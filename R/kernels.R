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
#'   Names embed the Output_ID label: 'l_t_<o>', 'S_t_<o>', 'l_u_t' (single
#' @return The covariance matrix or its partial derivative.
#'
#' @export
convolution_kernel <- function(x,
                                  y,
                                  hp,
                                  vectorized = FALSE,
                                  deriv = NULL) {

  # Identify the input dimension from the numeric coordinate columns.
  if (is.data.frame(x)) {
    coord_cols <- identify_coord_cols(x)
    input_dim <- length(coord_cols)
  } else {
    input_dim <- length(x)
  }

  ## Non-vectorized single-pair evaluation (kept for completeness).
  if (!vectorized) {
    l_i <- exp(hp[["lengthscale_output1"]])
    l_j <- exp(hp[["lengthscale_output2"]])
    l_u <- exp(hp[["lengthscale_u"]])
    S_i <- exp(hp[["variance_output1"]])
    S_j <- exp(hp[["variance_output2"]])
    distance_sq <- sum((x - y)^2)
    denominator_sum <- l_i + l_j + l_u
    K <- (S_i * S_j) / ((2 * pi * denominator_sum)^(input_dim / 2)) *
      exp(-0.5 * distance_sq / denominator_sum)
    if (!is.null(deriv)) {
      stop("Derivatives of the convolution kernel are only implemented in ",
           "vectorized mode.")
    }
    return(K)
  }

  ## Vectorized evaluation: build the full multi-output covariance as a SUM
  ## over Q latent processes (Alvarez & Lawrence convolution processes). A
  ## single latent (no 'Latent_ID' column) reduces to the original kernel.
  if (!"Output_ID" %in% names(x) || !"Output_ID" %in% names(hp)) {
    stop("'input' and 'hp' must contain an 'Output_ID' column for ",
         "vectorized mode.")
  }

  x_ids <- as.character(x$Output_ID)
  y_ids <- as.character(y$Output_ID)
  x_coords <- as.matrix(x[, coord_cols, drop = FALSE])
  y_coords <- as.matrix(y[, coord_cols, drop = FALSE])
  distance_sq <- cpp_dist(x_coords, y_coords)
  N <- nrow(x)
  M <- nrow(y)

  if (!"Latent_ID" %in% names(hp)) hp$Latent_ID <- "1"
  hp_lat <- as.character(hp$Latent_ID)
  lat_levels <- unique(hp_lat)
  multi_latent <- length(lat_levels) > 1

  K <- matrix(0, nrow = N, ncol = M)
  denom_list <- K_list <- l1_list <- l2_list <-
    stats::setNames(vector("list", length(lat_levels)), lat_levels)

  for (q in lat_levels) {
    rows_q <- hp_lat == q
    out_q <- as.character(hp$Output_ID[rows_q])
    ix <- match(x_ids, out_q)
    iy <- match(y_ids, out_q)
    l1 <- exp(hp$l_t[rows_q][ix])
    l2 <- exp(hp$l_t[rows_q][iy])
    S1 <- exp(hp$S_t[rows_q][ix])
    S2 <- exp(hp$S_t[rows_q][iy])
    l_u_q <- exp(unique(hp$l_u_t[rows_q]))

    denom_q <- outer(l1, l2, FUN = "+") + l_u_q
    S_prod <- outer(S1, S2, FUN = "*")
    K_q <- S_prod / ((2 * pi * denom_q)^(input_dim / 2)) *
      exp(-0.5 * distance_sq / denom_q)

    K <- K + K_q
    denom_list[[q]] <- denom_q; K_list[[q]] <- K_q
    l1_list[[q]] <- l1; l2_list[[q]] <- l2
  }

  if (is.null(deriv)) {
    return(K)
  }

  ## Derivative names embed the Output_ID / Latent_ID labels directly:
  ##   single latent : 'l_t_{o}', 'S_t_{o}', 'l_u_t'
  ##   multi-latent  : 'l_t_o{o}_l{q}', 'S_t_o{o}_l{q}', 'l_u_t_l{q}'
  ## (latent labels are internal integers). Labels are matched, not decoded.
  if (deriv == "l_u_t" || grepl("^l_u_t_l", deriv)) {
    lat_lab <- if (deriv == "l_u_t") lat_levels[1] else sub("^l_u_t_l", "", deriv)
    denom_q <- denom_list[[lat_lab]]; K_q <- K_list[[lat_lab]]
    common <- K_q * ((-0.5 * input_dim / denom_q) +
                       (0.5 * distance_sq / (denom_q^2)))
    return(common * exp(unique(hp$l_u_t[hp_lat == lat_lab])))

  } else if (startsWith(deriv, "l_t")) {
    if (multi_latent) {
      lat_lab <- sub("^.*_l", "", deriv)
      out_lab <- sub("^l_t_o(.*)_l[0-9]+$", "\\1", deriv)
    } else {
      out_lab <- sub("^l_t_", "", deriv); lat_lab <- lat_levels[1]
    }
    denom_q <- denom_list[[lat_lab]]; K_q <- K_list[[lat_lab]]
    common <- K_q * ((-0.5 * input_dim / denom_q) +
                       (0.5 * distance_sq / (denom_q^2)))
    l_i_mat <- matrix(l1_list[[lat_lab]], nrow = N, ncol = M)
    l_j_mat <- matrix(l2_list[[lat_lab]], nrow = N, ncol = M, byrow = TRUE)
    mask_i <- outer(x_ids == out_lab, rep(TRUE, M))
    mask_j <- outer(rep(TRUE, N), y_ids == out_lab)
    return(common * (l_i_mat * mask_i + l_j_mat * mask_j))

  } else if (startsWith(deriv, "S_t")) {
    if (multi_latent) {
      lat_lab <- sub("^.*_l", "", deriv)
      out_lab <- sub("^S_t_o(.*)_l[0-9]+$", "\\1", deriv)
    } else {
      out_lab <- sub("^S_t_", "", deriv); lat_lab <- lat_levels[1]
    }
    K_q <- K_list[[lat_lab]]
    mask_i <- outer(x_ids == out_lab, rep(TRUE, M))
    mask_j <- outer(rep(TRUE, N), y_ids == out_lab)
    return(K_q * (mask_i + mask_j))

  } else {
    stop("Invalid 'deriv' argument: ", deriv)
  }
}


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
  vectorized = FALSE
) {
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
    stop(
      "Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'"
    )
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
  vectorized = FALSE
) {
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
    (exp(hp[["perio_variance"]]) *
      exp(-2 * exp(-hp[["perio_lengthscale"]]) * perio_term) *
      2 *
      exp(-hp[["perio_lengthscale"]]) *
      sum_deriv) %>%
      return()
  } else if (deriv == "perio_lengthscale") {
    (exp(hp[["perio_variance"]]) *
      2 *
      perio_term *
      exp(-hp[["perio_lengthscale"]]) *
      exp(-2 * perio_term * exp(-hp[["perio_lengthscale"]]))) %>%
      return()
  } else {
    stop(
      "Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'"
    )
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
  vectorized = FALSE
) {
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
    (exp(hp[["rq_variance"]]) *
      term^(-hp[["rq_scale"]]) *
      (distance *
        exp(-hp[["rq_lengthscale"]]) /
        (2 * hp[["rq_scale"]] * term) -
        log(term))) %>%
      return()
  } else if (deriv == "rq_lengthscale") {
    (exp(hp[["rq_variance"]]) *
      distance *
      0.5 *
      exp(-hp[["rq_lengthscale"]]) *
      term^(-hp[["rq_scale"]] - 1)) %>%
      return()
  } else {
    stop(
      "Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'"
    )
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
  vectorized = FALSE
) {
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
    stop(
      "Please enter a valid hyper-parameter's name or NULL for the
    argument 'deriv'"
    )
  }
}


#' Generate initial hyper-parameters for GP and multi-output kernels
#'
#' `hp()` draws a tibble of random initial hyper-parameters with the correct
#' structure for a given kernel. It supports two regimes:
#' * Single-output kernels given as a character string (e.g. "SE",
#'   "SE + PERIO"): returns one row of the kernel's named hyper-parameters.
#' * The multi-output convolution kernel given as a function
#'   ([convolution_kernel()]): returns a tibble keyed by `Task_ID` and
#'   `Output_ID`, with the smoothing-kernel parameters `l_t`, `S_t`, the
#'   shared latent-process lengthscale `l_u_t`, and optionally `noise`.
#'
#' @param kern A character string naming a single-output kernel (one of
#'   "SE", "LIN", "PERIO", "RQ", or a whitespace-separated combination
#'   such as "SE + PERIO"), or a kernel *function* for the multi-output case
#'   (typically [convolution_kernel()]).
#' @param list_task_ID A vector of task identifiers. Required for the
#'   multi-output convolution kernel.
#' @param list_output_ID A vector of output identifiers (labels; may be arbitrary strings, e.g. 'weight', 'height').
#'   Required for the multi-output convolution kernel.
#' @param shared_hp_tasks A logical value. If TRUE (default), hyper-parameters
#'   are shared across tasks; otherwise they are drawn independently per task.
#' @param shared_hp_outputs A logical value. If TRUE, the smoothing-kernel
#'   hyper-parameters are shared across outputs; otherwise (default) they are
#'   output-specific.
#' @param n_latent An integer. The number of latent processes Q in the
#'   convolution kernel (default 1). When greater than 1, hyper-parameters are
#'   drawn per latent process and a 'Latent_ID' column is added.
#' @param noise A logical value. If TRUE, a 'noise' hyper-parameter is added.
#' @param hp_config An optional tibble of per-output generation bounds for the
#'   convolution kernel, with an 'Output_ID' column and columns
#'   'lt_min'/'lt_max', 'St_min'/'St_max', 'lu_min'/'lu_max',
#'   'noise_min'/'noise_max'. If NULL (default), sensible bounds are used.
#'
#' @return A `tibble` of initial hyper-parameters.
#' @export
#'
#' @examples
#' ## Single-output kernel
#' hp("SE", noise = TRUE)
#'
#' ## Multi-output convolution kernel: 2 tasks, 3 outputs
#' hp(convolution_kernel, list_task_ID = c("a", "b"),
#'    list_output_ID = 1:3, noise = TRUE)

hp <- function(kern = "SE",
               list_task_ID = NULL,
               list_output_ID = NULL,
               shared_hp_tasks = TRUE,
               shared_hp_outputs = FALSE,
               n_latent = 1,
               noise = FALSE,
               hp_config = NULL) {
  ## Initiate interval boundaries
  min_val_l <- -3
  max_val_l <- 0
  min_val_S <- log(0.5)
  max_val_S <- log(10)
  min_noise <- -5
  max_noise <- -2
  min_val <- -3
  max_val <- 10

  ## ---- Multi-output convolution kernel -----------------------------------
  if (is.function(kern)) {
    if (is.null(list_task_ID) || is.null(list_output_ID)) {
      stop("For a multi-output convolution kernel, both 'list_task_ID' and ",
           "'list_output_ID' must be provided.", call. = FALSE)
    }

    ## Output/Task IDs are identifier labels used throughout; they may be
    ## arbitrary strings (e.g. 'weight', 'height').
    out_ids  <- sort(unique(as.character(list_output_ID)))
    task_ids <- unique(as.character(list_task_ID))
    task_fct <- task_ids
    out_fct  <- out_ids

    ## Per-output generation bounds (log-scale); overridable via 'hp_config'.
    config <- .hp_conv_config(hp_config, out_ids)

    ## ---- Multiple latent processes -------------------------------------
    ## Q > 1 latent convolution processes: draw l_t/S_t per (output, latent),
    ## l_u_t per latent, and a per-output noise shared across latents.
    if (n_latent > 1) {
      grid <- tidyr::crossing(Task_ID = task_fct, Output_ID = out_fct,
                              Latent_ID = as.character(seq_len(n_latent)))
      gi_out  <- match(as.character(grid$Output_ID), out_ids)
      gi_task <- match(as.character(grid$Task_ID), task_ids)
      gi_lat  <- as.character(grid$Latent_ID)
      draw_group <- function(param, use_output, use_latent) {
        lo <- config[[paste0(param, "_min")]]
        hi <- config[[paste0(param, "_max")]]
        key <- paste(
          if (!shared_hp_tasks) gi_task else 0L,
          if (use_output && !shared_hp_outputs) gi_out else 0L,
          if (use_latent) gi_lat else 0L,
          sep = "_"
        )
        val <- numeric(length(key))
        for (k in unique(key)) {
          idx <- which(key == k)
          oi <- if (use_output && !shared_hp_outputs) gi_out[idx[1]] else 1L
          val[idx] <- stats::runif(1, lo[oi], hi[oi])
        }
        val
      }
      final_hp <- grid
      final_hp$l_t   <- draw_group("lt", TRUE,  TRUE)
      final_hp$S_t   <- draw_group("St", TRUE,  TRUE)
      final_hp$l_u_t <- draw_group("lu", FALSE, TRUE)
      if (noise) final_hp$noise <- draw_group("noise", TRUE, FALSE)
      return(final_hp)
    }

    ## Full Task x Output grid onto which every hyper-parameter is broadcast.
    grid <- tidyr::crossing(Task_ID = task_fct, Output_ID = out_fct)

    ## Draw a smoothing-kernel parameter ('lt', 'St' or 'noise') at the
    ## granularity implied by the sharing flags, then broadcast onto 'grid'.
    draw_smoothing <- function(param) {
      lo <- config[[paste0(param, "_min")]]
      hi <- config[[paste0(param, "_max")]]
      if (shared_hp_outputs) {
        if (shared_hp_tasks) {
          rep(stats::runif(1, lo[1], hi[1]), nrow(grid))
        } else {
          vals <- stats::runif(length(task_ids), lo[1], hi[1])
          vals[match(grid$Task_ID, task_fct)]
        }
      } else {
        if (shared_hp_tasks) {
          vals <- stats::runif(length(out_ids), lo, hi)
          vals[match(grid$Output_ID, out_fct)]
        } else {
          idx <- match(grid$Output_ID, out_fct)
          stats::runif(nrow(grid), lo[idx], hi[idx])
        }
      }
    }

    final_hp <- grid %>%
      dplyr::mutate(l_t = draw_smoothing("lt"), S_t = draw_smoothing("St"))
    if (noise) final_hp$noise <- draw_smoothing("noise")

    ## A single latent process is shared by all outputs, so 'l_u_t' never
    ## varies by output. It varies by task only when HPs are task-specific.
    if (shared_hp_tasks) {
      final_hp$l_u_t <- stats::runif(1, config$lu_min[1], config$lu_max[1])
    } else {
      lu <- stats::runif(length(task_ids), config$lu_min[1], config$lu_max[1])
      final_hp$l_u_t <- lu[match(final_hp$Task_ID, task_fct)]
    }

    return(final_hp)

  } else {
    # Case 2: Kernel is provided as a string (e.g., "SE", "SE + PERIO")
    base_ids <- NULL
    n_draws <- 1

    if (!is.null(list_task_ID) && !is.null(list_output_ID)) {
      base_ids <- tidyr::crossing(Task_ID = as.factor(list_task_ID),
                                  Output_ID = as.factor(list_output_ID))
      if (!shared_hp_tasks && !shared_hp_outputs) {
        n_draws <- nrow(base_ids)
      } else if (!shared_hp_tasks && shared_hp_outputs) {
        n_draws <- length(list_task_ID)
      } else if (shared_hp_tasks && !shared_hp_outputs) {
        n_draws <- length(list_output_ID)
      } # else n_draws remains 1
    } else if (!is.null(list_output_ID)) {
      base_ids <- tibble::tibble(Output_ID = as.factor(list_output_ID))
      n_draws <- if (shared_hp_outputs) 1 else nrow(base_ids)
    } else if (!is.null(list_task_ID)) {
      base_ids <- tibble::tibble(Task_ID = as.factor(list_task_ID))
      n_draws <- if (shared_hp_tasks) 1 else nrow(base_ids)
    }

    ## 'kern = NULL' is used to request a HP tibble with no base-kernel
    ## parameters (e.g. to append only a 'noise' column to an existing tibble).
    str_kern <- if (is.null(kern)) character(0) else strsplit(kern, " +")[[1]]

    # Generate the required number of unique HP sets
    generated_hps <- tibble::tibble(.rows = n_draws)
    for (i in str_kern) {
      temp_hp <- switch(i,
                        "SE" = tibble::tibble(
                          se_variance = stats::runif(n_draws, min_val_S, max_val_S),
                          se_lengthscale = stats::runif(n_draws, min_val_l, max_val_l)
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
                        NULL
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
        task_hps <- dplyr::bind_cols(Task_ID = as.factor(list_task_ID),
                                     generated_hps)
        final_hp <- dplyr::left_join(base_ids, task_hps, by = "Task_ID")
      } else if (shared_hp_tasks && !shared_hp_outputs) {
        output_hps <- dplyr::bind_cols(Output_ID = as.factor(list_output_ID),
                                       generated_hps)
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

#' Build the per-output generation bounds for the convolution kernel
#'
#' @param hp_config A user-supplied tibble of bounds, or NULL for defaults.
#' @param out_ids An integer vector of (sorted, unique) output IDs.
#' @return A tibble of bounds with one row per output, in `out_ids` order.
#' @keywords internal
.hp_conv_config <- function(hp_config, out_ids) {
  default <- tibble::tibble(
    Output_ID = out_ids,
    lt_min = log(1 / 1000), lt_max = log(1 / 100),
    St_min = log(0.4),      St_max = log(20),
    lu_min = log(1 / 1000), lu_max = log(1 / 100),
    noise_min = -5,         noise_max = -2
  )
  if (is.null(hp_config)) return(default)
  ## Accept the historical lower-case 'output_id' spelling.
  if ("output_id" %in% names(hp_config) && !("Output_ID" %in% names(hp_config))) {
    hp_config <- dplyr::rename(hp_config, Output_ID = "output_id")
  }
  needed <- setdiff(names(default), "Output_ID")
  miss <- setdiff(needed, names(hp_config))
  if (length(miss) > 0) {
    stop("'hp_config' is missing column(s): ", paste(miss, collapse = ", "),
         ".", call. = FALSE)
  }
  hp_config %>%
    dplyr::mutate(Output_ID = as.character(.data$Output_ID)) %>%
    dplyr::arrange(match(.data$Output_ID, out_ids))
}

#' Flatten a multi-output hyper-parameter tibble into a named vector
#'
#' Converts the canonical multi-output hyper-parameter tibble into the flat
#' named numeric vector consumed by `stats::optim()`. Names are the derivative
#' identifiers understood by `convolution_kernel()`. Names embed the Output_ID
#' label: 'l_t_<o>', 'S_t_<o>', 'l_u_t', 'noise_<o>' (single latent); with a
#' 'Latent_ID' column: 'l_t_o<o>_l<q>', 'S_t_o<o>_l<q>', 'l_u_t_l<q>',
#' 'noise_o<o>'.
#'
#' @param hp A hyper-parameter tibble (single individual, no 'Task_ID').
#' @return A named numeric vector.
#' @keywords internal
flatten_hp_mo <- function(hp) {
  if ("Latent_ID" %in% names(hp) &&
      length(unique(as.character(hp$Latent_ID))) > 1) {
    hp <- hp %>% dplyr::arrange(as.character(.data$Latent_ID),
                                as.character(.data$Output_ID))
    ol <- as.character(hp$Output_ID); ql <- as.character(hp$Latent_ID)
    vec <- c(
      stats::setNames(hp$l_t, paste0("l_t_o", ol, "_l", ql)),
      stats::setNames(hp$S_t, paste0("S_t_o", ol, "_l", ql))
    )
    lu <- hp %>% dplyr::distinct(.data$Latent_ID, .data$l_u_t) %>%
      dplyr::arrange(as.character(.data$Latent_ID))
    if (anyDuplicated(as.character(lu$Latent_ID))) {
      stop("flatten_hp_mo(): 'l_u_t' must be constant within each latent.",
           call. = FALSE)
    }
    vec <- c(vec, stats::setNames(lu$l_u_t,
      paste0("l_u_t_l", as.character(lu$Latent_ID))))
    if ("noise" %in% names(hp)) {
      ns <- hp %>% dplyr::distinct(.data$Output_ID, .data$noise) %>%
        dplyr::arrange(as.character(.data$Output_ID))
      vec <- c(vec, stats::setNames(ns$noise,
        paste0("noise_o", as.character(ns$Output_ID))))
    }
    return(vec)
  }
  hp <- hp %>% dplyr::arrange(as.character(.data$Output_ID))
  lab <- as.character(hp$Output_ID)
  lu <- unique(hp$l_u_t)
  if (length(lu) != 1) {
    stop("flatten_hp_mo() expects a single shared 'l_u_t' value.",
         call. = FALSE)
  }
  vec <- c(
    stats::setNames(hp$l_t, paste0("l_t_", lab)),
    stats::setNames(hp$S_t, paste0("S_t_", lab)),
    stats::setNames(lu, "l_u_t")
  )
  if ("noise" %in% names(hp)) {
    vec <- c(vec, stats::setNames(hp$noise, paste0("noise_", lab)))
  }
  vec
}

#' Rebuild a multi-output hyper-parameter tibble from a flat named vector
#'
#' Inverse of [flatten_hp_mo()]. Multi-latent names (containing a latent
#' suffix '_l<q>') reconstruct a 'Latent_ID' column.
#'
#' @param hp A named numeric vector as produced by [flatten_hp_mo()].
#' @return A hyper-parameter tibble.
#' @keywords internal
unflatten_hp_mo <- function(hp) {
  nm <- names(hp)
  if (any(grepl("^l_u_t_l", nm))) {
    lt <- hp[grepl("^l_t_o", nm)]
    st <- hp[grepl("^S_t_o", nm)]
    lu <- hp[grepl("^l_u_t_l", nm)]
    ns <- hp[grepl("^noise_o", nm)]
    d <- sub("^l_t_o(.*)_l[0-9]+$", "\\1", names(lt))
    q <- sub("^.*_l", "", names(lt))
    lat <- sub("^l_u_t_l", "", names(lu))
    lu_map <- stats::setNames(as.numeric(lu), lat)
    out <- tibble::tibble(
      Output_ID = d,
      Latent_ID = q,
      l_t = as.numeric(lt),
      S_t = as.numeric(st[paste0("S_t_o", d, "_l", q)]),
      l_u_t = as.numeric(lu_map[q])
    )
    if (length(ns) > 0) {
      nd <- sub("^noise_o", "", names(ns))
      ns_map <- stats::setNames(as.numeric(ns), nd)
      out$noise <- as.numeric(ns_map[d])
    }
    return(out)
  }
  lt <- hp[grepl("^l_t_", nm)]
  st <- hp[grepl("^S_t_", nm)]
  ns <- hp[grepl("^noise_", nm)]
  lu <- hp[nm == "l_u_t"]
  lab <- sub("^l_t_", "", names(lt))
  out <- tibble::tibble(
    Output_ID = lab,
    l_t = as.numeric(lt),
    S_t = as.numeric(st[paste0("S_t_", lab)]),
    l_u_t = as.numeric(lu)
  )
  if (length(ns) > 0) {
    out$noise <- as.numeric(ns[paste0("noise_", lab)])
  }
  out
}

#' Test whether a hyper-parameter object is a flat multi-output vector
#'
#' @param hp An object to test.
#' @return A logical value.
#' @keywords internal
is_flat_mo_hp <- function(hp) {
  if (!is.numeric(hp) || is.null(names(hp)) || length(hp) == 0) return(FALSE)
  ## Multi-output flat names: 'l_t_<o>', 'S_t_<o>', 'noise_<o>' (labels may be
  ## arbitrary strings), 'l_u_t' or 'l_u_t_l<q>' (single/multi-latent).
  all(grepl("^(l_t_|S_t_|noise_).+$|^l_u_t(_l.+)?$", names(hp)))
}
