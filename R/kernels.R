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

#' Test whether a kernel is a genuinely custom kernel function
#'
#' The built-in [convolution_kernel()] is a function, but (unlike a
#' user-supplied custom kernel) [hp()] knows its parameter structure and can
#' therefore auto-generate initial hyper-parameters for it, just like the
#' character-string kernels. This helper distinguishes the two cases, used
#' throughout the training/prediction functions to decide whether providing
#' 'ini_hp'/'hp' is mandatory.
#'
#' @param kern A kernel, either a character string or a function.
#' @return A logical value.
#' @keywords internal
.is_custom_kernel_fn <- function(kern) {
  is.function(kern) && !identical(kern, convolution_kernel)
}

