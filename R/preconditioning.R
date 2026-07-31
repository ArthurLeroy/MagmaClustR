#' Precondition hyper-parameters via a pooled-likelihood pre-fit
#'
#' Computes a data-driven, reasonably-good starting point for the
#' hyper-parameters used to initialise the EM/VEM algorithms in
#' [train_magma()] and [train_magmaclust()] (an 'iteration 0'). Rather than a
#' purely random draw, (a subset of) the tasks in 'data' are used to fit a
#' single Gaussian Process with hyper-parameters shared across tasks (via
#' [train_shared_gp()]), pooling information across tasks to stabilise the
#' subsequent EM/VEM optimisation and reduce the risk of pathological
#' training runs.
#'
#' @param data A tibble, in the internal working format (column 'ID',
#'   'Input', 'Output', and, for the multi-output case, 'Output_ID').
#' @param kern A kernel, character string or [convolution_kernel()].
#' @param precondition_tasks A number in (0, 1], the fraction of tasks used
#'   for the preconditioning fit. Default 1 (all tasks). Values closer to 0
#'   reduce the preprocessing cost, at the risk of a less representative fit.
#' @param noise A logical value, whether to draw/fit a 'noise'
#'   hyper-parameter.
#'
#' @return A condensed hyper-parameter tibble (no 'ID' column): the result of
#'   [train_shared_gp()] fitted on the (sub)sampled tasks.
#' @keywords internal
precondition_hp <- function(data, kern, precondition_tasks = 1, noise = TRUE) {
  list_ID <- unique(data$ID)
  n_sample <- max(1, round(precondition_tasks * length(list_ID)))

  data_sub <- if (n_sample >= length(list_ID)) {
    data
  } else {
    data %>% dplyr::filter(.data$ID %in% sample(list_ID, n_sample))
  }

  ini <- draw_ini_hp(kern, data_sub, noise = noise)

  quiet(
    train_shared_gp(data = data_sub, ini_hp = ini, kern = kern)
  )
}

#' Broadcast a condensed hyper-parameter tibble onto every task
#'
#' If 'hp_tbl' already provides an 'ID' column (one row, or one set of rows,
#' per task), it is returned as-is. Otherwise, it is treated as a *condensed*
#' hyper-parameter set (shared across tasks, e.g. coming from
#' [precondition_hp()] or a user-supplied shared 'ini_hp') and duplicated for
#' every task in 'list_ID'.
#'
#' @param hp_tbl A hyper-parameter tibble, condensed or not.
#' @param list_ID A vector of task identifiers.
#' @return A hyper-parameter tibble with an 'ID' column covering every task
#'   in 'list_ID'.
#' @keywords internal
broadcast_hp_i <- function(hp_tbl, list_ID) {
  if ("ID" %in% names(hp_tbl)) {
    return(hp_tbl)
  }
  list_ID %>%
    purrr::map(~ dplyr::mutate(hp_tbl, ID = as.character(.x), .before = 1)) %>%
    dplyr::bind_rows()
}

#' Derive data-driven initialisation bounds for lengthscale/variance/noise
#'
#' Translates the observed scale of the inputs and outputs into sensible
#' bounds for the log-scale hyper-parameters drawn by [hp()], instead of
#' fixed defaults that implicitly assume inputs on an \[-1, 1\]-ish domain
#' and outputs of order 1. The lengthscale bound is a fraction of the
#' observed input range (median-heuristic flavoured); the variance and noise
#' bounds are centred on the observed output variance. Only parameters whose
#' interpretation directly scales with the data are affected; other
#' kernels' auxiliary parameters (e.g. 'rq_scale', 'period') are untouched.
#'
#' @param df A tibble or data.frame with columns 'Input' and 'Output'.
#' @return A named list of bounds ('min_val_l', 'max_val_l', 'min_val_S',
#'   'max_val_S', 'min_noise', 'max_noise'), on the internal log-scale shared
#'   by 'se_lengthscale'/'se_variance' and the convolution kernel's
#'   'l_t'/'l_u_t'/'S_t' (all of the form \code{l = 2 * log(lengthscale)} or
#'   \code{S = log(variance)}).
#' @keywords internal
.data_scale_bounds <- function(df) {
  input_range <- suppressWarnings(diff(range(df$Input, na.rm = TRUE)))
  if (!is.finite(input_range) || input_range <= 0) input_range <- 1

  output_var <- stats::var(df$Output, na.rm = TRUE)
  if (!is.finite(output_var) || output_var <= 0) output_var <- 1

  list(
    min_val_l = 2 * log(input_range / 20),
    max_val_l = 2 * log(input_range / 2),
    min_val_S = log(output_var / 4),
    max_val_S = log(output_var * 4),
    min_noise = log(output_var * 0.01),
    max_noise = log(output_var * 0.1)
  )
}

#' Build data-driven per-output generation bounds for the convolution kernel
#'
#' @param data A tibble or data.frame with columns 'Input', 'Output', and
#'   (optionally) 'Output_ID'.
#' @param out_ids A character vector of (sorted, unique) output IDs.
#' @return A tibble of bounds compatible with \code{hp_config}, one row per
#'   output, in \code{out_ids} order.
#' @keywords internal
.data_driven_hp_config <- function(data, out_ids) {
  ## The shared latent-process lengthscale bound is derived from the pooled
  ## (all-output) input range, since it is common to every output.
  global <- .data_scale_bounds(data)

  per_output <- if ("Output_ID" %in% names(data)) {
    split(data, as.character(data$Output_ID))
  } else {
    NULL
  }

  floop <- function(o) {
    df_o <- if (!is.null(per_output)) per_output[[o]] else NULL
    b <- if (is.null(df_o)) global else .data_scale_bounds(df_o)
    tibble::tibble(
      Output_ID = o,
      lt_min = b$min_val_l, lt_max = b$max_val_l,
      St_min = b$min_val_S, St_max = b$max_val_S,
      lu_min = global$min_val_l, lu_max = global$max_val_l,
      noise_min = b$min_noise, noise_max = b$max_noise
    )
  }
  purrr::map_dfr(out_ids, floop)
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
