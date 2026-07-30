laplace_matching = function(data, likelihood, eps = 1e-2) {

  if (likelihood == 'Bernoulli') {
    ## Compute alpha et beta parameters of the Beta pseudo-likelihood
    ## Then transform into Gaussian observations via Laplace-Matching
    pseudo_data = data %>%
      dplyr::mutate(
        alpha = .data$Output == 1,
        beta = .data$Output == 0
      ) %>%
      dplyr::mutate(Output = log((.data$alpha + eps) / (.data$beta + eps))) %>%
      dplyr::mutate(
        Var_output = (.data$alpha + .data$beta + 2 * eps) /
          ((.data$alpha + eps) * (.data$beta + eps))
      ) %>%
      dplyr::select(c(.data$ID, .data$Input, .data$Output, .data$Var_output))
  }

  return(pseudo_data)
}

revert_laplace_matching = function(sample, likelihood = 'Bernoulli') {
  if (likelihood == 'Bernoulli') {
    ## Apply the logistic transform to revert back to [0,1]
    revert_pred = sample %>%
      dplyr::mutate(Output = 1 / (1 + exp(-.data$Output)))
  }

  return(revert_pred)
}


#' Convert data to the legacy single-output format
#'
#' Before Multi-Output support was introduced, every MagmaClustR function
#' expected a single-output dataset with mandatory columns 'ID', 'Input', and
#' 'Output', where *any* additional column (whatever its name) was implicitly
#' treated as an extra input dimension. This legacy convention is not
#' directly compatible with the newer, unified long format
#' ('Task_ID'/'Input_ID'/'Input'/'Output_ID'/'Output'), which does not
#' constrain the number of outputs.
#'
#' `format_to_legacy()` bridges the two, to keep using functions that have
#' not yet been updated for Multi-Output data. Because the legacy format has
#' no notion of multiple outputs, a dataset with several distinct
#' `Output_ID` values is split into one legacy dataset per output.
#'
#' @param db A tibble or data.frame in the new long format, containing at
#'   least 'Task_ID', 'Input', and 'Output' columns (and, if relevant,
#'   'Input_ID' and 'Output_ID').
#' @param keep_extra_cols A logical value. Columns other than
#'   Task_ID/Input_ID/Input/Output_ID/Output (e.g. simulation hyperparameters
#'   such as those added by `simu_data(add_hp = TRUE)`) would be silently
#'   reinterpreted as extra input dimensions by legacy code. By default
#'   (`FALSE`), such columns are dropped (with a message). Set to `TRUE` to
#'   keep them anyway, at your own risk.
#'
#' @return If `db` contains a single output, a tibble with columns 'ID',
#'   'Input' (plus 'Input2', 'Input3', ... when several input dimensions are
#'   present), and 'Output'. If `db` contains several outputs, a *named list*
#'   of such tibbles (one per `Output_ID` level).
#'
#' @export
#'
#' @examples
#' db_so <- simu_data(n_outputs = 1)
#' format_to_legacy(db_so)
#'
#' db_mo <- simu_data(n_outputs = 2)
#' format_to_legacy(db_mo)
format_to_legacy <- function(db, keep_extra_cols = FALSE) {
  required_cols <- c("Task_ID", "Input", "Output")
  if (!all(required_cols %in% names(db))) {
    stop(
      "'db' should contain at least the columns: ",
      paste(required_cols, collapse = ", ")
    )
  }

  if (!"Input_ID" %in% names(db)) {
    db <- db %>% dplyr::mutate(Input_ID = factor(1))
  }
  if (!"Output_ID" %in% names(db)) {
    db <- db %>% dplyr::mutate(Output_ID = factor(1))
  }

  ## Any column beyond the 5 standardised ones would be silently treated as
  ## an extra input dimension by legacy code: drop them by default.
  core_cols <- c("Task_ID", "Input_ID", "Input", "Output_ID", "Output")
  extra_cols <- setdiff(names(db), core_cols)
  if (length(extra_cols) > 0) {
    if (keep_extra_cols) {
      message(
        "Keeping non-standard column(s) that will be treated as extra ",
        "input dimensions by legacy code: ", paste(extra_cols, collapse = ", ")
      )
    } else {
      message(
        "Dropping column(s) not part of the legacy format (they would ",
        "otherwise be treated as extra input dimensions): ",
        paste(extra_cols, collapse = ", "),
        ". Use 'keep_extra_cols = TRUE' to keep them."
      )
      db <- db %>% dplyr::select(dplyr::all_of(core_cols))
    }
  }

  ## Pivot each output separately, and use the legacy naming convention for
  ## additional input dimensions: 'Input', 'Input2', 'Input3', ...
  to_legacy_single_output <- function(sub_db) {
    input_ids <- sort(unique(as.integer(as.character(sub_db$Input_ID))))
    legacy_names <- c("Input", paste0("Input", input_ids[-1]))
    names(legacy_names) <- input_ids

    sub_db %>%
      dplyr::mutate(
        Input_ID = legacy_names[as.character(.data$Input_ID)]
      ) %>%
      dplyr::select(-.data$Output_ID) %>%
      tidyr::pivot_wider(names_from = "Input_ID", values_from = "Input") %>%
      dplyr::rename(ID = .data$Task_ID) %>%
      dplyr::select(.data$ID, dplyr::all_of(unname(legacy_names)), dplyr::everything())
  }

  output_levels <- unique(db$Output_ID)

  if (length(output_levels) == 1) {
    return(db %>% to_legacy_single_output())
  }

  db %>%
    split(db$Output_ID) %>%
    purrr::map(to_legacy_single_output)
}

