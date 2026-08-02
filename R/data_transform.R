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
      dplyr::select(c("ID", "Input", "Output", "Var_output"))
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
      dplyr::select(-"Output_ID") %>%
      tidyr::pivot_wider(names_from = "Input_ID", values_from = "Input") %>%
      dplyr::rename(ID = .data$Task_ID) %>%
      dplyr::select("ID", dplyr::all_of(unname(legacy_names)), dplyr::everything())
  }

  output_levels <- unique(db$Output_ID)

  if (length(output_levels) == 1) {
    return(db %>% to_legacy_single_output())
  }

  db %>%
    split(db$Output_ID) %>%
    purrr::map(to_legacy_single_output)
}

#' Convert the new long format to the internal multi-output working format
#'
#' Unlike [format_to_legacy()], which bridges the new long format to the
#' legacy single-output convention, this helper keeps the multi-output
#' structure ('Output_ID') intact, while renaming 'Task_ID' to 'ID' (the
#' naming convention expected by the internal training/prediction code) and
#' pivoting a multi-level 'Input_ID' to wide 'Input'/'Input2'/... columns
#' (mirroring [format_to_legacy()]'s treatment of multiple input dimensions).
#'
#' @param db A tibble or data.frame in the new long format.
#' @return A tibble in the internal multi-output working format (columns
#'   'ID', 'Input' (plus 'Input2', ... when relevant), 'Output_ID', 'Output').
#' @keywords internal
.mo_working_format <- function(db) {
  if (!"Input_ID" %in% names(db)) {
    return(db %>% dplyr::rename(ID = "Task_ID"))
  }

  input_ids <- sort(unique(as.integer(as.character(db$Input_ID))))
  if (length(input_ids) <= 1) {
    return(db %>% dplyr::select(-"Input_ID") %>% dplyr::rename(ID = "Task_ID"))
  }

  legacy_names <- c("Input", paste0("Input", input_ids[-1]))
  names(legacy_names) <- input_ids

  db %>%
    dplyr::mutate(Input_ID = legacy_names[as.character(.data$Input_ID)]) %>%
    tidyr::pivot_wider(names_from = "Input_ID", values_from = "Input") %>%
    dplyr::rename(ID = "Task_ID") %>%
    dplyr::select(
      "ID",
      dplyr::all_of(unname(legacy_names)),
      "Output_ID",
      dplyr::everything()
    )
}

#' Convert an internal-format object back to the new long-format convention
#'
#' Inverse of [format_to_legacy()]/[.mo_working_format()]: renames 'ID' to
#' 'Task_ID' (if present) and pivots 'Input'/'Input2'/... wide columns back
#' to long 'Input_ID'/'Input' rows (if present), so that objects derived from
#' training/prediction (e.g. \code{trained_model$hp_i},
#' \code{trained_model$ini_args$data}, \code{trained_model$hyperpost$mean})
#' can be displayed to users in the same convention they used as input,
#' without exposing the legacy internal representation. A 'Reference' column,
#' if present, is dropped (it is an internal bookkeeping artefact).
#'
#' Objects with neither an 'ID' nor an 'Input'/'Input2'/... column (e.g. a
#' hyper-parameter tibble already keyed by 'Output_ID' only) are returned
#' unchanged, aside from the optional 'Output_ID' addition below.
#'
#' @param x A tibble or data.frame, in the internal working format.
#' @param output_id An optional value (e.g. \code{"1"}) used to add an
#'   'Output_ID' column when \code{x} doesn't already have one (typically for
#'   single-output objects, to mirror the 'Output_ID' column present in the
#'   user's original new-format data).
#'
#' @return A tibble in the new long-format convention.
#' @export
#'
#' @examples
#' db <- simu_data(n_tasks = 3)
#' model <- train_magma(data = db, n_iter_max = 1)
#' format_to_new(model$hp_i)
format_to_new <- function(x, output_id = NULL) {
  if ("Reference" %in% names(x)) {
    x <- x %>% dplyr::select(-"Reference")
  }
  if ("ID" %in% names(x)) {
    x <- x %>% dplyr::rename(Task_ID = "ID")
  }

  input_cols <- grep("^Input[0-9]*$", names(x), value = TRUE)
  if (length(input_cols) == 0) {
    if (!is.null(output_id) && !"Output_ID" %in% names(x)) {
      x <- x %>% dplyr::mutate(Output_ID = as.character(output_id))
    }
    return(x)
  }

  input_ids <- ifelse(input_cols == "Input", "1", sub("^Input", "", input_cols))

  x <- x %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(input_cols),
      names_to = "Input_ID",
      values_to = "Input"
    ) %>%
    dplyr::mutate(
      Input_ID = input_ids[match(.data$Input_ID, input_cols)]
    ) %>%
    tidyr::drop_na("Input")

  if (!is.null(output_id) && !"Output_ID" %in% names(x)) {
    x <- x %>% dplyr::mutate(Output_ID = as.character(output_id))
  }

  front <- intersect(c("Task_ID", "Input_ID", "Input", "Output_ID"), names(x))
  x %>% dplyr::select(dplyr::all_of(front), dplyr::everything())
}

