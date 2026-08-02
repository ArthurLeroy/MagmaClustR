quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

match_closest <- function(x, table, tolerance = Inf, nomatch = NA_integer_) {
  lIdx <- findInterval(x, table, rightmost.closed = FALSE, all.inside = TRUE)
  rIdx <- lIdx + 1L

  lIdx[lIdx == 0L] <- 1L

  lDiff <- abs(table[lIdx] - x)
  rDiff <- abs(table[rIdx] - x)

  d <- which(lDiff >= rDiff)

  lIdx[d] <- rIdx[d]

  if (any(is.finite(tolerance))) {
    if (any(tolerance < 0L)) {
      warning(sQuote("tolerance"), " < 0 is meaningless. Set to zero.")
      tolerance[tolerance < 0L] <- 0L
    }

    if (length(nomatch) != 1L) {
      stop("Length of ", sQuote("nomatch"), " has to be one.")
    }

    tolerance <- rep_len(tolerance, length(table))

    lDiff[d] <- rDiff[d]
    lIdx[lDiff > tolerance[lIdx]] <- nomatch
  }

  lIdx
}

#' Detect the task-identifier column of a task-level tibble
#'
#' Task-level objects (raw data, per-task hyper-parameters, per-task mixture
#' probabilities, ...) are keyed either by the legacy 'ID' column or, in the
#' new long-format convention, by 'Task_ID'. This helper lets functions that
#' consume such objects (e.g. plotting functions, [data_allocate_cluster()])
#' remain agnostic to which convention is in use, rather than hardcoding
#' 'ID'.
#'
#' @param df A tibble or data.frame.
#' @return A character string, either \code{"ID"} or \code{"Task_ID"}.
#' @keywords internal
.task_id_col <- function(df) {
  if ("ID" %in% names(df)) {
    "ID"
  } else if ("Task_ID" %in% names(df)) {
    "Task_ID"
  } else {
    stop(
      "Could not find a task-identifier column ('ID' or 'Task_ID') in the ",
      "provided data.",
      call. = FALSE
    )
  }
}

#' Normalise a task-level tibble's identifier column back to 'ID'
#'
#' Inverse of the 'ID' -> 'Task_ID' renaming performed when converting
#' outputs to the new format (see [format_to_new()]). Used internally by
#' functions that *recompute* on hyper-parameters/data pulled from a
#' \code{trained_model} (e.g. [hyperposterior()], [pred_magma()]), so that
#' the deep computational machinery (joins by 'ID', etc.) never needs to be
#' aware of which convention the object is currently in.
#'
#' @param df A tibble or data.frame.
#' @return \code{df}, with 'Task_ID' renamed to 'ID' if present; unchanged
#'   otherwise.
#' @keywords internal
.to_legacy_id <- function(df) {
  if (is.null(df)) return(df)
  if ("Task_ID" %in% names(df) && !("ID" %in% names(df))) {
    df <- df %>% dplyr::rename(ID = "Task_ID")
  }
  df
}

#' Resolve single-output vs multi-output mode for an entry-point function
#'
#' Detects whether a training/prediction call should operate in single-output
#' (SO) or multi-output (MO) mode, from the provided kernel and the
#' cardinality of the 'Output_ID' column (if any), and converts 'data' to the
#' adequate internal working format. Detection relies on the *identity* of
#' the kernel with [convolution_kernel()] (a user-supplied custom kernel
#' function does not, by itself, trigger MO mode) and/or the presence of
#' several distinct 'Output_ID' values in 'data'.
#'
#' If 'multi_output' is explicitly provided and contradicts the detected
#' mode, the function stops with an explicit error rather than silently
#' overriding either signal.
#'
#' @param data A tibble or data.frame, in either the legacy or the new long
#'   format.
#' @param kern A kernel, either a character string or a kernel function.
#' @param multi_output A logical value or NULL (default). If NULL, the mode
#'   is auto-detected. Otherwise, it must agree with the detected mode.
#'
#' @return A list with elements 'data' (converted to the internal working
#'   format when relevant), 'mo' (a logical value indicating whether
#'   multi-output mode was resolved), and 'is_new_format' (a logical value
#'   indicating whether 'data' was originally provided in the new long
#'   format, useful for callers that want to convert their outputs back to
#'   the user's original convention via [format_to_new()]).
#' @keywords internal
resolve_mo_mode <- function(data, kern, multi_output = NULL) {
  is_new_format <- all(c("Task_ID", "Input", "Output") %in% names(data))
  has_multi_outputs <- is_new_format && "Output_ID" %in% names(data) &&
    dplyr::n_distinct(data$Output_ID) > 1
  ## Only the built-in convolution kernel signals multi-output; a custom
  ## user-defined kernel function (single-output use case) must not trigger
  ## the MO branch.
  kern_is_mo <- identical(kern, convolution_kernel)
  detected_mo <- kern_is_mo || has_multi_outputs

  if (is.null(multi_output)) {
    mo <- detected_mo
  } else if (multi_output != detected_mo) {
    stop(
      "'multi_output = ", multi_output, "' is inconsistent with the ",
      "provided data/kernel (detected multi-output = ", detected_mo, "). ",
      "Please either adjust 'multi_output', the kernel, or the 'Output_ID' ",
      "column of 'data'.",
      call. = FALSE
    )
  } else {
    mo <- multi_output
  }

  if (is_new_format) {
    data <- if (mo) .mo_working_format(data) else format_to_legacy(data)
  }

  list(data = data, mo = mo, is_new_format = is_new_format)
}


