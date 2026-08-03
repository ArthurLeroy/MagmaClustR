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

#' Build the prediction grid of inputs (single/multi-output, single/multi-input)
#'
#' Centralises the resolution of the user-facing \code{grid_inputs} argument
#' shared by [pred_gp()], [pred_magma()], [pred_magmaclust()],
#' [hyperposterior()] and [hyperposterior_clust()]. It turns \code{grid_inputs}
#' (NULL, a numeric vector, or a data frame in the new front-end format) into a
#' tibble of prediction inputs expressed in the internal working convention,
#' with a 'Reference' key built exactly like the one attached to the (already
#' converted) \code{data}.
#'
#' Resolution rules:
#' - \code{NULL}: an automatic regular grid is built on the numeric coordinate
#'   dimension(s), spanning the observed range of each in \code{data}.
#' - a numeric vector: interpreted as the grid of 'Input' values; only valid
#'   when the data have a single coordinate dimension.
#' - a data frame: read from the new front-end format. If an 'Input_ID' column
#'   is present, its (sorted, integer) levels are mapped to 'Input', 'Input2',
#'   ... and the coordinate grid is the Cartesian product of the distinct
#'   'Input' values per level. If 'Input_ID' is absent, the data frame is
#'   assumed to already use the working/legacy coordinate columns (e.g. 'Input',
#'   or 'Input' plus a named covariate). Any other column is ignored.
#'
#' In multi-output mode, the coordinate grid is replicated across the
#' 'Output_ID' values supplied in \code{grid_inputs} if present, otherwise
#' across every output observed in \code{data}.
#'
#' @param grid_inputs NULL, a numeric vector, or a data frame (new front-end
#'   format).
#' @param data A tibble in the internal working format (already converted by
#'   [resolve_mo_mode()]), used both to define the default grid ranges and, in
#'   multi-output mode, the set of observed outputs.
#' @param names_col A character vector of the input column names of the working
#'   format (coordinate dimensions, plus 'Output_ID' in multi-output mode).
#' @param fn_name A character string naming the calling function, used only to
#'   craft error messages.
#'
#' @return A tibble with a 'Reference' column followed by the \code{names_col}
#'   columns, sorted by 'Reference'.
#' @keywords internal
.build_grid_inputs <- function(grid_inputs, data, names_col,
                               fn_name = "pred_magma") {
  coord_cols <- setdiff(names_col, "Output_ID")
  has_output_id <- "Output_ID" %in% names_col

  ## 0. Fast path: an already-resolved working grid (carries 'Reference' and
  ## every input column) is passed through unchanged. This covers internal
  ## recompute calls (e.g. pred_magma() -> hyperposterior() with the union of
  ## observed and prediction inputs) without re-crossing its rows.
  if (is.data.frame(grid_inputs) &&
      "Reference" %in% names(grid_inputs) &&
      all(names_col %in% names(grid_inputs))) {
    out <- grid_inputs %>%
      dplyr::select("Reference", tidyselect::all_of(names_col)) %>%
      dplyr::distinct()
    if (has_output_id) {
      out <- out %>% dplyr::mutate(Output_ID = as.character(.data$Output_ID))
    }
    return(out %>% dplyr::arrange(.data$Reference))
  }

  ## 1. Build the coordinate grid (columns == coord_cols)
  if (is.null(grid_inputs)) {
    size_grid <- if (length(coord_cols) == 1) {
      500
    } else {
      round(1000^(1 / length(coord_cols)))
    }
    seqs <- stats::setNames(
      lapply(coord_cols, function(cc) {
        seq(min(data[[cc]]), max(data[[cc]]), length.out = size_grid)
      }),
      coord_cols
    )
    coord_grid <- tibble::as_tibble(expand.grid(seqs))
  } else if (is.vector(grid_inputs)) {
    if (length(coord_cols) != 1) {
      stop(
        "'grid_inputs' was provided as a vector but the data have ",
        length(coord_cols), " input dimensions. Please provide a data frame ",
        "with 'Input_ID' and 'Input' columns instead (in '", fn_name, "').",
        call. = FALSE
      )
    }
    coord_grid <- tibble::tibble(Input = grid_inputs)
  } else if (is.data.frame(grid_inputs)) {
    if ("Input_ID" %in% names(grid_inputs)) {
      ## New long format -> working columns Input/Input2/... (Cartesian
      ## product of the distinct 'Input' values per 'Input_ID' level).
      ids <- sort(unique(as.integer(as.character(grid_inputs$Input_ID))))
      nms <- c("Input", paste0("Input", ids[-1]))
      seqs <- lapply(ids, function(i) {
        sort(unique(grid_inputs$Input[
          as.integer(as.character(grid_inputs$Input_ID)) == i
        ]))
      })
      coord_grid <- tibble::as_tibble(expand.grid(stats::setNames(seqs, nms)))
    } else {
      if (!all(coord_cols %in% names(grid_inputs))) {
        stop(
          "'grid_inputs' must provide the coordinate column(s): ",
          paste(coord_cols, collapse = ", "), " (in '", fn_name, "').",
          call. = FALSE
        )
      }
      coord_grid <- grid_inputs %>%
        dplyr::select(tidyselect::all_of(coord_cols)) %>%
        dplyr::distinct()
    }
    if (!setequal(names(coord_grid), coord_cols)) {
      stop(
        "The input dimensions of 'grid_inputs' do not match those of the ",
        "data (expected: ", paste(coord_cols, collapse = ", "), ").",
        call. = FALSE
      )
    }
  } else {
    stop(
      "'grid_inputs' should be NULL, a numeric vector or a data frame ",
      "(in '", fn_name, "').",
      call. = FALSE
    )
  }

  ## 2. Round for numerical stability / Reference consistency with 'data'
  coord_grid <- coord_grid %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), signif))

  ## 3. Replicate across outputs (multi-output)
  inputs_pred <- if (has_output_id) {
    outputs <- if (
      is.data.frame(grid_inputs) && "Output_ID" %in% names(grid_inputs)
    ) {
      unique(as.character(grid_inputs$Output_ID))
    } else {
      unique(as.character(data$Output_ID))
    }
    tidyr::crossing(Output_ID = outputs, coord_grid)
  } else {
    coord_grid
  }

  ## 4. Build the 'Reference' key exactly as done for 'data'
  inputs_pred %>%
    tidyr::unite(
      "Reference",
      tidyselect::all_of(names_col),
      sep = ":",
      remove = FALSE
    ) %>%
    dplyr::arrange(.data$Reference) %>%
    dplyr::select("Reference", tidyselect::all_of(names_col))
}


