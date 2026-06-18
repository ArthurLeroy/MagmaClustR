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


#' Pivot MagmaClustR data to long format
#'
#' Converts a simulated dataset from a wide format (e.g. containing an 'ID',
#' 'Input_1', 'Input_2', 'Output') to a standardized long format with
#' 'Task_ID', 'Input_ID', 'Input', 'Output_ID', and 'Output'.
#'
#' @param db A data frame or tibble containing the simulated data.
#'
#' @return A formatted \code{tibble} in long format.
#' @export
format_longer <- function(db) {
  # Rename "ID" in "Task_ID" if necessary
  if ("ID" %in% names(db) && !"Task_ID" %in% names(db)) {
    db <- db %>% dplyr::rename(Task_ID = .data$ID)
  }


  if (!"Output_ID" %in% names(db)) {
    db <- db %>% dplyr::mutate(Output_ID = as.factor(1))
  }

  input_cols <- grep("^Input_\\d+$", names(db), value = TRUE)

  if (length(input_cols) > 0) {
    if ("Input_ID" %in% names(db)) {
      db <- db %>% dplyr::select(-.data$Input_ID)
    }

    db <- db %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(input_cols),
        names_to = "Input_ID",
        values_to = "Input",
        names_prefix = "Input_"
      ) %>%
      dplyr::mutate(Input_ID = as.factor(.data$Input_ID))

  } else if ("Input" %in% names(db)) {
    if (!"Input_ID" %in% names(db)) {
      db <- db %>% dplyr::mutate(Input_ID = as.factor(1))
    } else {
      db <- db %>% dplyr::mutate(Input_ID = as.factor(.data$Input_ID))
    }
  }


  if ("Task_ID" %in% names(db)) {
    db <- db %>% dplyr::mutate(Task_ID = as.factor(.data$Task_ID))
  }
  db <- db %>% dplyr::mutate(Output_ID = as.factor(.data$Output_ID))

  cols_order <- c("Task_ID", "Input_ID", "Input", "Output_ID", "Output")
  present_order <- intersect(cols_order, names(db))
  other_cols <- setdiff(names(db), present_order)

  db %>% dplyr::select(dplyr::all_of(present_order), dplyr::all_of(other_cols))
}

#' Pivot MagmaClustR data to wide format
#'
#' Converts a simulated dataset from the standard long format back to a
#' wider format, extracting dimensions into 'Input_1', 'Input_2', etc.,
#' and renaming 'Task_ID' to 'ID'.
#'
#' @param db A data frame or tibble containing the simulated data in long format.
#'
#' @return A formatted \code{tibble} in wide format.
#' @export
format_wider <- function(db) {
  if ("Input_ID" %in% names(db) && "Input" %in% names(db)) {
    db <- db %>%
      tidyr::pivot_wider(
        names_from = "Input_ID",
        values_from = "Input",
        names_prefix = "Input_"
      )
  }

  if ("Task_ID" %in% names(db)) {
    db <- db %>% dplyr::rename(ID = .data$Task_ID)
  }

  if ("Output_ID" %in% names(db)) {
    if (length(unique(db$Output_ID)) == 1 && as.character(unique(db$Output_ID)[1]) == "1") {
      db <- db %>% dplyr::select(-.data$Output_ID)
    }
  }

  input_cols <- grep("^Input_\\d+$", names(db), value = TRUE)
  input_cols <- input_cols[order(as.numeric(gsub("Input_", "", input_cols)))]

  expected_cols <- c("ID", input_cols, "Output")
  present_cols <- intersect(expected_cols, names(db))
  other_cols <- setdiff(names(db), present_cols)

  db %>% dplyr::select(dplyr::all_of(present_cols), dplyr::all_of(other_cols))
}

