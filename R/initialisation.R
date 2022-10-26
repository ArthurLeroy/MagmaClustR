match_closest <- function(x, table, tolerance=Inf, nomatch=NA_integer_) {
  lIdx <- findInterval(x, table, rightmost.closed=FALSE, all.inside=TRUE)
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

#' Expand a grid of inputs
#'
#' @param Input A vector of inputs.
#' @param ... As many vector of covariates as desired. We advise to give
#'    explicit names when using the function.
#'
#' @return A tibble containing all the combination of values of the
#'    parameters.
#' @export
#'
#' @examples
#' TRUE
expand_grid_inputs <- function(Input, ...) {
  arguments <- list(Input, ...)

  if (!(arguments %>% purrr::every(is.numeric))) {
    stop("The arguments must all be numerical sequences.")
  }

  dim_all <- sapply(arguments, length) %>%
    prod()

  if (dim_all > 10000) {
    warning("The number of grid points is too high. Magma has a cubic ",
            "complexity, so the execution will be extremely long. ",
            "We advise to reduce the length of your grid of inputs."
    )
  }

  tidyr::expand_grid(Input, ...) %>% return()
}

#' Regularise a grid of inputs in a dataset
#'
#' Modify the original grid of inputs to make it more 'regular' (in the sense
#' that the interval between each observation is constant, or corresponds to a
#' specific pattern defined by the user). In particular, this function can also
#' be used to summarise several data points into one, at a specific location. In
#' this case, the output values are averaged according to the 'summarise_fct'
#' argument.
#'
#' @param data A tibble or data frame. Required columns: \code{ID},
#'    \code{Output}. The \code{ID} column contains the unique names/codes used
#'    to identify each individual/task (or batch of data). The \code{Output}
#'    column specifies the observed values (the response variable). The data
#'    frame can also provide as many inputs as desired, with no constraints
#'    on the column names.
#'
#' @param size_grid A number, which indicates how many points each column of the
#'    grid must contain. This number is used when no inputs grid is specified by
#'    the user. See 'grid_inputs' below for details. Default is 30.
#'
#' @param scale A boolean. If TRUE, the dataset is scaled to zero mean and
#'    unit variance in all column. Default is FALSE.
#'
#' @param grid_inputs A data frame, corresponding to a pre-defined grid of
#'    inputs according to which we want to regularise a dataset. If
#'    NULL (default), a dedicated grid of inputs is defined: for each
#'    input column, a regular sequence is created from the min of the input
#'    values to the max, with a step equal to the 'size_grid' parameter.
#'
#' @param summarise_fct A character string or a function. If several similar
#'    inputs are associated with different outputs, the user can choose the
#'    summarising function for the output among the following: min, max, mean,
#'    median. A custom function can be defined if necessary. Default is "mean".
#'
#' @return A data frame, where input columns have been regularised accordingly.
#' @export
#'
#' @examples
#' TRUE
regularize_data <- function(data,
                            size_grid = 30,
                            scale = FALSE,
                            grid_inputs = NULL,
                            summarise_fct = base::mean) {

  if (data %>% is.data.frame()) {
    if (!all(c("ID", "Output") %in% names(data))) {
      stop(
        "The 'data' argument should be a tibble or a data frame containing ",
        "at least the mandatory column names: 'ID', 'Output'"
      )
    }
  } else {
    stop(
      "The 'data' argument should be a tibble or a data frame containing ",
      "at least the mandatory column names: 'ID', 'Output'"
    )
  }
  ## summarize function for data on the same grid node
  if(is.character(summarise_fct)){
    if (summarise_fct == "mean"){
      summarise_fct <- base::mean
    } else if (summarise_fct == "min"){
      summarise_fct <- base::min
    } else if (summarise_fct == "max"){
      summarise_fct <- base::max
    } else if (summarise_fct == "median"){
      summarise_fct <- stats::median
    }
  } else if(!(is.function(summarise_fct))){
    stop("Incorrect type. summarise_fct argument must be either a character or",
         "a function."
    )
  }
  ## Function to scale columns
  scale_function <- function(data) {
    (data - base::mean(data)) / stats::sd(data) %>% return()
  }

  ## Get the Input columns names
  names_col <- data %>%
    dplyr::select(-.data$ID, -.data$Output) %>%
    names()

  ## Scale the Input columns if required
  if (scale) {
    data <- data %>%
      dplyr::mutate_at(tidyselect::all_of(names_col), scale_function)
  }

  if (is.null(grid_inputs)) {
    ## Put the data on a grid node
    fct_round <- function(data, size_grid) {
      round_step <- (base::max(data) - base::min(data)) / size_grid
      data <- data %>%
        plyr::round_any(round_step)
    }

    data %>%
      dplyr::mutate_at(tidyselect::all_of(names_col), fct_round, size_grid) %>%
      dplyr::group_by_at(c("ID", tidyselect::all_of(names_col))) %>%
      dplyr::summarise_all(summarise_fct) %>%
      dplyr::ungroup() %>%
      return()
  } else {
    if (!(setequal(names(grid_inputs), names_col))) {
      stop("Input column names in grid_inputs must be the same as in data.")
    } else {
      round_col <- function(col_name) {
        vector_input <- data %>% dplyr::pull(col_name)
        vector_grid_input <- grid_inputs %>%
          dplyr::pull(col_name) %>%
          unique() %>%
          sort()

        vector_grid_input[match_closest(vector_input, vector_grid_input)] %>%
          return()
      }

      inputs <- sapply(names_col, round_col) %>%
        tibble::as_tibble()
      tibble::tibble(ID = data$ID, inputs, Output = data$Output) %>%
        dplyr::group_by_at(c("ID", tidyselect::all_of(names_col))) %>%
        dplyr::summarise_all(summarise_fct) %>%
        dplyr::ungroup() %>%
        return()
    }
  }
}

#' Run a k-means algorithm to initialise clusters' allocation
#'
#' @param data A tibble containing common Input and associated Output values
#'   to cluster.
#' @param k A number of clusters assumed for running the kmeans algorithm.
#' @param nstart A number, indicating how many re-starts of kmeans are set.
#' @param summary A boolean, indicating whether we want an outcome summary
#'
#' @return A tibble containing the initial clustering obtained through kmeans.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
ini_kmeans <- function(data, k, nstart = 50, summary = FALSE) {
  # if (!identical(
  #   unique(data$Input),
  #   data %>%
  #     dplyr::filter(.data$ID == unique(data$ID)[[1]]) %>%
  #     dplyr::pull(.data$Input)
  # )) {
  floop <- function(i) {
    obs_i <- data %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Output)
    tibble::tibble(
      "ID" = i,
      "Input" = seq_len(3),
      "Output" = c(min(obs_i), mean(obs_i), max(obs_i))
    ) %>%
      return()
  }
  db_regular <- unique(data$ID) %>%
    lapply(floop) %>%
    dplyr::bind_rows() %>%
    dplyr::select(c(.data$ID, .data$Input, .data$Output))
  # } else {
  #   db_regular <- data %>% dplyr::select(c(.data$ID, .data$Input, .data$Output))
  # }

  res <- db_regular %>%
    tidyr::spread(key = .data$Input, value = .data$Output) %>%
    dplyr::select(-.data$ID) %>%
    stats::kmeans(centers = k, nstart = nstart)

  if (summary) {
    res %>% print()
  }

  broom::augment(
    res,
    db_regular %>% tidyr::spread(key = .data$Input, value = .data$Output)
  ) %>%
    dplyr::select(c(.data$ID, .data$.cluster)) %>%
    dplyr::rename(Cluster_ini = .data$.cluster) %>%
    dplyr::mutate(Cluster_ini = paste0("K", .data$Cluster_ini)) %>%
    return()
}


#' Mixture initialisation with kmeans
#'
#' Provide an initial kmeans allocation of the individuals/tasks in a dataset
#' into a definite number of clusters, and return the associated mixture
#' probabilities.
#'
#' @param data A tibble or data frame. Required columns: \code{ID}, \code{Input}
#'    , \code{Output}.
#' @param k A number, indicating the number of clusters.
#' @param name_clust A vector of characters. Each element should correspond to
#'    the name of one cluster.
#' @param nstart A number of restart used in the underlying kmeans algorithm
#'
#' @return A tibble indicating for each \code{ID} in which cluster it belongs
#'    after a kmeans initialisation.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
ini_mixture <- function(data, k, name_clust = NULL, nstart = 50) {
  db_ini <- ini_kmeans(data, k, nstart) %>%
    dplyr::mutate(value = 1) %>%
    tidyr::spread(key = .data$Cluster_ini, value = .data$value, fill = 0)

  if (!is.null(name_clust)) {
    names(db_ini) <- c("ID", name_clust)
  }

  return(db_ini)
}
