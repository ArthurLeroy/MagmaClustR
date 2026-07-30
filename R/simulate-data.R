#' Draw a number
#'
#' Draw uniformly a number within a specified interval
#'
#' @param int An interval of values we want to draw uniformly in.
#'
#' @return A 2-decimals-rounded random number
#'
#' @keywords internal
#'
#' @examples
#' TRUE
draw <- function(int) {
  int <- sort(int)

  stats::runif(1, int[1], int[2]) %>%
    round(2) %>%
    return()
}

#' Get the number of points in a grid object
#'
#' Evaluates the size of a grid, whether it is formatted as a simple vector
#' or a multidimensional data.frame.
#'
#' @param grid_obj A numeric vector, list, or data.frame representing the input grid.
#'
#' @return An integer corresponding to the length of the vector or the number of rows
#'   of the data.frame.
#'
#' @keywords internal
grid_n_points <- function(grid_obj) {
  if (is.data.frame(grid_obj)) {
    return(nrow(grid_obj))
  } else {
    return(length(grid_obj))
  }
}

#' Coerce a grid object into a tibble
#'
#' Ensures that the provided grid object is consistently formatted as a tibble.
#' If the input is a simple numeric vector, it is converted to a tibble with a single
#' column named 'Input'.
#'
#' @param grid_obj A numeric vector or data.frame representing the input grid.
#'
#' @return A \code{tibble} containing the grid data.
#'
#' @keywords internal
coerce_grid_to_df <- function(grid_obj) {
  if (is.data.frame(grid_obj)) {
    return(tibble::as_tibble(grid_obj))
  } else {
    return(tibble::tibble(Input = grid_obj))
  }
}

#' Sort a grid object
#'
#' Sorts the elements of a grid. If the grid is a data.frame, it orders the rows
#' based on the values of its columns (from the first to the last column).
#'
#' @param grid_obj A numeric vector or data.frame representing the input grid.
#'
#' @return The sorted grid object, maintaining its original structure (vector or data.frame).
#'
#' @keywords internal
sort_grid_object <- function(grid_obj) {
  if (is.data.frame(grid_obj)) {
    ordered_indices <- do.call(order, grid_obj)
    return(grid_obj[ordered_indices, , drop = FALSE])
  }

  sort(grid_obj)
}

#' Sample and sort points from a grid object
#'
#' Randomly selects a specified number of points from a grid without replacement,
#' and returns them in a sorted order. It automatically handles both 1D vectors
#' and multi-dimensional data.frames.
#'
#' @param grid_obj A numeric vector or data.frame representing the input grid.
#' @param size An integer specifying the number of points to randomly sample.
#'
#' @return A sorted grid object of the same class as \code{grid_obj}, containing
#'   the sampled points.
#'
#' @keywords internal
sample_grid_points <- function(grid_obj, size) {
  if (is.data.frame(grid_obj)) {
    sampled_indices <- sample(seq_len(nrow(grid_obj)), size, replace = FALSE)
    sampled_grid <- grid_obj[sampled_indices, , drop = FALSE]
    return(sort_grid_object(sampled_grid))
  }

  sample(grid_obj, size, replace = FALSE) %>% sort()
}

#' Build one candidate grid per output
#'
#' Provides a single, output-count-agnostic strategy for constructing the
#' working input grid of every output: for each output, \code{size} points are
#' sampled *without replacement* from a candidate pool of points, then sorted.
#' The very same discrete-sampling strategy is applied whether there is a
#' single output or several, and whether the grid is a 1D numeric vector or a
#' multidimensional data.frame, which keeps the number of points effectively
#' generated for the mean process and for each task consistent and
#' predictable (no separate, precision-based continuous re-sampling is used
#' for the multi-output case anymore).
#'
#' @param grid_input A numeric vector, a data.frame, or a list of such objects
#'   (one candidate pool per output). If a single vector/data.frame is
#'   provided while several outputs are requested, the same pool is reused
#'   (independently, unless `shared_grid_outputs = TRUE`) for every output.
#' @param points_per_output A numeric vector giving the number of points to
#'   sample for each output.
#' @param shared_grid_outputs A logical value. If `TRUE`, a single grid is
#'   sampled once and shared identically across all outputs.
#' @param precomputed_grid_list An optional list of precomputed grids (one per
#'   output), used verbatim instead of sampling.
#'
#' @return A list of grid objects (vectors or data.frames), one per output.
#'
#' @keywords internal
build_output_grid_list <- function(grid_input,
                                    points_per_output,
                                    shared_grid_outputs = FALSE,
                                    precomputed_grid_list = NULL) {
  num_outputs <- length(points_per_output)

  if (!is.null(precomputed_grid_list)) {
    if (length(precomputed_grid_list) != num_outputs) {
      stop("'precomputed_grid_list' length must match the number of outputs.")
    }
    return(precomputed_grid_list)
  }

  # A distinct candidate pool per output is only used when 'grid_input' is
  # itself a plain list (not a data.frame) of the expected length.
  is_per_output_pool <- is.list(grid_input) &&
    !is.data.frame(grid_input) &&
    length(grid_input) == num_outputs

  pool_list <- if (is_per_output_pool) {
    grid_input
  } else {
    rep(list(grid_input), num_outputs)
  }

  if (shared_grid_outputs) {
    if (length(unique(points_per_output)) > 1) {
      stop("For a shared grid across outputs, 'points_per_output' must be ",
           "identical for every output.")
    }
    shared_grid <- sample_grid_points(pool_list[[1]], points_per_output[1])
    return(rep(list(shared_grid), num_outputs))
  }

  purrr::map2(pool_list, points_per_output, sample_grid_points)
}

#' Build unique reference names for a grid of (possibly multi-output) inputs
#'
#' Creates a character vector uniquely identifying each row of an input
#' data.frame containing an `Output_ID` column, by combining the output
#' identifier with its numeric coordinates. Used to keep a consistent naming
#' convention between the mean process grid and any task-specific sub-grid,
#' regardless of the kernel used to generate the data (SE or convolution).
#'
#' @param input_df A data.frame containing numeric coordinate columns and an
#'   `Output_ID` column.
#'
#' @return A character vector of unique reference names, one per row.
#'
#' @keywords internal
build_reference_names <- function(input_df) {
  coord_cols <- names(input_df)[
    sapply(input_df, is.numeric) & names(input_df) != "Output_ID"
  ]
  pasted_coords <- do.call(paste, c(input_df[coord_cols], sep = ";"))
  paste0("o", input_df$Output_ID, ";", pasted_coords)
}

#' Pivot the wide coordinate columns of a simulated dataset into the
#' standardised long format
#'
#' Internally, the mean process and per-task realizations are built with one
#' column per input dimension ('Input', 'Input2', 'Input3', ...), because the
#' covariance computations need the coordinates side by side. This helper
#' performs the final pivot into the package's standardised long format
#' ('Input_ID' + 'Input'), dropping the placeholder `Input_ID` column set
#' during task generation and rebuilding it from the actual coordinate
#' columns present.
#'
#' @param db A tibble containing 'Task_ID', 'Output_ID', 'Output', and one or
#'   several coordinate columns named 'Input', 'Input2', 'Input3', ...
#'
#' @return A tibble in the standardised long format, with 'Input_ID' and
#'   'Input' columns instead of the wide coordinate columns.
#'
#' @keywords internal
pivot_inputs_longer <- function(db) {
  input_cols <- grep("^Input\\d*$", names(db), value = TRUE)

  input_ids <- ifelse(input_cols == "Input", "1", sub("^Input", "", input_cols))
  id_map <- stats::setNames(input_ids, input_cols)
  id_levels <- as.character(sort(as.integer(unique(input_ids))))

  db %>%
    dplyr::select(-dplyr::any_of("Input_ID")) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(input_cols),
      names_to = "Input_ID",
      values_to = "Input"
    ) %>%
    dplyr::mutate(
      Input_ID = factor(id_map[.data$Input_ID], levels = id_levels)
    ) %>%
    dplyr::select(
      .data$Task_ID, .data$Input_ID, .data$Input, .data$Output_ID, .data$Output,
      dplyr::everything()
    )
}

#' Generate the prior mean of the mean process, for every output
#'
#' If `prior_means` is provided, it is used as-is (one scalar, or one full
#' vector matching an output's grid length, per output). Otherwise, a linear
#' trend (randomly drawn intercept and slope, applied along the first
#' coordinate of each output's grid) is generated for every output, providing
#' a single, output-count-agnostic default.
#'
#' @param grid_list A list of grid objects (one per output), as returned by
#'   \code{build_output_grid_list()}.
#' @param prior_means An optional vector or list, providing one element per
#'   output (a scalar, or a vector matching the output's grid length).
#' @param m0_intercept A vector of 2 numbers. Interval for the intercept of
#'   the default linear trend.
#' @param m0_slope A vector of 2 numbers. Interval for the slope of the
#'   default linear trend.
#'
#' @return A list of length `length(grid_list)`, each element being the prior
#'   mean (scalar or vector) associated with the corresponding output.
#'
#' @keywords internal
generate_prior_means <- function(grid_list,
                                  prior_means = NULL,
                                  m0_intercept = c(-50, 50),
                                  m0_slope = c(-5, 5)) {
  num_outputs <- length(grid_list)

  if (!is.null(prior_means)) {
    if (length(prior_means) != num_outputs) {
      stop("'prior_means' must provide exactly one element per output.")
    }
    return(as.list(prior_means))
  }

  purrr::map(grid_list, function(grid) {
    coords <- coerce_grid_to_df(grid)
    draw(m0_intercept) + draw(m0_slope) * coords[[1]]
  })
}

#' Build a task or cluster identifier
#'
#' @param i An integer, the task index within its cluster.
#' @param k An integer, the cluster index.
#' @param n_clusters An integer, the total number of clusters being simulated.
#'
#' @return A character string identifying the task.
#'
#' @keywords internal
make_task_id <- function(i, k, n_clusters) {
  if (n_clusters > 1) {
    paste0("ID", i, "-Clust", k)
  } else {
    as.character(i)
  }
}

#' Post-process a single cluster's simulated dataset
#'
#' Optionally drops the hyperparameter columns and/or adds a 'Cluster' column,
#' factorising the two `add_hp`/`add_clust` steps that were otherwise repeated
#' identically for every code path in \code{simu_data()}.
#'
#' @param db A tibble, the simulated data for one cluster.
#' @param add_hp A logical value. If `FALSE`, hyperparameter columns are
#'   dropped.
#' @param add_clust A logical value. If `TRUE`, a 'Cluster' column is added.
#' @param k An integer, the cluster index used to populate the 'Cluster'
#'   column.
#' @param hp_cols A character vector naming the hyperparameter columns to
#'   drop when `add_hp = FALSE`.
#'
#' @return The post-processed tibble.
#'
#' @keywords internal
postprocess_cluster_db <- function(db, add_hp, add_clust, k, hp_cols) {
  if (!add_hp) {
    db <- db %>% dplyr::select(-dplyr::any_of(hp_cols))
  }
  if (add_clust) {
    db <- db %>% dplyr::mutate(Cluster = k)
  }
  db
}


#' @title Generate a Mean Process Realization

#' Creates a working grid for each output and generates a single smooth
#' underlying mean process, `mu_0`, by drawing a realization from a Gaussian
#' Process. A single output uses the Squared Exponential ("SE") kernel;
#' several outputs jointly use a multi-output convolution kernel. Apart from
#' this unavoidable kernel-choice fork (a single-output convolution kernel is
#' not equivalent to a standalone SE kernel, and \code{kern_to_cov()} itself
#' requires at least 2 distinct outputs to enable its vectorized multi-output
#' mode), every other step (grid sampling, prior mean, final tibble format) is
#' shared between the two cases.
#'
#' @param points_per_output A numeric vector specifying the number of points
#'    to sample for each output's working grid.
#' @param grid_input A numeric vector, a data.frame, or a list of such
#'    objects (one candidate pool per output), defining the input domain(s)
#'    from which the working grid(s) are sampled.
#' @param prior_means An optional vector or list providing one element (a
#'    scalar, or a full vector matching the output's grid length) per output.
#'    If `NULL` (default), a random linear trend is drawn for each output.
#' @param m0_intercept A vector of 2 numbers. Interval for the intercept of
#'    the default linear trend (used only when `prior_means` is `NULL`).
#' @param m0_slope A vector of 2 numbers. Interval for the slope of the
#'    default linear trend (used only when `prior_means` is `NULL`).
#' @param hp_config_mean_process A tibble configuring the hyperparameters for
#'    the mean process's convolution kernel. Required when there is more than
#'    one output.
#' @param noise_0 An optional numeric vector specifying the log-variance of
#'    the noise for each output (multi-output case only). If provided, its
#'    length must match the number of outputs. Default is `NULL` (no noise).
#' @param shared_grid_outputs A logical value. If `TRUE`, all outputs are
#'    defined on the exact same input grid.
#' @param shared_hp_outputs A logical value. If `TRUE`, all outputs share the
#'    same hyperparameters (multi-output case only).
#' @param precomputed_grid_list An optional list of precomputed grids, one per
#'    output.
#' @param v A number. The variance hyper-parameter of the SE kernel (single
#'    output case only).
#' @param l A number. The lengthscale hyper-parameter of the SE kernel
#'    (single output case only).
#' @param sigma A number. The noise variance added on the diagonal (single
#'    output case only).
#'
#' @return A list containing the realization (as both a named vector and a
#'   tibble), the grid used for each output, and the drawn hyperparameters.
#'
#' @keywords internal
#'
generate_mean_process <- function(
    points_per_output,
    grid_input,
    prior_means = NULL,
    m0_intercept = c(-50, 50),
    m0_slope = c(-5, 5),
    hp_config_mean_process = NULL,
    noise_0 = NULL,
    shared_grid_outputs = FALSE,
    shared_hp_outputs = FALSE,
    precomputed_grid_list = NULL,
    v = NULL,
    l = NULL,
    sigma = 0
) {
  num_outputs <- length(points_per_output)

  ## 1. Sample one working grid per output, using a single discrete-sampling
  ## strategy regardless of dimensionality or number of outputs.
  grid_list <- build_output_grid_list(
    grid_input = grid_input,
    points_per_output = points_per_output,
    shared_grid_outputs = shared_grid_outputs,
    precomputed_grid_list = precomputed_grid_list
  )

  ## 2. Build the prior mean vector, shared logic for any number of outputs.
  prior_means_list <- generate_prior_means(
    grid_list,
    prior_means = prior_means,
    m0_intercept = m0_intercept,
    m0_slope = m0_slope
  )
  mean_list <- purrr::map2(grid_list, prior_means_list, function(grid, pm) {
    if (length(pm) == 1) {
      return(rep(pm, grid_n_points(grid)))
    }
    if (length(pm) != grid_n_points(grid)) {
      stop("Each element of 'prior_means' must be a scalar or match its ",
           "output's grid length.")
    }
    pm
  })
  m0_mean_function <- unlist(mean_list)

  ## 3. Build the full input data.frame and the consistent reference names.
  input_df_mean_process <- purrr::imap_dfr(
    grid_list,
    ~coerce_grid_to_df(.x) %>% dplyr::mutate(Output_ID = as.factor(.y))
  )
  reference <- build_reference_names(input_df_mean_process)

  ## 4. Compute the covariance matrix. This is the only step where the
  ## single-output and multi-output cases genuinely differ: they rely on two
  ## different kernel families (see title note above).
  hp_draws <- NULL

  if (num_outputs == 1) {
    grid_coords <- input_df_mean_process %>% dplyr::select(-.data$Output_ID)
    K_theta0_X <- kern_to_cov(
      grid_coords,
      "SE",
      tibble::tibble(se_variance = v, se_lengthscale = l)
    ) + diag(sigma, nrow(input_df_mean_process))
  } else {
    if (is.null(hp_config_mean_process)) {
      stop("'hp_config_mean_process' is required when there is more than ",
           "one output.")
    }

    lu0_shared <- stats::runif(1,
                               hp_config_mean_process$lu0_min[1],
                               hp_config_mean_process$lu0_max[1])
    lu0_vals <- rep(lu0_shared, num_outputs)

    if (shared_hp_outputs) {
      l0_shared <- stats::runif(1,
                                hp_config_mean_process$l0_min[1],
                                hp_config_mean_process$l0_max[1])
      S0_shared <- stats::runif(1,
                                hp_config_mean_process$S0_min[1],
                                hp_config_mean_process$S0_max[1])
      l0_vals <- rep(l0_shared, num_outputs)
      S0_vals <- rep(S0_shared, num_outputs)
    } else {
      l0_vals <- purrr::map2_dbl(hp_config_mean_process$l0_min,
                                 hp_config_mean_process$l0_max,
                                 ~stats::runif(1, .x, .y))
      S0_vals <- purrr::map2_dbl(hp_config_mean_process$S0_min,
                                 hp_config_mean_process$S0_max,
                                 ~stats::runif(1, .x, .y))
    }

    hp_tibble_for_kernel <- tibble::tibble(
      Output_ID = as.factor(seq_len(num_outputs)),
      l_t       = l0_vals,
      S_t       = S0_vals,
      l_u_t     = lu0_vals
    )

    K_theta0_X <- suppressWarnings(kern_to_cov(
      input = input_df_mean_process,
      kern = convolution_kernel_KD,
      hp = hp_tibble_for_kernel
    ))

    if (!is.null(noise_0)) {
      if (length(noise_0) != num_outputs) {
        stop(sprintf(
          "'noise_0' length should match the number of outputs (%d).",
          num_outputs
        ))
      }
      diag(K_theta0_X) <- diag(K_theta0_X) +
        exp(noise_0[input_df_mean_process$Output_ID])
    }

    hp_draws <- list(l0 = l0_vals, S0 = S0_vals, lu0 = lu0_vals)
  }

  rownames(K_theta0_X) <- colnames(K_theta0_X) <- reference
  names(m0_mean_function) <- reference

  ## 5. Draw the mean process realization, shared logic for any number of
  ## outputs.
  mean_process_realization <- MASS::mvrnorm(
    n = 1,
    mu = m0_mean_function,
    Sigma = K_theta0_X
  )
  names(mean_process_realization) <- reference

  mean_process_df <- input_df_mean_process %>%
    dplyr::mutate(Output = as.numeric(mean_process_realization))

  list(
    mean_process_realization = mean_process_realization,
    mean_process_df = mean_process_df,
    grid_list = grid_list,
    points_per_output = points_per_output,
    hp_draws = hp_draws,
    noise_0 = noise_0
  )
}


#' @title Generate Data for a Single Task

#' Generates a simulated dataset for a single, specific task, centered around
#' the mean process. This function samples a sparse sub-grid from a provided
#' mean process grid, computes a task-specific covariance matrix (using
#' `kern_to_cov`), and then draws a realization from a multivariate normal
#' distribution. As in \code{generate_mean_process()}, the only genuine fork
#' between the single-output and multi-output cases is the choice of kernel.
#'
#' @param task_id An identifier (numeric or string) for the task being
#'   generated.
#' @param mean_process_info A list object containing information from the
#'   shared mean process (output from \code{generate_mean_process()}). Must
#'   contain `grid_list` (a list of input grids, one per output) and
#'   `mean_process_realization` (a named vector of the full mean process
#'   realization).
#' @param task_hp_tibble A `tibble` containing the specific hyperparameters
#'   for this task. Required in the multi-output case; passed directly to the
#'   `hp` argument of `kern_to_cov`.
#' @param n_points_per_output A numeric vector giving the number of points to
#'   sample, for each output, from that output's mean-process grid.
#' @param shared_grid_outputs A logical value. If `TRUE`, all outputs in the
#'   task are observed on the same sampled inputs.
#' @param precomputed_task_grid_list An optional list of precomputed task
#'   grids, one per output. If provided, the same sampled inputs are reused
#'   instead of drawing a new sub-grid for the task.
#' @param v A number. The variance hyper-parameter of the SE kernel (single
#'   output case only).
#' @param l A number. The lengthscale hyper-parameter of the SE kernel
#'   (single output case only).
#' @param sigma A number. The noise variance added on the diagonal (single
#'   output case only).
#'
#' @return
#' A `tibble` containing the simulated data for the single task. Includes
#' columns: `Task_ID`, `Input_ID`, `Input`, `Output_ID`, and `Output`.
#'
#' @keywords internal
generate_single_task_data <- function(
    task_id,
    mean_process_info,
    task_hp_tibble = NULL,
    n_points_per_output,
    shared_grid_outputs = FALSE,
    precomputed_task_grid_list = NULL,
    v = NULL,
    l = NULL,
    sigma = NULL
) {
  num_outputs <- length(mean_process_info$grid_list)

  ## 1. Sample a sparse sub-grid for this task (same strategy for any number
  ## of outputs).
  if (is.null(precomputed_task_grid_list)) {
    if (shared_grid_outputs) {
      shared_task_grid <- sample_grid_points(
        grid_obj = mean_process_info$grid_list[[1]],
        size = n_points_per_output[1]
      )
      task_grid_list <- rep(list(shared_task_grid), num_outputs)
    } else {
      task_grid_list <- purrr::map2(mean_process_info$grid_list,
                                    n_points_per_output,
                                    sample_grid_points)
    }
  } else {
    if (length(precomputed_task_grid_list) != num_outputs) {
      stop("'precomputed_task_grid_list' length must match the number of",
           " outputs.")
    }
    task_grid_list <- precomputed_task_grid_list
  }

  ## 2. Build the task input data.frame and its reference names.
  input_df <- purrr::imap_dfr(
    task_grid_list,
    ~coerce_grid_to_df(.x) %>% dplyr::mutate(Output_ID = as.factor(.y))
  )
  reference <- build_reference_names(input_df)

  ## 3. Compute the task-specific covariance matrix (kernel-choice fork).
  if (num_outputs == 1) {
    K_task_t <- kern_to_cov(
      input_df %>% dplyr::select(-.data$Output_ID),
      "SE",
      tibble::tibble(se_variance = v, se_lengthscale = l)
    ) + diag(sigma, nrow(input_df))
  } else {
    if (is.null(task_hp_tibble)) {
      stop("'task_hp_tibble' is required when there is more than one output.")
    }
    K_task_t <- kern_to_cov(
      input = input_df,
      kern = convolution_kernel_KD,
      hp = task_hp_tibble
    )
  }
  rownames(K_task_t) <- colnames(K_task_t) <- reference

  ## 4. Select the corresponding mean process subset and draw the task data.
  mu0_subset <- mean_process_info$mean_process_realization[reference]
  y_task <- MASS::mvrnorm(n = 1, mu = mu0_subset, Sigma = K_task_t)

  ## 5. Format the result, keeping track of the hyperparameters used so that
  ## they can optionally be exposed as columns (see 'add_hp' in simu_data()).
  result <- input_df %>%
    dplyr::mutate(
      Task_ID = factor(task_id),
      Input_ID = factor(1),
      Output = as.numeric(y_task)
    )

  if (num_outputs == 1) {
    result <- result %>%
      dplyr::mutate(se_variance = v, se_lengthscale = l, noise = sigma)
  } else {
    result <- result %>%
      dplyr::left_join(
        task_hp_tibble %>% dplyr::select(-dplyr::any_of("Task_ID")),
        by = "Output_ID"
      )
  }

  result %>%
    dplyr::select(.data$Task_ID, .data$Input_ID, dplyr::everything())
}



#' Simulate a dataset tailored for MagmaClustR
#'
#' Simulate a complete synthetic dataset, which may be representative of various
#' applications. Several flexible arguments allow adjustment of the number of
#' tasks, clusters, outputs, or dimension of observed inputs, and the values of
#' many parameters controlling the data generation.
#'
#' @param n_tasks An integer. The number of tasks (per cluster).
#' @param n_points An integer. The number of observations per task, for each
#'    output.
#' @param n_clusters An integer. The number of underlying clusters.
#' @param n_outputs An integer. The number of outputs.
#'    If n_outputs = 1, a standard SE kernel is used.
#'    If n_outputs > 1, the Multi-Output Convolution kernel structure is simulated.
#' @param input2 A logical value indicating whether the dataset should
#'    include an additional input named 'Input2'.
#' @param prior_means A vector of numbers. Prior means for each output's mean
#'    process. If `NULL` (default), a random linear trend is generated for
#'    each output (see `m0_intercept`/`m0_slope`).
#' @param grid_input A vector of numbers defining a grid of observations
#'    (i.e. the reference inputs).
#' @param grid_input2 A vector of numbers defining a grid of observations
#'    (i.e. the reference inputs).
#' @param shared_input A logical value indicating whether the reference inputs
#'    are shared across tasks.
#' @param shared_hp  A logical value indicating whether the hyper-parameters are
#'   shared across tasks. If TRUE and n_points >1, the hyper-parameters remain
#'   different between the clusters.
#' @param add_hp A logical value indicating whether the values of
#'    hyper-parameters should be added as columns in the dataset.
#' @param add_clust A logical value indicating whether the name of the
#'    clusters should be added as a column in the dataset.
#' @param int_mu_v A vector of 2 numbers. Interval for mean process variance.
#' @param int_mu_l A vector of 2 numbers. Interval for mean process lengthscale.
#' @param int_i_v A vector of 2 numbers. Interval for task process variance.
#' @param int_i_l A vector of 2 numbers. Interval for task process lengthscale.
#' @param int_i_sigma A vector of 2 numbers. Interval for noise hyper-parameter.
#' @param m0_slope A vector of 2 numbers. Slope interval for the mean
#'    process's default linear prior trend.
#' @param m0_intercept A vector of 2 numbers. Intercept interval for the mean
#'    process's default linear prior trend.
#' @param warning_grid A boolean, indicating whether data generation should be
#'    stopped by anticipation if the size of input grid is > 2000 points. To
#'    force data generation to complete, switch this argument to \code{FALSE}
#'
#' @return A full dataset of simulated training data.
#'
#' @export
#'
#' @examples
#'
#'
#' data = simu_data()
#'
#' data_MO = simu_data(n_outputs = 2)
#'
#' data_clust = simu_data(n_clusters = 3)
#'
#' data_2D = simu_data(
#'    input2 = TRUE,
#'    grid_input = seq(-1, 1, 0.1),
#'    grid_input2 = seq(-1, 1, 0.1))
simu_data <- function(
    n_tasks = 10,
    n_points = 10,
    n_clusters = 1,
    n_outputs = 1,
    input2 = FALSE,
    prior_means = NULL,
    grid_input = seq(-1, 1, 0.01),
    grid_input2 = seq(-1, 1, 0.5),
    shared_input = FALSE,
    shared_hp = TRUE,
    add_hp = FALSE,
    add_clust = FALSE,
    int_mu_v = c(0, log(3)),
    int_mu_l = c(log(1/300), log(1/200)),
    int_i_v = c(0, log(1.5)),
    int_i_l = c(log(1/300), log(1/200)),
    int_i_sigma = c(0, 0.2),
    m0_slope = c(-5, 5),
    m0_intercept = c(-50, 50),
    warning_grid = TRUE
) {

  ## Stop if the size of the grid is too large
  if ((grid_n_points(grid_input) > 2000) && warning_grid) {
    stop("The defined input grid is > 2000 points. There is a high risk ",
        "that data generation takes a long time to complete, so the process ",
        "has been preemptively stopped. If you understand the issue, and still ",
        "want to force data generation to complete, use 'warning_grid=FALSE'.")
  }

  ## Build the domain over which each output's mean process grid is sampled.
  ## A list of per-output candidate grids (MO-specific feature) is preserved
  ## as-is; a single vector/data.frame is expanded with 'Input2' when relevant.
  is_per_output_pool <- is.list(grid_input) && !is.data.frame(grid_input)

  if (input2) {
    if (is_per_output_pool) {
      domain_grid <- purrr::map2(
        grid_input, grid_input2,
        ~tidyr::expand_grid(Input = .x, Input2 = .y)
      )
    } else {
      domain_grid <- tidyr::expand_grid(Input = grid_input, Input2 = grid_input2)

      if ((nrow(domain_grid) > 2000) && warning_grid) {
        stop("The defined input grid is > 2000 points. There is a high risk ",
            "that data generation takes a long time to complete, so the process ",
            "has been preemptively stopped. If you understand the issue, and still ",
            "want to force data generation to complete, use 'warning_grid=FALSE'.")
      }
    }
  } else {
    domain_grid <- grid_input
  }

  domain_size_per_output <- if (is_per_output_pool) {
    purrr::map_dbl(domain_grid, grid_n_points)
  } else {
    rep(grid_n_points(domain_grid), n_outputs)
  }

  ## Hyperparameter bound configuration (multi-output case only).
  hp_config_mean_process <- NULL
  hp_config_tasks <- NULL

  if (n_outputs > 1) {
    hp_config_mean_process <- tibble::tibble(
      output_id = 1:n_outputs,
      l0_min = rep(int_mu_l[1], n_outputs), l0_max = rep(int_mu_l[2], n_outputs),
      S0_min = rep(int_mu_v[1], n_outputs), S0_max = rep(int_mu_v[2], n_outputs),
      lu0_min = rep(int_mu_l[1], n_outputs), lu0_max = rep(int_mu_l[2], n_outputs)
    )

    hp_config_tasks <- tibble::tibble(
      output_id = 1:n_outputs,
      lt_min = rep(int_i_l[1], n_outputs), lt_max = rep(int_i_l[2], n_outputs),
      St_min = rep(int_i_v[1], n_outputs), St_max = rep(int_i_v[2], n_outputs),
      noise_min = rep(int_i_sigma[1], n_outputs),
      noise_max = rep(int_i_sigma[2], n_outputs),
      lu_min = rep(int_i_l[1], n_outputs), lu_max = rep(int_i_l[2], n_outputs)
    )
  }

  hp_cols_to_drop <- c(
    "se_variance", "se_lengthscale", "noise", "l_t", "S_t", "l_u_t"
  )

  floop_k <- function(k) {
    n_points_task <- rep(n_points, n_outputs)

    ## --- Mean process realization for this cluster (kernel-choice fork) ---
    if (n_outputs == 1) {
      mp_info_k <- generate_mean_process(
        points_per_output = domain_size_per_output,
        grid_input = domain_grid,
        prior_means = prior_means,
        m0_intercept = m0_intercept,
        m0_slope = m0_slope,
        v = draw(int_mu_v),
        l = draw(int_mu_l)
      )
    } else {
      mp_info_k <- generate_mean_process(
        points_per_output = domain_size_per_output,
        grid_input = domain_grid,
        prior_means = prior_means,
        m0_intercept = m0_intercept,
        m0_slope = m0_slope,
        hp_config_mean_process = hp_config_mean_process,
        shared_hp_outputs = FALSE
      )
    }

    ## --- Task hyperparameters shared across tasks within this cluster ---
    if (shared_hp) {
      if (n_outputs == 1) {
        shared_task_hp <- list(v = draw(int_i_v), l = draw(int_i_l),
                                sigma = draw(int_i_sigma))
      } else {
        shared_task_hp_tibble <- hp(
          kern = convolution_kernel_KD,
          list_task_ID = as.factor(1),
          list_output_ID = as.factor(1:n_outputs),
          shared_hp_tasks = TRUE,
          shared_hp_outputs = FALSE,
          noise = TRUE,
          hp_config = hp_config_tasks
        )
      }
    }

    ## When inputs are shared across tasks, the same sub-grid is sampled once
    ## per cluster and reused by every task; otherwise, each task samples its
    ## own sub-grid independently (handled inside generate_single_task_data).
    shared_task_grid_list <- if (shared_input) {
      purrr::map2(mp_info_k$grid_list, n_points_task, sample_grid_points)
    } else {
      NULL
    }

    floop_i <- function(i) {
      task_id <- make_task_id(i, k, n_clusters)

      if (n_outputs == 1) {
        task_hp <- if (shared_hp) {
          shared_task_hp
        } else {
          list(v = draw(int_i_v), l = draw(int_i_l), sigma = draw(int_i_sigma))
        }

        generate_single_task_data(
          task_id = task_id,
          mean_process_info = mp_info_k,
          n_points_per_output = n_points_task,
          precomputed_task_grid_list = shared_task_grid_list,
          v = task_hp$v,
          l = task_hp$l,
          sigma = task_hp$sigma
        )
      } else {
        task_hp_tibble <- if (shared_hp) {
          shared_task_hp_tibble
        } else {
          hp(
            kern = convolution_kernel_KD,
            list_task_ID = as.factor(1),
            list_output_ID = as.factor(1:n_outputs),
            shared_hp_tasks = TRUE,
            shared_hp_outputs = FALSE,
            noise = TRUE,
            hp_config = hp_config_tasks
          )
        }

        generate_single_task_data(
          task_id = task_id,
          mean_process_info = mp_info_k,
          task_hp_tibble = task_hp_tibble,
          n_points_per_output = n_points_task,
          precomputed_task_grid_list = shared_task_grid_list
        )
      }
    }

    cluster_db <- purrr::map_dfr(seq_len(n_tasks), floop_i)

    postprocess_cluster_db(
      cluster_db,
      add_hp = add_hp,
      add_clust = add_clust,
      k = k,
      hp_cols = hp_cols_to_drop
    )
  }

  final_db <- purrr::map_dfr(seq_len(n_clusters), floop_k)

  pivot_inputs_longer(final_db)
}

#' @rdname simu_data
#' @export
simu_db <- simu_data
