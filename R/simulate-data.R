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


#' @title Generate a Mean Process Realization (Convolution Kernel Version)

#' Creates a working grid for each output and generates a single smooth
#' underlying mean process, `mu_0`, by drawing a realization from a multi-output
#' Gaussian Process using a convolution kernel.
#'
#' @param points_per_output A numeric vector specifying the number of points for
#'    each output's working grid.
#' @param grid_input A list of numeric vectors, where each defines the
#'    input domain for the corresponding output.
#' @param hp_config_mean_process A tibble configuring the hyperparameters for
#'    the mean process's convolution kernel.
#' @param prior_means A vector containing the values of the prior mean of each
#'    output.
#' @param noise_0 An optional numeric vector specifying the standard deviation
#'    of the noise for each output. If provided, its length must
#'    match the number of outputs. Default is `NULL` (no noise).
#' @param shared_grid_outputs A logical value. If `TRUE`, all outputs are
#'    defined on the exact same input grid.
#' @param shared_hp_outputs A logical value. If `TRUE`, all outputs share the
#'    same hyperparameters.
#' @param shared_grid_clusters A logical value. If `TRUE`, clusters share the same grid.
#' @param precomputed_grid_list An optional list of precomputed grids, one per output.
#' @param v A number. The variance hyper-parameter of the SE kernel.
#' @param l A number. The lengthscale hyper-parameter of the SE kernel.
#' @param sigma A number. The noise hyper-parameter.
#'
#' @return A list containing the realization, grid, and hyperparameters.
#'
#' @keywords internal
#'

generate_mean_process <- function(
    points_per_output,
    grid_input,
    hp_config_mean_process,
    prior_means,
    noise_0 = NULL,
    shared_grid_outputs,
    shared_hp_outputs,
    shared_grid_clusters = FALSE,
    precomputed_grid_list = NULL,
    v = NULL,
    l = NULL,
    sigma = 0
) {
  # Extract the number of outputs
  num_outputs <- length(points_per_output)

  if (num_outputs == 1) {
    # Generate mean process
    ID <- "0"

    if (length(prior_means) == 1) {
      mean_vec <- rep(prior_means, grid_n_points(grid_input))
    } else {
      mean_vec <- prior_means
    }

    if (is.data.frame(grid_input)) {
      db <- tibble::tibble(
        "ID" = ID,
        grid_input,
        "Output" = mvtnorm::rmvnorm(
          1,
          mean = mean_vec,
          sigma = kern_to_cov(
            grid_input,
            "SE",
            tibble::tibble(se_variance = v, se_lengthscale = l)
          ) + diag(sigma, nrow(grid_input), nrow(grid_input))
        ) %>% as.vector(),
        "se_variance" = v,
        "se_lengthscale" = l,
        "noise" = sigma
      )
    } else {
      db <- tibble::tibble(
        "ID" = ID,
        "Input" = grid_input,
        "Output" = mvtnorm::rmvnorm(
          1,
          mean = mean_vec,
          sigma = kern_to_cov(
            grid_input,
            "SE",
            tibble::tibble(se_variance = v, se_lengthscale = l)
          ) + diag(sigma, length(grid_input), length(grid_input))
        ) %>% as.vector(),
        "se_variance" = v,
        "se_lengthscale" = l,
        "noise" = sigma
      )
    }

    return(list(
      mean_process_realization = db,
      grid_list = list(grid_input)))

  } else {
    # Multi-Output case
    precision <- 1e6
    # Detect if we have a unique dataframe or a list of dataframes
    is_multidim <- is.data.frame(grid_input) || (is.list(grid_input) &&
                                             all(purrr::map_lgl(grid_input, is.data.frame)))

    # Grid generation
    if (is_multidim) {
      if (!is.null(precomputed_grid_list)) {
        grid_list <- precomputed_grid_list
      } else if (shared_grid_outputs) {
        # First dataframe of the list
        grid_to_sample <- if (is.data.frame(grid_input)) grid_input else grid_input[[1]]
        sampled_idx <- sample(seq_len(nrow(grid_to_sample)),
                              size = points_per_output[1],
                              replace = FALSE)
        shared_grid_input <- grid_to_sample[sampled_idx, , drop = FALSE]
        shared_grid_input <- shared_grid_input [do.call(order, shared_grid_input ), , drop = FALSE]
        grid_list <- rep(list(shared_grid_input ), num_outputs)
      } else {
        if (is.data.frame(grid_input)) {
          # If only one dataframe is provided
          grid_list <- purrr::map(points_per_output, ~{
            sampled_idx <- sample(seq_len(nrow(grid_input)), size = .x, replace = FALSE)
            sub_grid <- grid[sampled_idx, , drop = FALSE]
            sub_grid[do.call(order, sub_grid), , drop = FALSE]
          })
        } else {
          # If a list of dataframes is provided (MO case)
          grid_list <- purrr::map2(grid_input, points_per_output, ~{
            sampled_idx <- sample(seq_len(nrow(.x)), size = .y, replace = FALSE)
            sub_grid <- .x[sampled_idx, , drop = FALSE]
            sub_grid[do.call(order, sub_grid), , drop = FALSE]
          })
        }
      }
    } else {
      if (!is.list(grid_input)) {
        grid_ranges <- rep(list(range(grid_input)), num_outputs)
      } else {
        grid_ranges <- purrr::map(grid_input, range)
      }

      if (length(grid_ranges) != num_outputs) {
        stop("The number of grid ranges must match the number of outputs.")
      }

      if (!is.null(precomputed_grid_list)) {
        if (length(precomputed_grid_list) != num_outputs) {
          stop("'precomputed_grid_list' length must match the number of outputs.")
        }
        grid_list <- precomputed_grid_list

      } else if (shared_grid_outputs) {
        if (length(unique(points_per_output)) > 1 ||
            length(unique(lapply(grid_ranges, as.character))) > 1) {
          stop("For a shared grid_input, 'points_per_output' and 'grid_ranges' must be identical.")
        }
        possible_points <- seq(from = grid_ranges[[1]][1],
                               to = grid_ranges[[1]][2],
                               by = 1/precision)
        shared_grid_input  <- sort(sample(possible_points,
                                   size = points_per_output[1],
                                   replace = FALSE))
        grid_list <- rep(list(shared_grid_input ), num_outputs)

      } else {
        grid_list <- purrr::map2(
          points_per_output,
          grid_ranges,
          function(n_pts, range) {
            possible_points <- seq(from = range[1], to = range[2], by = 1/precision)
            sort(sample(possible_points, size = n_pts, replace = FALSE))
          }
        )
      }
    }

    mean_list <- purrr::pmap(
      .l = list(grid_input  = grid_list, prior_means = prior_means),
      .f = function(grid_input, prior_means) {
        if (length(prior_means) == 1) {
          if (is.data.frame(grid_input)) {
            return(rep(prior_means, nrow(grid_input)))
          }
          return(rep(prior_means, length(grid_input)))
        }

        if (length(prior_means) != grid_n_points(grid_input)) {
          stop("Each element of 'prior_means' must be a scalar or match the grid_input  length.")
        }
        prior_means
      }
    )
    m0_mean_function <- unlist(mean_list)

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

    input_df_mean_process <- purrr::imap_dfr(grid_list,
                                             ~coerce_grid_to_df(.x) %>%
                                               dplyr::mutate(Output_ID = as.factor(.y))
    )

    hp_tibble_for_kernel <- tibble::tibble(
      Output_ID = as.factor(1:num_outputs),
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
          "'noise_0' length should match with the number of outputs (%d).",
          length(noise_0), num_outputs
        ))
      }
      noise_variances <- exp(noise_0[input_df_mean_process$Output_ID])
      diag(K_theta0_X) <- diag(K_theta0_X) + noise_variances
    }

    all_point_names <- rownames(K_theta0_X)

    mean_process_realization <- MASS::mvrnorm(
      n = 1,
      mu = m0_mean_function,
      Sigma = K_theta0_X
    )
    names(mean_process_realization) <- all_point_names

    mean_process_df <- input_df_mean_process %>%
      dplyr::mutate(Output = as.numeric(mean_process_realization))

    return(list(
      mean_process_realization = mean_process_realization,
      grid_list = grid_list,
      points_per_output = points_per_output,
      list_l0 = l0_vals,
      list_S0 = S0_vals,
      list_lu0 = lu0_vals,
      noise_0 = noise_0,
      mean_process_df = mean_process_df
    ))
  }
}


#' @title Generate Data for a Single Task

#' Generates a simulated dataset for a single, specific task, centered around
#' the mean process. This function samples a sparse sub-grid from a provided
#' mean process grid, computes a task-specific covariance matrix (using
#' `kern_to_cov`), and then draws a realization from a multivariate normal
#' distribution.
#'
#' @param task_id An identifier (numeric or string) for the task being
#'   generated.
#' @param mean_process_info A list object containing information from the
#'   shared mean process (e.g., output from
#'   `generate_mean_process`). Must contain at least:
#'   \itemize{
#'     \item `grid_list`: A list of input grids, one for each output.
#'     \item `mean_process_realization`: A named vector of the full mean
#'       process realization.
#'   }
#' @param task_hp_tibble A `tibble` containing the specific hyperparameters
#'   for this task. This tibble is passed directly to the `hp` argument
#'   of `kern_to_cov`.
#' @param n_points_range A numeric vector of length 2, `c(min, max)`.
#'   The function currently uses `n_points_range[2]` (the max value) as the
#'   fixed number of points to sample for each output.
#' @param shared_grid_outputs A logical value. If `TRUE`, all outputs in the
#'   task are observed on the same sampled inputs.
#' @param precomputed_task_grid_list An optional list of precomputed task grids,
#'   one per output. If provided, the same sampled inputs are reused instead of
#'   drawing a new sub-grid for the task.
#' @param input_i A vector or data.frame representing the inputs for the task.
#' @param mean_i A vector representing the mean process for the task.
#' @param v A number. The variance hyper-parameter of the SE kernel.
#' @param l A number. The lengthscale hyper-parameter of the SE kernel).
#' @param sigma A number. The noise hyper-parameter.
#'
#' @return
#' A `tibble` containing the simulated data for the single task. Includes
#' columns: `Task_ID`, `Input_ID`, `Input`, `Output_ID`, and `Output`.
#'
#' @keywords internal


generate_single_task_data <- function(
    task_id,
    mean_process_info,
    task_hp_tibble,
    n_points_range,
    shared_grid_outputs = FALSE,
    precomputed_task_grid_list = NULL,
    input_i = NULL,
    mean_i = NULL,
    v = NULL,
    l = NULL,
    sigma = NULL
) {
  # Extract the number of outputs
  num_outputs <- length(mean_process_info$grid_list)

  if (num_outputs == 1) {
    # Generate task data
    if (is.data.frame(input_i)) {
      db <- tibble::tibble(
        "ID" = task_id,
        input_i,
        "Output" = mvtnorm::rmvnorm(
          1,
          mean = mean_i,
          sigma = kern_to_cov(
            input_i,
            "SE",
            tibble::tibble(se_variance = v, se_lengthscale = l)
          ) + diag(sigma, nrow(input_i), nrow(input_i))
        ) %>% as.vector(),
        "se_variance" = v,
        "se_lengthscale" = l,
        "noise" = sigma
      )
    } else {
      db <- tibble::tibble(
        "ID" = task_id,
        "Input" = input_i,
        "Output" = mvtnorm::rmvnorm(
          1,
          mean = mean_i,
          sigma = kern_to_cov(
            input_i,
            "SE",
            tibble::tibble(se_variance = v, se_lengthscale = l)
          ) + diag(sigma, length(input_i), length(input_i))
        ) %>% as.vector(),
        "se_variance" = v,
        "se_lengthscale" = l,
        "noise" = sigma
      )
    }
    return(db)

  } else {
    # Sample a sparse sub-grid for this task.
    # Each task has the same number of observations
    n_points_per_output_task <- rep(n_points_range[2], num_outputs)

    if (is.null(precomputed_task_grid_list)) {
      if (shared_grid_outputs) {
        shared_task_grid <- sample_grid_points(
          grid_obj = mean_process_info$grid_list[[1]],
          size = n_points_per_output_task[1]
        )
        task_grid_list <- rep(list(shared_task_grid), num_outputs)
      } else {
        task_grid_list <- purrr::map2(mean_process_info$grid_list,
                                      n_points_per_output_task,
                                      ~sample_grid_points(grid_obj = .x, size = .y))
      }
    } else {
      if (length(precomputed_task_grid_list) != num_outputs) {
        stop("'precomputed_task_grid_list' length must match the number of",
             " outputs.")
      }
      task_grid_list <- precomputed_task_grid_list
    }

    ## Build task-specific covariance and add noise
    # Create the unique data frame input from the sampled grid list
    input_df <- purrr::imap_dfr(task_grid_list,
                                ~coerce_grid_to_df(.x) %>%
                                  dplyr::mutate(Output_ID = as.factor(.y)))

    # Call kern_to_cov(), passing the received hp_tibble directly
    K_task_t <- kern_to_cov(
      input = input_df,
      kern = convolution_kernel_KD,
      hp = task_hp_tibble
    )

    ## Select mean process subset and generate final data
    target_point_names <- rownames(K_task_t)
    mu0_subset <- mean_process_info$mean_process_realization[target_point_names]

    # Draw a sample for task_id
    y_task <- MASS::mvrnorm(n = 1, mu = mu0_subset, Sigma = K_task_t)

    # Format result into a clean tibble
    input_df %>%
      dplyr::mutate(
        Task_ID = factor(task_id),
        Input_ID = factor(1),
        Output = as.numeric(y_task)
      ) %>%
      dplyr::select(.data$Task_ID, .data$Input_ID, dplyr::everything()) %>%
      return()
  }
}



#' Simulate a dataset tailored for MagmaClustR
#'
#' Simulate a complete synthetic dataset, which may be representative of various
#' applications. Several flexible arguments allow adjustment of the number of
#' tasks, clusters, outputs, or dimension of observed inputs, and the values of
#' many parameters controlling the data generation.
#'
#' @param n_tasks An integer. The number of tasks (per cluster).
#' @param n_points An integer. The number of observations per task.
#' @param n_clusters An integer. The number of underlying clusters.
#' @param n_outputs An integer. The number of outputs.
#'    If n_outputs = 1, standard SE kernels are used.
#'    If n_outputs > 1, the Multi-Output Convolution kernel structure is simulated.
#' @param input2 A logical value indicating whether the dataset should
#'    include an additional input named 'Input2'.
#' @param prior_means A vector of numbers. Prior means for the MO processes
#'    (n_outputs > 1).
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
#' @param lambda_int A vector of 2 numbers. Lambda parameter of 2D exponential.
#' @param m_int A vector of 2 numbers. Mean of 2D exponential.
#' @param lengthscale_int A vector of 2 numbers. Lengthscale of 2D exponential.
#' @param m0_slope A vector of 2 numbers. Slope interval for m0.
#' @param m0_intercept A vector of 2 numbers. Intercept interval for m0.
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
    grid_input2 = seq(-1, 1, 0.1),
    shared_input = TRUE,
    shared_hp = TRUE,
    add_hp = FALSE,
    add_clust = FALSE,
    int_mu_v = c(0, log(3)),
    int_mu_l = c(log(1/300), log(1/200)),
    int_i_v = c(0, log(1.5)),
    int_i_l = c(log(1/300), log(1/200)),
    int_i_sigma = c(0, 0.2),
    lambda_int = c(30, 40),
    m_int = c(0, 10),
    lengthscale_int = c(30, 40),
    m0_slope = c(-5, 5),
    m0_intercept = c(-50, 50),
    warning_grid = TRUE
) {

  ## Stop if the size of the grid is too large
  if( (length(grid_input) > 2000) && warning_grid){
    stop("The defined input grid is > 2000 points. There is a high risk ",
        "that data generation takes a long time to complete, so the process ",
        "has been preemptively stopped. If you understand the issue, and still ",
        "want to force data generation to complete, use 'warning_grid=FALSE'.")
  }

  if (n_outputs == 1) {
    if (input2) {
      t_0_tibble <- tidyr::expand_grid(Input = grid_input , Input2 = grid_input2) %>%
        purrr::modify(signif) %>%
        tidyr::unite("Reference", sep = ":", remove = FALSE) %>%
        dplyr::arrange(.data$Reference)

      ## Stop if the size of the grid is too large
      if( (nrow(t_0_tibble) > 2000) && warning_grid){
        stop("The defined input grid is > 2000 points. There is a high risk ",
        "that data generation takes a long time to complete, so the process ",
        "has been preemptively stopped. If you understand the issue, and still ",
        "want to force data generation to complete, use 'warning_grid=FALSE'.")
      }

      t_0 <- t_0_tibble %>% dplyr::pull(.data$Reference)

      if (shared_input) {
        t_i_input <- sample(grid_input , n_points, replace = F)
        t_i_input2 <- sample(grid_input2, n_points, replace = F)
        t_i_tibble <- tibble::tibble(
          Input = t_i_input,
          Input2 = t_i_input2
        ) %>%
          purrr::modify(signif) %>%
          tidyr::unite("Reference", sep = ":", remove = FALSE) %>%
          dplyr::arrange(.data$Reference)

        t_i <- t_i_tibble %>% dplyr::pull(.data$Reference)
      }

      floop_k <- function(k) {
        mu_v <- draw(int_mu_v)
        mu_l <- draw(int_mu_l)

        db_0_info <- generate_mean_process(
          points_per_output = c(nrow(t_0_tibble)),
          grid_input = t_0_tibble %>% dplyr::select(.data$Input, .data$Input2),
          hp_config_mean_process = NULL,
          prior_means = draw(m_int),
          noise_0 = NULL,
          shared_grid_outputs= FALSE,
          shared_hp_outputs = FALSE,
          shared_grid_clusters = FALSE,
          precomputed_grid_list = NULL,
          v = mu_v,
          l = mu_l
        )

        m_0 <- db_0_info$mean_process_realization %>%
          dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), signif)) %>%
          tidyr::unite(
            "Reference",
            -dplyr::any_of(c("ID", "Output", "se_variance", "se_lengthscale", "noise")),
            sep = ":",
            remove = FALSE
          )

        if (shared_hp) {
          i_v <- draw(int_i_v)
          i_l <- draw(int_i_l)
          i_sigma <- draw(int_i_sigma)
        }

        floop_i <- function(i, k) {
          if (!shared_input) {
            t_i_input <- sample(grid_input , n_points, replace = F)
            t_i_input2 <- sample(grid_input2, n_points, replace = F)
            t_i_tibble <- tibble::tibble(
              Input = t_i_input,
              Input2 = t_i_input2
            ) %>%
              purrr::modify(signif) %>%
              tidyr::unite("Reference", sep = ":", remove = FALSE) %>%
              dplyr::arrange(.data$Reference)

            t_i <- t_i_tibble %>% dplyr::pull(.data$Reference)
          }
          if (!shared_hp) {
            i_v <- draw(int_i_v)
            i_l <- draw(int_i_l)
            i_sigma <- draw(int_i_sigma)
          }
          mean_i <- m_0 %>%
            dplyr::filter(.data$Reference %in% t_i) %>%
            dplyr::pull(.data$Output)

          if (n_clusters > 1) {
            ID <- paste0("ID", i, "-Clust", k)
          } else {
            ID <- as.character(i)
          }

          generate_single_task_data(
            task_id = ID,
            mean_process_info = db_0_info,
            task_hp_tibble = NULL,
            n_points_range = c(5, 20),
            shared_grid_outputs = FALSE,
            precomputed_task_grid_list = NULL,
            input_i = t_i_tibble %>% dplyr::select(.data$Input, .data$Input2),
            mean_i = mean_i,
            v = i_v,
            l = i_l,
            sigma = i_sigma
          ) %>% return()
        }

        db <- sapply(seq_len(n_tasks), floop_i, k = k, simplify = F, USE.NAMES = T) %>%
          dplyr::bind_rows()

        if (!add_hp) {
          db <- db %>%
            dplyr::select(-c("se_variance", "se_lengthscale", "noise"))
        }

        if (add_clust) {
          db <- db %>% dplyr::mutate("Cluster" = k)
        }
        return(db)
      }

      final_db <- lapply(seq_len(n_clusters), floop_k) %>% dplyr::bind_rows()

    } else {
      # Single Output & no input2
      if (shared_input) {
        t_i <- sample(grid_input , n_points, replace = F) %>% sort()
      }

      floop_k <- function(k) {
        m_0 <- draw(m0_intercept) + draw(m0_slope) * grid_input
        mu_v <- draw(int_mu_v)
        mu_l <- draw(int_mu_l)

        db_0 <- generate_mean_process(
          points_per_output = length(grid_input),
          grid_input = grid_input,
          hp_config_mean_process = NULL,
          prior_means = m_0,
          noise_0 = NULL,
          shared_grid_outputs= FALSE,
          shared_hp_outputs = FALSE,
          shared_grid_clusters = FALSE,
          precomputed_grid_list = NULL,
          v = mu_v,
          l = mu_l
        )

        if (shared_hp) {
          i_v <- draw(int_i_v)
          i_l <- draw(int_i_l)
          i_sigma <- draw(int_i_sigma)
        }

        floop_i <- function(i, k) {
          if (!shared_input) {
            t_i <- sample(grid_input, n_points, replace = F) %>% sort()
          }
          if (!shared_hp) {
            i_v <- draw(int_i_v)
            i_l <- draw(int_i_l)
            i_sigma <- draw(int_i_sigma)
          }
          mean_i <- db_0$mean_process_realization %>%
            dplyr::filter(.data$Input %in% t_i) %>%
            dplyr::pull(.data$Output)

          if (n_clusters > 1) {
            ID <- paste0("ID", i, "-Clust", k)
          } else {
            ID <- as.character(i)
          }

          generate_single_task_data(
            task_id = ID,
            mean_process_info = db_0,
            task_hp_tibble = NULL,
            n_points_range = c(5, 20),
            shared_grid_outputs = FALSE,
            precomputed_task_grid_list = NULL,
            input_i = t_i,
            mean_i = mean_i,
            v = i_v,
            l = i_l,
            sigma = i_sigma
          ) %>% return()
        }
        db <- sapply(seq_len(n_tasks), floop_i, k = k, simplify = F, USE.NAMES = T) %>%
          dplyr::bind_rows()

        if (!add_hp) {
          db <- db %>% dplyr::select(-c("se_variance", "se_lengthscale", "noise"))
        }

        if (add_clust) {
          db <- db %>% dplyr::mutate("Cluster" = k)
        }
        return(db)
      }

      final_db <- lapply(seq_len(n_clusters), floop_k) %>% dplyr::bind_rows()
    }

  } else {
    # Multi-Output case
    hp_config_mean_process <- tibble::tibble(
      output_id = 1:n_outputs ,
      l0_min = rep(int_mu_l[1], n_outputs),l0_max = rep(int_mu_l[2], n_outputs),
      S0_min = rep(int_mu_v[1], n_outputs),S0_max = rep(int_mu_v[2], n_outputs),
      lu0_min = rep(int_mu_l[1], n_outputs),lu0_max = rep(int_mu_l[2],n_outputs)
    )

    hp_config_tasks <- tibble::tibble(
      output_id = 1:n_outputs ,
      lt_min = rep(int_i_l[1], n_outputs), lt_max = rep(int_i_l[2], n_outputs),
      St_min = rep(int_i_v[1], n_outputs), St_max = rep(int_i_v[2], n_outputs),
      noise_min = rep(int_i_sigma[1], n_outputs),
      noise_max = rep(int_i_sigma[2], n_outputs),
      lu_min = rep(int_i_l[1], n_outputs), lu_max = rep(int_i_l[2], n_outputs)
    )

    if (!is.null(hp_config_mean_process) &&
        nrow(hp_config_mean_process) == n_outputs && n_clusters > 1) {
      hp_config_mean_process <- purrr::map_dfr(seq_len(n_clusters), ~hp_config_mean_process)
    }
    if (!is.null(hp_config_tasks) && nrow(hp_config_tasks) == n_outputs &&
        n_clusters > 1) {
      hp_config_tasks <- purrr::map_dfr(seq_len(n_clusters), ~hp_config_tasks)
    }

    list_config_mp <- vector("list", n_clusters)
    list_config_tasks <- vector("list", n_clusters)

    for (k in seq_len(n_clusters)) {
      start_idx <- (k - 1) * n_outputs + 1
      end_idx   <- k * n_outputs

      if (!is.null(hp_config_mean_process)) {
        list_config_mp[[k]] <- hp_config_mean_process %>%
          dplyr::slice(start_idx:end_idx)
      }
      if (!is.null(hp_config_tasks)) {
        list_config_tasks[[k]] <- hp_config_tasks %>%
          dplyr::slice(start_idx:end_idx)
      }
    }

    floop_k_MO <- function(k) {
      if (input2) {
        if (is.list(grid_input) && !is.data.frame(grid_input)) {
          working_grid_input <- purrr::map2(grid_input, grid_input2,
                                      ~tidyr::expand_grid(Input = .x, Input2 = .y))
        } else {
          working_grid_input <- tidyr::expand_grid(Input = grid_input, Input2 = grid_input2)
        }
      } else {
        working_grid_input <- grid_input
      }

      if (is.null(prior_means)) {
        prior_means_k <- purrr::map_dbl(1:n_outputs , ~ draw(m0_intercept))
      } else {
        prior_means_k <- prior_means
      }

      mp_info_k <- generate_mean_process(
        points_per_output = rep(200, n_outputs),
        grid_input = working_grid_input,
        hp_config_mean_process = list_config_mp[[k]],
        prior_means = prior_means_k,
        noise_0 = NULL,
        shared_grid_outputs= FALSE,
        shared_hp_outputs = FALSE,
        shared_grid_clusters = FALSE,
        precomputed_grid_list = NULL
      )

      floop_i_MO <- function(i, k) {
        if (n_clusters > 1) {
          task_id <- paste0("ID", i, "-Clust", k)
        } else {
          task_id <- as.character(i)
        }

        hp_i <- hp(
          kern = convolution_kernel_KD,
          list_task_ID = as.factor(1),
          list_output_ID = as.factor(1:n_outputs),
          shared_hp_tasks = TRUE,
          shared_hp_outputs = FALSE,
          noise = TRUE,
          hp_config = list_config_tasks[[k]]
        )

        generate_single_task_data(
          task_id = task_id,
          mean_process_info = mp_info_k,
          task_hp_tibble = hp_i,
          n_points_range = c(5, 20),
          shared_grid_outputs = FALSE,
          precomputed_task_grid_list = NULL
        )
      }

      cluster_data <- sapply(seq_len(n_tasks), floop_i_MO, k = k, simplify = FALSE, USE.NAMES = TRUE) %>%
        dplyr::bind_rows()

      if (add_clust) {
        cluster_data <- cluster_data %>% dplyr::mutate("Cluster" = k)
      }

      return(cluster_data)
    }

    final_db <- lapply(seq_len(n_clusters), floop_k_MO) %>% dplyr::bind_rows()

    if (!add_hp) {
      cols_to_remove <- c("l_t", "S_t", "l_u_t", "se_variance", "se_lengthscale", "noise")
      cols_to_remove <- intersect(cols_to_remove, names(final_db))
      if (length(cols_to_remove) > 0) {
        final_db <- final_db %>% dplyr::select(-dplyr::all_of(cols_to_remove))
      }
    }
  }

  if (input2) {
    final_db <- final_db %>%
      dplyr::rename(Input_1 = .data$Input, Input_2 = .data$Input2)
  } else {
    final_db <- final_db %>%
      dplyr::rename(Input_1 = .data$Input)
  }

  # Appel de la fonction de pivot universelle
  final_db <- format_longer(final_db)

  return(final_db)
}

#' @rdname simu_data
#' @export
simu_db <- simu_data
