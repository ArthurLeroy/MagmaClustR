#' Simulate a batch of data
#'
#' Simulate a batch of output data, corresponding to one individual, coming from
#' a GP with a the Squared Exponential kernel as covariance structure, and
#' specified hyper-parameters and input.
#'
#' @param ID An identification code, whether numeric or character.
#' @param input A vector of numbers. The input variable that is used as
#'  'reference' for input and outputs.
#' @param mean A vector of numbers. Prior mean values of the GP.
#' @param v A number. The variance hyper-parameter of the SE kernel.
#' @param l A number. The lengthscale hyper-parameter of the SE kernel.
#' @param sigma A number. The noise hyper-parameter.
#'
#' @return A tibble containing a batch of output data along with input and
#'  additional information for a simulated individual.
#'
#' @importFrom rlang .data
#'
#' @keywords internal
#'
#' @examples
#' TRUE
simu_indiv_se <- function(ID, input, mean, v, l, sigma) {
  if(is.data.frame(input)){
    db <- tibble::tibble(
      "ID" = ID,
      input,
      "Output" = mvtnorm::rmvnorm(
        1,
        mean,
        kern_to_cov(
          input,
          "SE",
          tibble::tibble(se_variance = v, se_lengthscale = l)
        ) +
          diag(sigma, length(input$Reference), length(input$Reference))
      ) %>% as.vector(),
      "se_variance" = v,
      "se_lengthscale" = l,
      "noise" = sigma
    )
  }else{
    db <- tibble::tibble(
      "ID" = ID,
      "Input" = input,
      "Output" = mvtnorm::rmvnorm(
        1,
        mean,
        kern_to_cov(
          input,
          "SE",
          tibble::tibble(se_variance = v, se_lengthscale = l)
        ) +
          diag(sigma, length(input), length(input))
      ) %>% as.vector(),
      "se_variance" = v,
      "se_lengthscale" = l,
      "noise" = sigma
    )
  }
  return(db)
}

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
  stats::runif(1, int[1], int[2]) %>%
    round(2) %>%
    return()
}


#' Simulate a dataset tailored for MagmaClustR
#'
#' Simulate a complete training dataset, which may be representative of various
#' applications. Several flexible arguments allow adjustment of the number of
#' individuals, of observed inputs, and the values of many parameters
#' controlling the data generation.
#'
#' @param M An integer. The number of individual per cluster.
#' @param N An integer. The number of observations per individual.
#' @param K An integer. The number of underlying clusters.
#' @param covariate A logical value indicating whether the dataset should
#'    include an additional input covariate named 'Covariate'.
#' @param grid A vector of numbers defining a grid of observations
#'    (i.e. the reference inputs).
#' @param grid_cov A vector of numbers defining a grid of observations
#'    (i.e. the covariate reference inputs).
#' @param common_input A logical value indicating whether the reference inputs
#'    are common to all individual.
#' @param common_hp  A logical value indicating whether the hyper-parameters are
#'   common to all individual. If TRUE and K>1, the hyper-parameters remain
#'   different between the clusters.
#' @param add_hp A logical value indicating whether the values of
#'    hyper-parameters should be added as columns in the dataset.
#' @param add_clust A logical value indicating whether the name of the
#'    clusters should be added as a column in the dataset.
#' @param int_mu_v A vector of 2 numbers, defining an interval of admissible
#'    values for the variance hyper-parameter of the mean process' kernel.
#' @param int_mu_l A vector of 2 numbers, defining an interval of admissible
#'    values for the lengthscale hyper-parameter of the mean process' kernel.
#' @param int_i_v A vector of 2 numbers, defining an interval of admissible
#'    values for the variance hyper-parameter of the individual process' kernel.
#' @param int_i_l A vector of 2 numbers, defining an interval of admissible
#'    values for the lengthscale hyper-parameter of the individual process'
#'    kernel.
#' @param int_i_sigma A vector of 2 numbers, defining an interval of admissible
#'    values for the noise hyper-parameter.
#' @param lambda_int A vector of 2 numbers, defining an interval of admissible
#'    values for the lambda parameter of the 2D exponential.
#' @param m_int A vector of 2 numbers, defining an interval of admissible
#'    values for the mean of the 2D exponential.
#' @param lengthscale_int A vector of 2 numbers, defining an interval of
#'    admissible values for the lengthscale parameter of the 2D exponential.
#' @param m0_intercept A vector of 2 numbers, defining an interval of admissible
#'    values for the intercept of m0.
#' @param m0_slope A vector of 2 numbers, defining an interval of admissible
#'    values for the slope of m0.
#'
#' @return A full dataset of simulated training data.
#' @export
#'
#' @examples
#' ## Generate a dataset with 3 clusters of 4 individuals, observed at 10 inputs
#' data = simu_db(M = 4, N = 10, K = 3)
#'
#' ## Generate a 2-D dataset with an additional input 'Covariate'
#' data = simu_db(covariate = TRUE)
#'
#' ## Generate a dataset where input locations are different among individuals
#' data = simu_db(common_input = FALSE)
#'
#' ## Generate a dataset with an additional column indicating the true clusters
#' data = simu_db(K = 3, add_clust = TRUE)
simu_db <- function(M = 10,
                    N = 10,
                    K = 1,
                    covariate = FALSE,
                    grid = seq(0, 10, 0.05),
                    grid_cov = seq(0, 10, 0.5),
                    common_input = TRUE,
                    common_hp = TRUE,
                    add_hp = FALSE,
                    add_clust = FALSE,
                    int_mu_v = c(4, 5),
                    int_mu_l = c(0, 1),
                    int_i_v = c(1, 2),
                    int_i_l = c(0, 1),
                    int_i_sigma = c(0, 0.2),
                    lambda_int = c(30,40),
                    m_int = c(0,10),
                    lengthscale_int = c(30, 40),
                    m0_slope = c(-5, 5),
                    m0_intercept = c(-50, 50)) {
  if(covariate){
    exponential_mean <- function(x, y, lambda, m1, m2, lengthscale){
      return (lambda * base::exp(-((x-m1)^2 + (y-m2)^2)/lengthscale))
    }
    t_0_tibble <- tidyr::expand_grid(Input = grid, Covariate = grid_cov) %>%
        purrr::modify(signif) %>%
        tidyr::unite("Reference", sep = ":", remove = FALSE) %>%
        dplyr::arrange(.data$Reference)

    t_0 <- t_0_tibble %>% dplyr::pull(.data$Reference)

    if (common_input) {
        t_i_input <- sample(grid, N, replace = F)
        t_i_covariate <- sample(grid_cov, N, replace = F)
        t_i_tibble <- tibble::tibble(Input = t_i_input,
                                     Covariate = t_i_covariate) %>%
          purrr::modify(signif) %>%
          tidyr::unite("Reference", sep = ":", remove = FALSE) %>%
          dplyr::arrange(.data$Reference)

      t_i <- t_i_tibble %>% dplyr::pull(.data$Reference)

    }

    floop_k <- function(k) {

      m_0 <- tibble::tibble(
        Input = t_0_tibble$Input,
        Covariate = t_0_tibble$Covariate) %>%
        purrr::modify(signif) %>%
        tidyr::unite("Reference", sep = ":", remove = FALSE) %>%
        dplyr::mutate(Output =  exponential_mean(.data$Input,
                                          .data$Covariate,
                                          lambda = draw(lambda_int),
                                          m1 = draw(m_int),
                                          m2 = draw(m_int),
                                          lengthscale = draw(lengthscale_int)))

        if (common_hp) {
          i_v <- draw(int_i_v)
          i_l <- draw(int_i_l)
          i_sigma <- draw(int_i_sigma)
        }

        floop_i <- function(i, k) {
          if (!common_input) {
            t_i_input <- sample(grid, N, replace = F)
            t_i_covariate <- sample(grid_cov, N, replace = F)
            t_i_tibble <- tibble::tibble(Input = t_i_input,
                                         Covariate = t_i_covariate) %>%
              purrr::modify(signif) %>%
              tidyr::unite("Reference", sep = ":", remove = FALSE) %>%
              dplyr::arrange(.data$Reference)

            t_i <- t_i_tibble %>% dplyr::pull(.data$Reference)
          }
          if (!common_hp) {
            i_v <- draw(int_i_v)
            i_l <- draw(int_i_l)
            i_sigma <- draw(int_i_sigma)
          }
          mean_i <- m_0 %>%
            dplyr::filter(.data$Reference %in% t_i) %>%
            dplyr::pull(.data$Output)

          if (K > 1) {
            ID <- paste0("ID", i, "-Clust", k)
          } else {
            ID <- as.character(i)
          }

          simu_indiv_se(
            ID = ID,
            input = t_i_tibble,
            mean = mean_i,
            v = i_v,
            l = i_l,
            sigma = i_sigma
          ) %>%
            return()
        }

        db <- sapply(seq_len(M), floop_i, k = k,simplify = F, USE.NAMES = T) %>%
          dplyr::bind_rows() %>%
          dplyr::select(-.data$Reference)
        if (!add_hp) {
          db <- db %>%
            dplyr::select(-c("se_variance", "se_lengthscale", "noise"))
        }

        if (add_clust) {
          db <- db %>% dplyr::mutate("Cluster" = k)
        }
        return(db)
      }

      lapply(seq_len(K), floop_k) %>%
        dplyr::bind_rows() %>%
        return()
  }else{
    if (common_input) {
      t_i <- sample(grid, N, replace = F) %>% sort()
    }

    floop_k <- function(k) {
      m_0 <- draw(m0_intercept) + draw(m0_slope) * grid
      mu_v <- draw(int_mu_v)
      mu_l <- draw(int_mu_l)

      db_0 <- simu_indiv_se(
        ID = "0",
        input = grid,
        mean = m_0,
        v = mu_v,
        l = mu_l,
        sigma = 0
      )
      if (common_hp) {
        i_v <- draw(int_i_v)
        i_l <- draw(int_i_l)
        i_sigma <- draw(int_i_sigma)
      }

      floop_i <- function(i, k) {
        if (!common_input) {
          t_i <- sample(grid, N, replace = F) %>% sort()
        }
        if (!common_hp) {
          i_v <- draw(int_i_v)
          i_l <- draw(int_i_l)
          i_sigma <- draw(int_i_sigma)
        }
        mean_i <- db_0 %>%
          dplyr::filter(.data$Input %in% t_i) %>%
          dplyr::pull(.data$Output)

        if (K > 1) {
          ID <- paste0("ID", i, "-Clust", k)
        } else {
          ID <- as.character(i)
        }

        simu_indiv_se(
          ID = ID,
          input = t_i,
          mean = mean_i,
          v = i_v,
          l = i_l,
          sigma = i_sigma
        ) %>%
          return()
      }
      db <- sapply(seq_len(M), floop_i, k = k, simplify = F, USE.NAMES = T) %>%
        dplyr::bind_rows()

      if (!add_hp) {
        db <- db %>% dplyr::select(-c("se_variance", "se_lengthscale", "noise"))
      }

      if (add_clust) {
        db <- db %>% dplyr::mutate("Cluster" = k)
      }
      return(db)
    }

    lapply(seq_len(K), floop_k) %>%
      dplyr::bind_rows() %>%
      return()
  }
}


#' @title Generate a Mean Process Realization (Convolution Kernel Version)

#' @description Creates a working grid for each output and generates a single
#'    smooth underlying mean process, `mu_0`, by drawing a
#'    realization from a multi-output Gaussian Process using a
#'    convolution kernel. This GP is centered around a simple linear function.
#'
#' @param points_per_output A numeric vector specifying the number of points for
#'    each output's working grid.
#' @param grid_ranges A list of numeric vectors, where each defines the
#'    `c(min, max)` input domain for the corresponding output.
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
#'
#' @return A list containing the realization, grid, and hyperparameters.
#'
generate_mean_process_convol <- function(
    points_per_output,
    grid_ranges,
    hp_config_mean_process,
    prior_means,
    noise_0 = NULL,
    shared_grid_outputs,
    shared_hp_outputs
) {
  # Extract the number of outputs
  num_outputs <- length(points_per_output)

  if (length(grid_ranges) != num_outputs) {
    stop("The number of grid ranges must match the number of outputs.")
  }

  if (shared_grid_outputs) {
    if (length(unique(points_per_output)) > 1 ||
        length(unique(lapply(grid_ranges, as.character))) > 1) {
      stop("For a shared grid, 'points_per_output' and 'grid_ranges' must be identical.")
    }
    shared_grid <- sort(runif(n = points_per_output[1],
                              min = grid_ranges[[1]][1],
                              max = grid_ranges[[1]][2]))
    grid_list <- rep(list(shared_grid), num_outputs)
  } else {
    grid_list <- purrr::map2(
      points_per_output,
      grid_ranges,
      ~ sort(runif(n = .x, min = .y[1], max = .y[2]))
    )
  }

  # Use pmap to iterate over the grids and prior means coefficients
  # simultaneously. This creates a list where each element is the mean vector
  # for the corresponding output.
  mean_list <- purrr::pmap(
    .l = list(grid = grid_list, prior_means = prior_means),
    .f = function(grid, prior_means) {
      grid + prior_means
    }
  )
  m0_mean_function <- unlist(mean_list)

  # Generate HPs from the lower and upper bounds given by the user
  lu0_shared <- runif(1,
                      hp_config_mean_process$lu0_min[1],
                      hp_config_mean_process$lu0_max[1])
  lu0_vals <- rep(lu0_shared, num_outputs)

  if (shared_hp_outputs) {
    l0_shared <- runif(1,
                       hp_config_mean_process$l0_min[1],
                       hp_config_mean_process$l0_max[1])
    S0_shared <- runif(1,
                       hp_config_mean_process$S0_min[1],
                       hp_config_mean_process$S0_max[1])

    l0_vals <- rep(l0_shared, num_outputs)
    S0_vals <- rep(S0_shared, num_outputs)
  } else {
    l0_vals <- purrr::map2_dbl(hp_config_mean_process$l0_min,
                               hp_config_mean_process$l0_max,
                               ~runif(1, .x, .y))
    S0_vals <- purrr::map2_dbl(hp_config_mean_process$S0_min,
                               hp_config_mean_process$S0_max,
                               ~runif(1, .x, .y))
  }

  # Create a complete Input dataframe, where each Output_ID is associated to
  # its own grid of inputs
  input_df_mean_process <- purrr::imap_dfr(grid_list,
                                           ~dplyr::tibble(Input = .x, Output_ID = .y)
  )

  # Create an HPs tibble for the mean process
  hp_tibble_for_kernel <- tibble::tibble(
    Output_ID = as.factor(1:num_outputs),
    l_t       = l0_vals,
    S_t       = S0_vals,
    l_u_t     = lu0_vals
  )

  # Call kern_to_cov() to builld the covariance matrix of the mean process
  # We use suppressWarnings() to avoid raising a warning on the fact that
  # kern_to_cov() is used with a set of HPs that does not contain noise
  K_theta0_X <- suppressWarnings(kern_to_cov(
    input = input_df_mean_process,
    kern = convolution_kernel,
    hp = hp_tibble_for_kernel
  ))


  # If we want to put some noise on the mean process.
  # WARNING: this SHOULD NOT be used in a Magma MO context, but is necessary
  # when the user wants to simulate MO data
  if (!is.null(noise_0)) {
    if (length(noise_0) != num_outputs) {
      stop(sprintf(
        "'noise_0' length should match with the number of outputs (%d).",
        length(noise_0), num_outputs
      ))
    }

    # Create a vector of noises for each point of the inputs grid
    noise_variances <- exp(noise_0[input_df_mean_process$Output_ID])

    # Add these noises to the covariance matrix diagonal
    diag(K_theta0_X) <- diag(K_theta0_X) + noise_variances
  }

  all_point_names <- paste0("o", input_df_mean_process$Output_ID,
                            "_pt_", input_df_mean_process$Input)
  rownames(K_theta0_X) <- colnames(K_theta0_X) <- all_point_names

  # Finally draw the mean process realisation
  mean_process_realization <- MASS::mvrnorm(
    n = 1,
    mu = m0_mean_function,
    Sigma = K_theta0_X
  )
  names(mean_process_realization) <- all_point_names

  # Format the result into a tibble
  mean_process_df <- tibble::tibble(
    Input = unlist(grid_list),
    Output_ID = factor(rep(
      1:length(points_per_output),
      times = points_per_output
    )),
    Output = as.numeric(mean_process_realization)
  )

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


#' @title Generate Data for a Single Task
#'
#' @description
#' Generates a simulated dataset for a single, specific task. This function
#' samples a sparse sub-grid from a provided mean process grid, computes a
#' task-specific covariance matrix (using `kern_to_cov`), and then draws a
#' realization from a multivariate normal distribution.=
#'
#' @param task_id An identifier (numeric or string) for the task being
#'   generated.
#' @param mean_process_info A list object containing information from the
#'   shared mean process (e.g., output from
#'   `generate_mean_process_convol`). Must contain at least:
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
#'
#' @return
#' A `tibble` containing the simulated data for the single task. Includes
#' columns: `Task_ID`, `Input_ID`, `Input`, `Output_ID`, and `Output`.


generate_single_task_data <- function(
    task_id,
    mean_process_info,
    task_hp_tibble,
    n_points_range
) {

  # Extract the number of outputs
  num_outputs <- length(mean_process_info$grid_list)

  # Sample a sparse sub-grid for this task.
  # Each task has the same number of observations
  n_points_per_output_task <- rep(n_points_range[2], num_outputs)

  task_grid_list <- purrr::map2(mean_process_info$grid_list,
                                n_points_per_output_task,
                                ~sample(x = .x, size = .y) %>% sort())
  task_grid_list <- purrr::map(task_grid_list, as.matrix)

  ## Build task-specific covariance and add noise
  # Create the unique data.frame input from the sampled grid list
  input_df <- purrr::imap_dfr(task_grid_list,
                              ~dplyr::tibble(Input = .x[,1], Output_ID = .y)
  )

  input_df$Output_ID <- as.factor(input_df$Output_ID)

  # Call kern_to_cov(), passing the received hp_tibble directly
  K_task_t <- kern_to_cov(
    input = input_df,
    kern = convolution_kernel,
    hp = task_hp_tibble
  )

  ## Select mean process subset and generate final data
  target_point_names <- unlist(lapply(1:num_outputs,
                                      function(i) paste0("o", i,
                                                         "_pt_", task_grid_list[[i]])))
  mu0_subset <- mean_process_info$mean_process_realization[target_point_names]

  # Draw a sample for task_id
  y_task <- MASS::mvrnorm(n = 1, mu = mu0_subset, Sigma = K_task_t)

  # Format result into a clean tibble
  tibble::tibble(
    "Task_ID"   = factor(task_id),
    "Input_ID"  = factor(1),
    "Input"     = unlist(task_grid_list),
    "Output_ID" = factor(rep(1:num_outputs, times = n_points_per_output_task)),
    "Output"    = as.numeric(y_task)
  ) %>% return()
}



#' @title Simulate a Multi-Output, Multi-Task Dataset
#' @description Main orchestrator function that generates a complete dataset.
#'
#' @param num_tasks Total number of tasks.
#' @param points_per_output_grid Points for dense grid per output.
#' @param grid_ranges Input domain for each output.
#' @param prior_means A vector containing the values of the prior mean of each
#'    output.
#' @param hp_config_mean_process Configuration tibble for mean process GP HPs.
#' @param hp_config_tasks Tibble for task-specific HPs.
#' @param n_points_per_task_range Min/max points per task.
#' @param shared_hp_tasks If TRUE, tasks share HPs.
#' @param shared_hp_outputs If TRUE, outputs share HPs.
#' @param shared_grid_outputs If TRUE, outputs share the same grid.
#'
#' @return A list with simulated data, mean process, and HPs.
#'
#' @example
#'
#' simulation_results <- simulate_multi_output_data_convol(
#' num_tasks = 30,
#' points_per_output_grid = c(200, 200),
#' grid_ranges = list(c(-1, 1), c(-1, 1)),
#' prior_means = c(0, 0),
#' hp_config_mean_process = tibble::tibble(
#'   output_id = 1:2,
#'   l0_min = c(log(1/600), log(1/200)), l0_max = c(log(1/600), log(1/200)),
#'   S0_min = c(log(1), log(5)), S0_max = c(log(1), log(5)),
#'   lu0_min = c(log(1/200), log(1/200)), lu0_max = c(log(1/200), log(1/200)),
#' ),
#'
#' hp_config_tasks = tibble::tibble(
#'   output_id = 1:2,
#'   lt_min = c(log(1/600), log(1/200)), lt_max = c(log(1/600), log(1/200)),
#'   St_min = c(log(1), log(5)), St_max = c(log(1), log(5)),
#'   noise_min = c(-3, -3), noise_max = c(-3, -3),
#'   lu_min = c(log(1/200), log(1/200)), lu_max = c(log(1/200), log(1/200))
#' ),
#' n_points_per_task_range = c(20, 20),
#' shared_hp_tasks = TRUE,
#' shared_hp_outputs = FALSE,
#' shared_grid_outputs = FALSE,
#' )


simulate_multi_output_data_convol <- function(
    num_tasks = 50,
    points_per_output_grid = c(500, 150),
    grid_ranges = list(c(0, 10), c(0, 10)),
    prior_means = c(0, 0),
    hp_config_mean_process = tibble::tibble(
      output_id = 1:2,
      l0_min = c(-3, -3), l0_max = c(3, 3),
      S0_min = c(-3, -3), S0_max = c(3, 3),
      lu0_min = c(-1, -1), lu0_max = c(3, 3)
    ),

    hp_config_tasks = tibble::tibble(
      output_id = 1:2,
      lt_min = c(-3, -3), lt_max = c(3, 3),
      St_min = c(-3, -3), St_max = c(3, 3),
      noise_min = c(-3, -3), noise_max = c(0, 0),
      lu_min = c(-1, -1), lu_max = c(3, 3)
    ),
    n_points_per_task_range = c(5, 20),
    shared_hp_tasks = FALSE,
    shared_hp_outputs = FALSE,
    shared_grid_outputs = FALSE
) {

  # Generate mean process with convolution kernel
  mean_process_info <- generate_mean_process_convol(
    points_per_output = points_per_output_grid,
    grid_ranges = grid_ranges,
    hp_config_mean_process = hp_config_mean_process,
    prior_means = prior_means,
    noise_0 = NULL,
    shared_grid_outputs = shared_grid_outputs,
    shared_hp_outputs = shared_hp_outputs
  )

  # Generate one realisation for each task of the dataset
  cat("Generating hyperparameters for all tasks...\n")
  all_tasks_hps <- hp(
    kern = convolution_kernel,
    list_task_ID = as.factor(1:num_tasks),
    list_output_ID = as.factor(1:nrow(hp_config_tasks)),
    shared_hp_tasks = shared_hp_tasks,
    shared_hp_outputs = shared_hp_outputs,
    noise = TRUE,
    hp_config = hp_config_tasks
  )

  task_dfs_list <- purrr::map(1:num_tasks, ~{
    current_task_hp <- all_tasks_hps %>%
      dplyr::filter(Task_ID == .x)
    generate_single_task_data(
      task_id = .x,
      mean_process_info = mean_process_info,
      task_hp_tibble = current_task_hp,
      n_points_range = n_points_per_task_range
    )
  })

  # Concatenate all realisations
  simulated_data_df <- dplyr::bind_rows(task_dfs_list)
  mean_process_df <- tibble::tibble(
    Input = unlist(mean_process_info$grid_list),
    Output_ID = factor(rep(
      1:length(mean_process_info$points_per_output),
      times = mean_process_info$points_per_output
    )),
    mean_value = as.numeric(mean_process_info$mean_process_realization)
  )

  cat(sprintf(
    "Simulation complete. Generated %d points for %d tasks.\n",
    nrow(simulated_data_df), num_tasks
  ))

  return(list(
    simulated_data_df = simulated_data_df,
    mean_process_df = mean_process_df,
    hyperparameters = list(
      mean_process = mean_process_info,
      tasks = all_tasks_hps
    )
  ))
}


#' Simulate a Multi-Output, Multi-Task, Clustered Dataset (Convolution Kernel)
#'
#' @description
#' Generates a synthetic dataset structured for MagmaClust models.
#' Returns a tibble of task hyperparameters that includes a 'Cluster_ID' column.
#'
#' @param nb_clusters Integer. The number of underlying clusters to generate.
#' @param tasks_per_cluster Integer. The number of tasks to simulate per cluster.
#' @param points_per_output_grid Numeric vector. Points for the dense grid per output.
#' @param grid_ranges List of numeric vectors (length 2). Input domain c(min, max).
#' @param prior_means Vector or List. Prior means for the mean processes.
#' @param shared_hp_clusts Logical. If TRUE, clusters share the same HP values.
#' @param hp_config_mean_process Tibble. HP configuration for mean processes.
#' @param hp_config_tasks Tibble. HP configuration for tasks.
#' @param n_points_per_task_range Numeric vector c(min, max).
#' @param shared_hp_tasks Logical. If TRUE, tasks share the same HP values.
#' @param shared_hp_outputs Logical. If TRUE, outputs share the same HP values.
#' @param shared_grid_outputs Logical. If TRUE, outputs share the same grid outputs.
#'
#' @return A list containing simulated data and hyperparameters (with Cluster_ID).
#' @export

simulate_magmaclust_data_convol <- function(
    nb_clusters = 3,
    tasks_per_cluster = 10,
    points_per_output_grid = c(500, 150),
    grid_ranges = list(c(0, 10), c(0, 10)),
    prior_means = c(0, 0),
    shared_hp_clusts = FALSE,
    hp_config_mean_process = tibble::tibble(
      output_id = 1:2,
      l0_min = c(-3, -3), l0_max = c(3, 3),
      S0_min = c(-3, -3), S0_max = c(3, 3),
      lu0_min = c(-1, -1), lu0_max = c(3, 3)
    ),
    hp_config_tasks = tibble::tibble(
      output_id = 1:2,
      lt_min = c(-3, -3), lt_max = c(3, 3),
      St_min = c(-3, -3), St_max = c(3, 3),
      noise_min = c(-3, -3), noise_max = c(0, 0),
      lu_min = c(-1, -1), lu_max = c(3, 3)
    ),
    n_points_per_task_range = c(5, 20),
    shared_hp_tasks = FALSE,
    shared_hp_outputs = FALSE,
    shared_grid_outputs = FALSE
) {

  num_outputs <- length(points_per_output_grid)

  # Control step
  expected_rows <- if(shared_hp_clusts) num_outputs else (num_outputs * nb_clusters)

  if(nrow(hp_config_mean_process) != expected_rows) stop("Dimension mismatch in hp_config_mean_process.")
  if(nrow(hp_config_tasks) != expected_rows) stop("Dimension mismatch in hp_config_tasks.")

  # Prior means
  list_prior_means <- list()
  if (is.list(prior_means)) {
    list_prior_means <- prior_means
  } else if (is.vector(prior_means)) {
    if (length(prior_means) == num_outputs * nb_clusters) {
      for (k in 1:nb_clusters) {
        indices_k <- seq(from = k, to = length(prior_means), by = nb_clusters)
        list_prior_means[[k]] <- prior_means[indices_k]
      }
    } else if (length(prior_means) == num_outputs) {
      for (k in 1:nb_clusters) list_prior_means[[k]] <- prior_means
    } else {
      stop("Invalid prior_means length.")
    }
  }

  # Hyper-parameters definition
  list_config_mp <- vector("list", nb_clusters)
  list_config_tasks <- vector("list", nb_clusters)

  if (shared_hp_clusts) {
    # HPs are shared among clusters
    paste("Hyper-parameters are supposed to be the same for each cluster within ",
          "each Output_ID.")

    fix_values <- function(config) {
      new_config <- config
      params <- names(config) %>% grep("_min$", ., value = TRUE) %>% gsub("_min$", "", .)
      for(p in params) {
        vals <- runif(nrow(config), config[[paste0(p, "_min")]], config[[paste0(p, "_max")]])
        new_config[[paste0(p, "_min")]] <- vals
        new_config[[paste0(p, "_max")]] <- vals
      }
      return(new_config)
    }
    fixed_mp_config <- fix_values(hp_config_mean_process)
    fixed_task_config <- fix_values(hp_config_tasks)

    for(k in 1:nb_clusters) {
      list_config_mp[[k]] <- fixed_mp_config
      list_config_tasks[[k]] <- fixed_task_config
    }
  } else {
    # HPs are distinct for each cluster
    paste("Hyper-parameters are supposed to be cluster-specific within each ",
          "Output_ID.")
    for (k in 1:nb_clusters) {
      indices_k <- seq(from = k, to = nrow(hp_config_mean_process), by = nb_clusters)
      list_config_mp[[k]] <- hp_config_mean_process %>% dplyr::slice(indices_k)
      list_config_tasks[[k]] <- hp_config_tasks %>% dplyr::slice(indices_k)
    }
  }

  # Loop over clusters
  cat(sprintf("Simulation started: %d Clusters x %d Outputs.\n", nb_clusters, num_outputs))

  cluster_results <- purrr::map(1:nb_clusters, function(k) {

    config_mp_k <- list_config_mp[[k]]
    config_tasks_k <- list_config_tasks[[k]]

    # Reinitialise HPs
    config_mp_k$output_id <- as.factor(1:num_outputs)
    config_tasks_k$output_id <- as.factor(1:num_outputs)

    # Generate mean processes (one for each cluster)
    mp_info_k <- generate_mean_process_convol(
      points_per_output = points_per_output_grid,
      grid_ranges = grid_ranges,
      hp_config_mean_process = config_mp_k,
      prior_means = list_prior_means[[k]],
      noise_0 = NULL,
      shared_grid_outputs = shared_grid_outputs,
      shared_hp_outputs = shared_hp_outputs
    )

    # Formate data of the mean process
    mp_df_k <- mp_info_k$mean_process_df %>%
      dplyr::mutate(Cluster_ID = paste0("K", factor(k)))

    mp_hps_k <- tibble::tibble(
      Cluster_ID = paste0("K", factor(k)),
      Output_ID = as.factor(1:num_outputs),
      l0  = mp_info_k$list_l0,
      S0  = mp_info_k$list_S0,
      lu0 = mp_info_k$list_lu0
    )

    # Generate Task Hyperparameters
    start_task_id <- (k - 1) * tasks_per_cluster + 1
    end_task_id   <- k * tasks_per_cluster
    task_ids_k    <- start_task_id:end_task_id

    tasks_hps_k <- hp(
      kern = convolution_kernel,
      list_task_ID = as.factor(task_ids_k),
      list_output_ID = as.factor(1:num_outputs),
      shared_hp_tasks = shared_hp_tasks,
      shared_hp_outputs = shared_hp_outputs,
      noise = TRUE,
      hp_config = config_tasks_k
    ) %>%
      dplyr::mutate(Cluster_ID = paste0("K",factor(k)))

    # Generate Task Data
    tasks_data_list <- purrr::map(task_ids_k, ~{
      current_task_hp <- tasks_hps_k %>% dplyr::filter(Task_ID == .x)
      generate_single_task_data(
        task_id = .x,
        mean_process_info = mp_info_k,
        task_hp_tibble = current_task_hp,
        n_points_range = n_points_per_task_range
      )
    })

    cluster_tasks_df <- dplyr::bind_rows(tasks_data_list) %>%
      dplyr::mutate(Cluster_ID = factor(k))

    return(list(
      data = cluster_tasks_df,
      mp_df = mp_df_k,
      mp_info = mp_info_k,
      mp_hps = mp_hps_k,
      tasks_hps = tasks_hps_k
    ))
  })

  # Aggregate
  cat("Simulation completed.\n")

  return(list(
    simulated_data_df = purrr::map_dfr(cluster_results, "data"),
    mean_processes_df = purrr::map_dfr(cluster_results, "mp_df"),
    hyperparameters = list(
      mean_process_hps = purrr::map_dfr(cluster_results, "mp_hps"),
      tasks_hps = purrr::map_dfr(cluster_results, "tasks_hps")
    )
  ))
}



#' @title Splits a dataframe into training and test sets.
#' @description This function splits a dataframe using one of two methods:
#'   1.  **Zone-based split:** If `zones_definition` is provided, it splits data
#'                             based on geographical zones,
#'                             ensuring no training points appear in test zones.
#'   2.  **Proportional split:** If `zones_definition` is NULL, it performs a
#'                               simple random split of the entire
#'                               dataset based on the `total_test_prop`.
#'
#' @param df_to_split The full dataframe to be split.
#' @param zones_definition A list defining train/test zones. If NULL, a simple
#'  proportional split is performed. Defaults to NULL.
#' @param total_test_prop Used only when `zones_definition` is NULL.
#'  The target proportion of the dataset to be in the test set.
#' @param train_zone_leak_prop Used only when `zones_definition` is provided.
#'  The proportion of points from "Train Zones" to be moved to the test set.
#'
#' @return A list containing two dataframes: `train` and `test`.

split_data_train_test <- function(df_to_split,
                                  zones_definition = NULL,
                                  total_test_prop = 0.25,
                                  train_zone_leak_prop = 0.10,
                                  seed = 123) {

  # METHOD 1: Proportional split when no zones are defined
  if (is.null(zones_definition)) {

    cat("Mode: Simple Proportional Split (no zones defined) \n")

    n_total <- nrow(df_to_split)
    n_test <- round(n_total * total_test_prop)

    test_indices <- sample(1:n_total, size = n_test, replace = FALSE)

    df_test <- df_to_split[test_indices, ]
    df_train <- df_to_split[-test_indices, ]

    # METHOD 2: Zone-based split when zones are provided
  } else {

    cat("Mode: Zone-Based Split\n")

    df_with_zones <- df_to_split %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        Zone = {
          output_id_str <- as.character(Output_ID)
          point_intervals <- zones_definition[[output_id_str]]
          is_in_test <- any(purrr::map_lgl(point_intervals$test,
                                    ~Input >= .x[1] && Input < .x[2]))
          is_in_train <- any(purrr::map_lgl(point_intervals$train,
                                     ~Input >= .x[1] && Input <= .x[2]))
          if (is_in_test) "TestZone" else if (is_in_train) "TrainZone" else "None"
        }
      ) %>%
      dplyr::ungroup()

    pool_train_zone <- df_with_zones %>%
      dplyr::filter(Zone == "TrainZone")
    pool_test_zone <- df_with_zones %>%
      dplyr::filter(Zone == "TestZone")

    n_to_leak <- round(nrow(pool_train_zone) * train_zone_leak_prop)
    leaked_to_test <- pool_train_zone %>%
      dplyr::sample_n(size = n_to_leak, replace = FALSE)

    df_test <- dplyr::bind_rows(pool_test_zone, leaked_to_test) %>%
      dplyr::select(-Zone)

    df_train <- dplyr::anti_join(pool_train_zone,
                          leaked_to_test,
                          by = names(df_to_split)) %>%
      dplyr::select(-Zone)
  }

  # Return the list of dataframes
  return(list(train = df_train, test = df_test))
}
