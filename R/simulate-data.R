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


#' @title Generate a Mean Process Realization
#' @description Creates a working grid for each output and generates a
#'              single smooth underlying mean process, `mu_0`, by drawing a
#'              realization from a Gaussian Process. This GP is centered around
#'              a simple linear function.
#'
#' @param points_per_output A numeric vector specifying the number of points for
#'                          each output's working grid.
#' @param grid_ranges A list of numeric vectors, where each defines the
#'                    `c(min, max)` input domain for the corresponding output.
#' @param hp_config_mean_process A tibble configuring the hyperparameters for
#'                               the mean process's covariance kernel.
#' @param shared_grid_outputs A logical value. If `TRUE`, all outputs are
#'                            defined on the exact same input grid. This
#'                            requires `points_per_output` and `grid_ranges`
#'                            to be identical for all outputs.
#' @param shared_hp_outputs A logicial value. If `TRUE`, all outputs share the
#'                          same hyperparameters.
#'
#' @return A list containing `mean_process_realization` (the `mu_0` values),
#'         `grid_list`, and `points_per_output`.
#'
generate_mean_process <- function(
    points_per_output,
    grid_ranges,
    hp_config_mean_process,
    shared_grid_outputs,
    shared_hp_outputs
) {
  num_outputs <- length(points_per_output)

  if (length(grid_ranges) != num_outputs) {
    stop("The number of grid ranges must match the number of outputs specified
         in 'points_per_output'.")
  }

  # STEP 1: Generate a dense grid for each output
  if (shared_grid_outputs) {
    # Validation: For a shared grid, the number of points and the ranges
    # must be identical for all outputs.
    if (length(unique(points_per_output)) > 1) {
      stop("For a shared grid, all elements of 'points_per_output' must be identical.")
    }
    # Trick: compare ranges by converting them to character strings
    if (length(unique(lapply(grid_ranges, as.character))) > 1) {
      stop("For a shared grid, all elements of 'grid_ranges' must be identical.")
    }

    # Generate the grid once using the configuration of the first output
    shared_grid <- sort(runif(n = points_per_output[1],
                              min = grid_ranges[[1]][1],
                              max = grid_ranges[[1]][2]))

    # Replicate this grid for each output in the list
    grid_list <- rep(list(shared_grid), num_outputs)

  } else {
    grid_list <- purrr::map2(
      points_per_output,
      grid_ranges,
      ~ sort(runif(n = .x, min = .y[1], max = .y[2]))
    )
  }
  X_unlist <- unlist(grid_list)

  # STEP 2: Define the GP's prior mean function, m_0(.)
  a <- runif(1, -0.5, 2)
  b <- runif(1, 0, 10)
  m0_mean_function <- a * X_unlist + b

  # STEP 3: Draw hyper-parameters for the mean process's GP covariance
  if (shared_hp_outputs) {
    # CASE 1: Shared HPs. Draw ONCE using the config of the first output.
    s0_shared <- runif(1, hp_config_mean_process$s0_min[1], hp_config_mean_process$s0_max[1])
    l0_shared <- runif(1, hp_config_mean_process$l0_min[1], hp_config_mean_process$l0_max[1])

    # Repeat the shared values for all outputs
    s0_variances <- rep(s0_shared, num_outputs)
    l0_lengthscales <- rep(l0_shared, num_outputs)

  } else {
    # CASE 2: Non-shared HPs. Draw one for each output
    s0_variances <- purrr::map2_dbl(
      hp_config_mean_process$s0_min,
      hp_config_mean_process$s0_max,
      ~runif(1, .x, .y)
    )
    l0_lengthscales <- purrr::map2_dbl(
      hp_config_mean_process$l0_min,
      hp_config_mean_process$l0_max,
      ~runif(1, .x, .y)
    )
  }

  # STEP 4: Build covariance and draw the mean process realization, mu_0
  cov_block_list <- lapply(1:num_outputs, function(i) {
    hp_i <- tibble::tibble("se_variance" = s0_variances[i],
                           "se_lengthscale" = l0_lengthscales[i])
    kern_to_cov(input = grid_list[[i]], kern = "SE", hp = hp_i)
  })
  K_theta0_X <- as.matrix(Matrix::bdiag(cov_block_list))

  prefix <- "pt_"
  all_point_names <- unlist(lapply(1:num_outputs,
                                   function(i) paste0("o", i, "_", prefix, grid_list[[i]])))
  rownames(K_theta0_X) <- colnames(K_theta0_X) <- all_point_names

  # Draw the realization from the GP
  mean_process_realization <- MASS::mvrnorm(n = 1,
                                            mu = m0_mean_function,
                                            Sigma = K_theta0_X)
  names(mean_process_realization) <- all_point_names

  return(list(
    mean_process_realization = mean_process_realization,
    grid_list = grid_list,
    points_per_output = points_per_output,
    list_s0 = s0_variances,
    list_l0 = l0_lengthscales
  ))
}


#' @title Generate Data for a Single Task
#' @description Simulates observed data points for one task. This version
#'              expects hyperparameters to be provided as a tibble.
#'
#' @param task_id A unique identifier for the task.
#' @param mean_process_info A list object from `generate_mean_process()`.
#' @param task_hp_tibble A tibble of hyperparameters for this specific task,
#'                       with one row per output. Must contain columns `l_t`,
#'                       `S_t`, `l_u_t`, and `noise`.
#' @param n_points_range A numeric vector `c(min, max)` for the number of points.
#'
#' @return A tibble containing the simulated data for the single task.

generate_single_task_data <- function(
    task_id,
    mean_process_info,
    task_hp_tibble,
    n_points_range
) {
  num_outputs <- length(mean_process_info$grid_list)

  # --- Sample a sparse sub-grid for this task ---
  n_points_per_output_task <- sample(n_points_range[1]:n_points_range[2],
                                     size = num_outputs, replace = TRUE)
  task_grid_list <- purrr::map2(mean_process_info$grid_list,
                                n_points_per_output_task,
                                ~sample(x = .x, size = .y) %>% sort())
  task_grid_list <- purrr::map(task_grid_list, as.matrix)

  # --- Build task-specific covariance and add noise ---

  # 1. Create the unique data.frame input from the sampled grid list
  input_df <- purrr::imap_dfr(task_grid_list,
                              ~dplyr::tibble(Input = .x[,1], Output_ID = .y)
  )

  # 2. Call kern_to_cov, passing the received hp_tibble directly
  K_task_t <- kern_to_cov(
    input = input_df,
    kern = convolution_kernel,
    hp = task_hp_tibble
  )

  # --- Select mean process subset and generate final data ---
  target_point_names <- unlist(lapply(1:num_outputs,
                                      function(i) paste0("o", i, "_pt_", task_grid_list[[i]])))
  mu0_subset <- mean_process_info$mean_process_realization[target_point_names]
  y_task <- MASS::mvrnorm(n = 1, mu = mu0_subset, Sigma = K_task_t)

  # --- Format into a clean tibble ---
  tibble::tibble(
    "Task_ID"   = factor(task_id),
    "Input_ID"  = factor(1),
    "Input"     = unlist(task_grid_list),
    "Output_ID" = factor(rep(1:num_outputs, times = n_points_per_output_task)),
    "Output"    = as.numeric(y_task)
  )
}


#' @title Simulate a Multi-Output, Multi-Task Dataset (Final Version)
#' @description Main orchestrator function that generates a complete dataset.
#'
#' @param num_tasks The total number of tasks to simulate.
#' @param points_per_output_grid A numeric vector for points per output's dense grid.
#' @param grid_ranges A list of c(min, max) vectors for each output's input domain.
#' @param hp_config_mean_process A tibble to configure the mean process's GP hyperparameters.
#' @param hp_config_tasks A tibble to configure the task-specific hyperparameters.
#' @param n_points_per_task_range A c(min, max) vector for points per task.
#' @param shared_hp_tasks If TRUE, all tasks share the same hyperparameter values.
#' @param shared_hp_outputs If TRUE, all outputs share the same HPs, both for the mean process and for the tasks.
#' @param shared_grid_outputs A logical value. If TRUE, all outputs are defined on the exact same input grid.
#' @param seed An optional integer for reproducibility.
#'
#' @return A list containing `simulated_data_df` and `mean_process_df`.

simulate_multi_output_data <- function(
    num_tasks = 50,
    points_per_output_grid = c(500, 150),
    grid_ranges = list(c(0, 10), c(0, 10)),
    hp_config_mean_process = tibble::tibble(
      output_id = 1:2, l0_min = c(-2, -2), l0_max = c(-2, -2),
      s0_min = c(1, 1), s0_max = c(1, 1)
    ),
    hp_config_tasks = tibble::tibble(
      output_id = 1:2, lt_min = c(-2, -2), lt_max = c(-2, -2),
      St_min = c(1, 1), St_max = c(1, 1)
    ),
    n_points_per_task_range = c(5, 20),
    shared_hp_tasks = FALSE,
    shared_hp_outputs = TRUE, # <<< ARGUMENT UNIQUE
    shared_grid_outputs = FALSE,
    seed = 456
) {
  if (!is.null(seed)) { set.seed(seed) }

  # === STEP 1: generate the mean process ===
  mean_process_info <- generate_mean_process(
    points_per_output = points_per_output_grid,
    grid_ranges = grid_ranges,
    hp_config_mean_process = hp_config_mean_process,
    shared_grid_outputs = shared_grid_outputs,
    shared_hp_outputs = shared_hp_outputs
  )

  # === STEP 2: generate HPs for all tasks simultaneously ===
  cat("Generating hyperparameters for all tasks...\n")
  all_tasks_hps <- hp(
    kern = convolution_kernel,
    list_task_ID = 1:num_tasks,
    list_output_ID = 1:nrow(hp_config_tasks),
    shared_hp_tasks = shared_hp_tasks,
    shared_hp_outputs = shared_hp_outputs,
    noise = TRUE
  )

  # === STEP 3: generate observed data for all tasks ===
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

  # === STEP 4: combined results for all tasks and format the final tibble ===
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
    mean_process_df = mean_process_df
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
#' @param seed A random seed for reproducible splitting.
#'
#' @return A list containing two dataframes: `train` and `test`.

split_data_train_test <- function(df_to_split,
                                  zones_definition = NULL,
                                  total_test_prop = 0.25,
                                  train_zone_leak_prop = 0.10,
                                  seed = 123) {

  set.seed(seed)

  # --- METHOD 1: Proportional split when no zones are defined ---
  if (is.null(zones_definition)) {

    cat("--- Mode: Simple Proportional Split (no zones defined) ---\n")

    n_total <- nrow(df_to_split)
    n_test <- round(n_total * total_test_prop)

    test_indices <- sample(1:n_total, size = n_test, replace = FALSE)

    df_test <- df_to_split[test_indices, ]
    df_train <- df_to_split[-test_indices, ]

    # --- METHOD 2: Zone-based split when zones are provided ---
  } else {

    cat("--- Mode: Zone-Based Split ---\n")

    df_with_zones <- df_to_split %>%
      rowwise() %>%
      mutate(
        Zone = {
          output_id_str <- as.character(Output_ID)
          point_intervals <- zones_definition[[output_id_str]]
          is_in_test <- any(map_lgl(point_intervals$test,
                                    ~Input >= .x[1] && Input < .x[2]))
          is_in_train <- any(map_lgl(point_intervals$train,
                                     ~Input >= .x[1] && Input <= .x[2]))
          if (is_in_test) "TestZone" else if (is_in_train) "TrainZone" else "None"
        }
      ) %>%
      ungroup()

    pool_train_zone <- df_with_zones %>% filter(Zone == "TrainZone")
    pool_test_zone <- df_with_zones %>% filter(Zone == "TestZone")

    n_to_leak <- round(nrow(pool_train_zone) * train_zone_leak_prop)
    leaked_to_test <- pool_train_zone %>%
      sample_n(size = n_to_leak, replace = FALSE)

    df_test <- bind_rows(pool_test_zone, leaked_to_test) %>%
      dplyr::select(-Zone)

    df_train <- anti_join(pool_train_zone,
                          leaked_to_test,
                          by = names(df_to_split)) %>%
      dplyr::select(-Zone)
  }

  # --- Return the list of dataframes ---
  return(list(train = df_train, test = df_test))
}


#' @title Initialize the Inverse Prior Covariance Matrix for the Mean Process
#'
#' @description
#' Calculates the inverse prior covariance matrix (K₀⁻¹) for the mean
#' Gaussian process (μ₀) in a multi-outputs model.
#' The function assumes that the different outputs are independent a priori,
#' which results in a sparse, block-diagonal matrix.
#' The block-diagonal structure of the returned matrix is a direct
#' consequence of the assumption of output independence in the prior. Each block
#' on the diagonal corresponds to the inverse covariance matrix for a specific
#' output, calculated over the set of all unique input points from `db`.
#'
#' @param db A `tibble` or `data.frame` containing the training data.
#'   Must include the columns `Task_ID`, `Output`, `Output_ID`, `Reference`,
#'   as well as the input coordinates (e.g., `Input_1`, `Input_2`, ...).
#' @param kern_0 The kernel to use for the mean process μ₀. Can be a
#'   character string (e.g., `"SE"`) or a kernel function object.
#' @param hp An object (e.g., a `tibble` or named vector) containing the
#'   hyperparameters required by `kern_0`.
#' @param pen_diag A numeric value. A "jitter" term added to the diagonal
#'   of the covariance matrix before its inversion to ensure numerical stability.
#'
#' @return A sparse, block-diagonal `Matrix` object, representing the
#'   inverse of the prior covariance matrix for the mean process.
#'
#' @export
ini_inverse_prior_cov <- function(db,
                                  kern_0 = "SE",
                                  hp,
                                  pen_diag) {
  # browser()
  ## Compute the prior inverse covariance for the mean process (mu_0)
  # Get the union of all unique input points from the training data
  all_inputs <- db %>%
    dplyr::select(-c(Task_ID, Output)) %>%
    unique() %>%
    dplyr::arrange(Reference)

  # Compute the inverse covariance matrix for each output block of the mean process
  # This assumes the prior on mu_0 treats outputs as independent GPs.
  list_inv_0 <- list_outputs_blocks_to_inv(
    db = db,
    kern = kern_0,
    hp = hp, # Note: Ensure that hp_0 in your original code corresponds to the hp argument.
    pen_diag = pen_diag
  )

  # Create the full block-diagonal inverse covariance matrix for mu_0
  inv_0 <- Matrix::bdiag(list_inv_0)

  # Set the row and column names of inv_0
  all_references <- unlist(lapply(list_inv_0, rownames), use.names = FALSE)
  dimnames(inv_0) <- list(all_references, all_references)
  inv_0 <- as.matrix(inv_0)
  return(inv_0)
}
