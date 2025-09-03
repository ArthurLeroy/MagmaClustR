#' E-Step of the EM algorithm
#'
#' Expectation step of the EM algorithm to compute the parameters of the
#' hyper-posterior Gaussian distribution of the mean process in Magma.
#'
#' @param db A tibble or data frame. Columns required: ID, Input, Output.
#'    Additional columns for covariates can be specified.
#' @param m_0 A vector, corresponding to the prior mean of the mean GP.
#' @param kern_0 A kernel function, associated with the mean GP.
#' @param kern_t A kernel function, associated with the task GPs.
#' @param hp_0 A tibble or data frame of hyper-parameters
#'    associated with \code{kern_0}.
#' @param hp_t A tibble or data frame of hyper-parameters
#'    associated with \code{kern_t}.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#' @return A named list, containing the elements \code{mean}, a tibble
#' containing the Input and associated Output of the hyper-posterior's mean
#' parameter, and \code{cov}, the hyper-posterior's covariance matrix.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
#'

e_step <- function(db,
                   m_0,
                   kern_0,
                   kern_t,
                   hp_0,
                   hp_t,
                   pen_diag) {
  # browser()
  # Get the union of all unique input points from the training data
  all_inputs <- db %>%
    dplyr::select(-c(Task_ID, Output)) %>%
    unique() %>%
    dplyr::arrange(Reference)

  ## Compute the inverse covariance of the mean process
  # inv_0 <- ini_inverse_prior_cov(db, kern_0, hp_0, pen_diag)
  inv_0 <- matrix(0, nrow = nrow(all_inputs), ncol = nrow(all_inputs)) %>%
    `rownames<-`(all_inputs$Reference) %>%
    `colnames<-`(all_inputs$Reference)

  list_inv_t <- list()
  list_ID_task <- unique(db$Task_ID)

  list_output_ID <-  db$Output_ID %>% unique()

  # For each task, compute its full multi-output inverse covariance matrix
  for (t in list_ID_task) {
    # browser()
    # Isolate the data and HPs for the current task
    db_t <- db %>% dplyr::filter(Task_ID == t) %>%
                   dplyr::select(-c(Output, Task_ID))
    hp_t_indiv <- hp_t %>% dplyr::filter(Task_ID == t)

    if(length(list_output_ID) > 1 && !(kern_t %>% is.character())){
      # MO case with dependent outputs
      # Call kern_to_cov directly.
      # It will handle the multi-output structure and the noise addition internally.
      # 'kern_t' is expected to be the 'convolution_kernel' function.
      K_task_t <- kern_to_cov(
        input = db_t,
        kern = kern_t,
        hp = hp_t_indiv
      )
    } else if (length(list_output_ID) > 1 && kern_t %>% is.character()){
      # MO case with independent outputs

      hp_t_indiv <- hp_t_indiv %>% dplyr::select(-c(Task_ID, Output_ID))

      inv_t <- kern_to_cov(
        input = db_t,
        kern = kern_t,
        hp = hp_t_indiv
      )
    } else{
      # Extract all_inputs to call kern_to_cov() on the single output case
      all_inputs_t <- db %>%
        dplyr::filter(Task_ID == t) %>%
        dplyr::select(-c(Task_ID, Output, Output_ID)) %>%
        unique() %>%
        dplyr::arrange(Reference)

      K_task_t <- kern_to_cov(
        input = all_inputs_t,
        kern = kern_t,
        hp = hp_t_indiv
      )
    }

    # Store the correct row/column names before they are lost during inversion
    task_references <- rownames(K_task_t)

    # Invert the covariance matrix (this strips the names)
    K_inv_t <- K_task_t %>% chol_inv_jitter(pen_diag = pen_diag)
    # paste0('Det de K_inv_t, avec t = ', t, ' : ', det(K_inv_t))

    # Re-apply the stored names to the inverted matrix
    dimnames(K_inv_t) <- list(task_references, task_references)

    # Add the inverted matrix to the list
    # The rownames are already correctly set by kern_to_cov
    list_inv_t[[t]] <- K_inv_t
  }

  ## Update the posterior distribution
  # Create a named list of output values, split by task
  list_output_t <- base::split(db$Output, list(db$Task_ID))

  ##--------------- Update Posterior Inverse Covariance ---------------##
  post_inv <- inv_0
  for (inv_t in list_inv_t) {
    # browser()
    # Find the common input points between the mean process and the current task
    co_input <- intersect(row.names(inv_t), row.names(post_inv))

    # Add the task's contribution to the posterior inverse covariance
    post_inv[co_input, co_input] <- post_inv[co_input, co_input] +
      inv_t[co_input, co_input]
  }

  post_cov <- post_inv %>%
    chol_inv_jitter(pen_diag = pen_diag) %>%
    `rownames<-`(all_inputs %>% dplyr::pull(.data$Reference)) %>%
    `colnames<-`(all_inputs %>% dplyr::pull(.data$Reference))

  ##--------------------- Update Posterior Mean ---------------------##
  weighted_0 <- inv_0 %*% m_0

  for (t in names(list_inv_t)) {
    # Compute the weighted mean for the t-th task
    weighted_t <- list_inv_t[[t]] %*% list_output_t[[t]]

    # Find the common input points between the mean process and the current task
    co_input <- intersect(row.names(weighted_t), row.names(weighted_0))

    # Add the task's contribution to the posterior weighted mean
    weighted_0[co_input, ] <- weighted_0[co_input, ] +
      weighted_t[co_input, ]
  }

  # Compute the final posterior mean
  post_mean <- post_cov %*% weighted_0 %>% as.vector()

  ## Part 4: Format and return the results
  # Format the posterior mean into a tibble
  tib_mean <- tibble::tibble(all_inputs,
                             "Output" = post_mean)

  return(list(
    "mean" = tib_mean,
    "cov" = post_cov
  ))
}


#' @title Maximisation Step of the EM Algorithm (Simplified)
#' @description Computes the optimal hyper-parameters of the kernels involved in the model.
#'
#' @param db A tibble or data frame. Required columns: Task_ID, Output_ID, Input_1, Output, Reference.
#' @param m_0 Prior mean vector for the mean GP.
#' @param kern_0 Kernel function for the mean GP.
#' @param kern_t Kernel function for the individual task GPs.
#' @param old_hp_0 Hyper-parameters for the mean GP from the previous step.
#' @param old_hp_t Hyper-parameters for the task GPs from the previous step.
#' @param post_mean Posterior mean from the E-step.
#' @param post_cov Posterior covariance from the E-step.
#' @param shared_hp_tasks If TRUE, all tasks share the same hyper-parameter values.
#' @param pen_diag A jitter term for matrix inversion.
#'
#' @return A tibble containing the updated hyperparameters`hp_t`.
#'

m_step <- function(db,
                   m_0,
                   kern_0,
                   kern_t,
                   old_hp_0,
                   old_hp_t,
                   post_mean,
                   post_cov,
                   shared_hp_tasks,
                   pen_diag) {
  # browser()

  list_hp_0 <- old_hp_0 %>% names()
  list_ID_task <- unique(db$Task_ID)
  output_ids_vector <- unique(db$Output_ID)

  if (length(output_ids_vector) > 1){
    # MO case
    # Optimise hyper-parameters of the mean process
    # new_hp_0 <- stats::optim(
    #   par = old_hp_0,
    #   fn = logL_GP_mod,
    #   gr = gr_GP_mod,
    #   db = post_mean,
    #   mean = m_0,
    #   kern = kern_0,
    #   post_cov = post_cov,
    #   pen_diag = pen_diag,
    #   hp_col_names = list_hp_0,
    #   output_ids = output_ids_vector,
    #   method = "L-BFGS-B",
    #   control = list(factr = 1e13, maxit = 25)
    # )$par %>%
    #   tibble::as_tibble_row()

    new_hp_0 <- old_hp_0

  } else {
    # Single output case
    # Optimise hyper-parameters of the mean process
    # old_hp_0 <- old_hp_0 %>%
    #   dplyr::select(-Output_ID) %>%
    #   as.numeric()
    #
    # list_hp_0 <- list_hp_0[list_hp_0 != "Output_ID"]
    #
    # new_hp_0 <- stats::optim(
    #   par = old_hp_0,
    #   fn = logL_GP_mod,
    #   gr = gr_GP_mod,
    #   db = post_mean,
    #   mean = m_0,
    #   kern = kern_0,
    #   post_cov = post_cov,
    #   pen_diag = pen_diag,
    #   hp_col_names = list_hp_0,
    #   output_ids = output_ids_vector,
    #   method = "L-BFGS-B",
    #   control = list(factr = 1e13, maxit = 25)
    # )$par %>% setNames(list_hp_0) %>%
    #   tibble::as_tibble_row()
    #
    # # To give the right format of HPs to the next E-step
    # new_hp_0$Output_ID <- "1"

    new_hp_0 <- old_hp_0
  }

  # =================================================================== #
  # Case 1: HPs are shared across tasks -> one optimisation over all data
  # =================================================================== #
  if (shared_hp_tasks) {
    cat("HPs are shared across tasks...\n")

    if (length(output_ids_vector) > 1){
      # Prepare parameters for optim() in MO case
      hp_per_output <- old_hp_t %>%
        dplyr::group_by(Output_ID) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Task_ID, -l_u_t) %>%
        tidyr::pivot_longer(cols = -Output_ID, names_to = "hp_name", values_to = "value") %>%
        dplyr::mutate(specific_name = paste(hp_name, Output_ID, sep = "_")) %>%
        dplyr::select(specific_name, value) %>%
        tibble::deframe()

      shared_hp_l_u_t <- old_hp_t$l_u_t[1]
      names(shared_hp_l_u_t) <- "l_u_t"
      par <- c(hp_per_output, shared_hp_l_u_t)
      hp_col_names <- names(par)
    } else {
      # Prepare parameters for optim() in single output case
      par <- old_hp_t %>%
        dplyr::select(-c(Task_ID, Output_ID)) %>%
        dplyr::slice(1)
      hp_col_names <- names(par)
    }

    # Optimise hyper-parameters of the task process
    result_optim <- stats::optim(
      par               = par,
      fn                = logL_GP_mod_shared_tasks,
      gr                = gr_GP_mod_shared_tasks,
      db                = db,
      mean              = post_mean,
      kern              = kern_t,
      post_cov          = post_cov,
      pen_diag          = pen_diag,
      hp_col_names      = hp_col_names,
      output_ids        = output_ids_vector,
      method            = "L-BFGS-B",
      control           = list(factr = 1e13, maxit = 25)
    )$par %>%
      tibble::as_tibble_row()

    # Reshape results
    if (length(output_ids_vector) > 1){
      # MO case
      new_hp_t <- result_optim %>%
        tidyr::pivot_longer(
          cols = -dplyr::any_of("l_u_t"),
          names_to = c(".value", "Output_ID"),
          names_pattern = "(.+)_(\\d+)$"
        ) %>%
        tidyr::crossing(Task_ID = list_ID_task, .)
    } else {
      # Single output case
      result_optim$Output_ID <- "1"
      result_optim <- result_optim %>%
        mutate(Task_ID = list(list_ID_task)) %>%
        unnest(Task_ID)
      new_hp_t <- result_optim
    }

    # =================================================================== #
    # Case 2: HPs are task-specific -> loop over tasks
    # =================================================================== #
  } else {
    cat("HPs are task-specific.\n")

    floop <- function(t) {
      cat(paste0("Optimizing for task ", t, "...\n"))
      # Filter data for the current task
      db_t <- db %>% dplyr::filter(Task_ID == t) %>%
                     dplyr::select(-Task_ID)

      input_t <- db_t %>% dplyr::pull(Reference)

      post_mean_t <- post_mean %>%
        dplyr::filter(Reference %in% input_t) %>%
        dplyr::pull(Output)

      post_cov_t <- post_cov[as.character(input_t), as.character(input_t)]

      old_hp_t_task <- old_hp_t %>% dplyr::filter(Task_ID == t)

      # Prepare parameters for optim() for this specific task
      if (length(output_ids_vector) > 1){
        hp_per_output <- old_hp_t_task %>%
          dplyr::select(-Task_ID, -l_u_t) %>%
          tidyr::pivot_longer(cols = -Output_ID, names_to = "hp_name", values_to = "value") %>%
          dplyr::mutate(specific_name = paste(hp_name, Output_ID, sep = "_")) %>%
          dplyr::select(specific_name, value) %>%
          tibble::deframe()

        shared_hp <- old_hp_t_task$l_u_t[1]
        names(shared_hp) <- "l_u_t"

        par_t <- c(hp_per_output, shared_hp)
        hp_col_names <- names(par_t)
      } else {
        ## Extract the hyper-parameters associated with the i-th individual
        par_t <- old_hp_t %>%
          dplyr::filter(Task_ID == t) %>%
          dplyr::select(-c(Task_ID, Output_ID))
        hp_col_names <- names(par_t)
      }

      # Optimisation for the single task 't'
      result_optim <- stats::optim(
        par               = par_t,
        fn                = logL_GP_mod,
        gr                = gr_GP_mod,
        db                = db_t,
        mean              = post_mean_t,
        kern              = kern_t,
        post_cov          = post_cov_t,
        pen_diag          = pen_diag,
        hp_col_names      = hp_col_names,
        output_ids        = output_ids_vector,
        method            = "L-BFGS-B",
        control           = list(factr = 1e13, maxit = 25)
      )$par %>%
        tibble::as_tibble_row()

      return(result_optim)
    }

    # browser()
    # Collect results from all tasks
    optim_results_by_task <- sapply(list_ID_task, floop, simplify = FALSE, USE.NAMES = TRUE) %>%
      tibble::enframe(name = "Task_ID") %>%
      tidyr::unnest(cols = value)

    # Reshape results
    if(length(output_ids_vector) > 1){
      # MO case
      new_hp_t <- optim_results_by_task %>%
        tidyr::pivot_longer(
          cols = -c(Task_ID, dplyr::any_of("l_u_t")),
          names_to = c(".value", "Output_ID"),
          names_pattern = "(.+)_(\\d+)$"
        )

      # --- Final standard formatting for the output tibble ---
      final_hp_names <- old_hp_t %>% dplyr::select(-Task_ID, -Output_ID) %>% names()
      new_hp_t <- new_hp_t %>%
        dplyr::select(Task_ID, Output_ID, dplyr::all_of(final_hp_names)) %>%
        dplyr::mutate(across(c(Task_ID, Output_ID), as.character)) %>%
        dplyr::arrange(as.numeric(Task_ID), as.numeric(Output_ID))
    } else {
      # Single output case
      optim_results_by_task$Output_ID = "1"
      new_hp_t <- optim_results_by_task
    }
  }

  list(
    "hp_0" = new_hp_0,
    "hp_t" = new_hp_t
  ) %>%
    return()
}
