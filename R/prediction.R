#' Gaussian Process prediction
#'
#' Compute the posterior distribution of a standard GP, using the formalism of
#' Magma. By providing observed data, the prior mean and covariance
#' matrix (by defining a kernel and its associated hyper-parameters), the mean
#' and covariance parameters of the posterior distribution are computed on the
#' grid of inputs that has been specified. This predictive distribution can be
#' evaluated on any arbitrary inputs since a GP is an infinite-dimensional
#' object.
#'
#' @param data A tibble or data frame. Required columns:
#'    \code{Input_ID}, \code{Input}, \code{Output_ID}, \code{Output}. Additional
#'    columns for covariates can be specified.
#'    The \code{Input_ID} contains the unique names/codes used to identify each
#'    explanatory variable.
#'    The \code{Input} column should define the value of the explanatory
#'    variables that are used as reference for the observations (e.g. time for
#'    longitudinal data).
#'    The \code{Output_ID} contains the unique names/codes used to identify each
#'    output (response) variable.
#'    The \code{Output} column specifies the observed values (the response
#'    variables).
#' @param grid_inputs The grid of inputs (reference Input and covariates) values
#'    on which the GP should be evaluated. Ideally, this argument should be a
#'    tibble or a data frame, providing the same columns as \code{data}, except
#'    'Output'. Nonetheless, in cases where \code{data} provides only one
#'    'Input' column, the \code{grid_inputs} argument can be NULL (default) or a
#'    vector. This vector would be used as reference input for prediction and if
#'    NULL, a vector of length 500 is defined, ranging between the min and max
#'    Input values of \code{data}.
#' @param mean Mean parameter of the GP. This argument can be specified under
#'    various formats, such as:
#'    - NULL (default). The mean would be set to 0 everywhere.
#'    - A number. The mean would be a constant function.
#' @param hp A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern}. The columns/elements should be named
#'    according to the hyper-parameters that are used in \code{kern}. If NULL
#'    (default), the function \code{\link{train_gp}} is called with random
#'    initial values for learning maximum-likelihood estimators of the
#'    hyper-parameters associated with \code{kern}.
#' @param kern A kernel function, defining the covariance structure of the GP.
#'    Several popular kernels
#'    (see \href{https://www.cs.toronto.edu/~duvenaud/cookbook/}{The Kernel
#'    Cookbook}) are already implemented and can be selected within the
#'    following list:
#'    - "SE": (default value) the Squared Exponential Kernel (also called
#'        Radial Basis Function or Gaussian kernel),
#'    - "LIN": the Linear kernel,
#'    - "PERIO": the Periodic kernel,
#'    - "RQ": the Rational Quadratic kernel.
#'    - convolution_kernel : the Convolution kernel used to manage MO scenario.
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#' @param get_full_cov A logical value, indicating whether the full posterior
#'    covariance matrix should be returned.
#' @param plot A logical value, indicating whether a plot of the results is
#'    automatically displayed.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#'
#' @return A list, containing:
#'    - a tibble, representing the GP predictions as two column 'Mean' and
#'   'Var', evaluated on the \code{grid_inputs}. The column 'Input' and
#'   additional covariates columns are associated to each predicted values.
#'    If the \code{get_full_cov} argument is TRUE, the function returns a list,
#'    in which the tibble described above is defined as 'pred' and the full
#'    posterior covariance matrix is defined as 'cov';
#'    - a matrix, corresponding to t(cov_crossed) %*% inv_obs
#' @export
#'
#' @examples
#' TRUE
pred_gp <- function(data = NULL,
                    grid_inputs = NULL,
                    mean = NULL,
                    hp = NULL,
                    kern = "SE",
                    get_full_cov = FALSE,
                    plot = TRUE,
                    pen_diag = 1e-10) {

  ## Create a dummy dataset if no data is provided
  if(data %>% is.null()){
    data = tibble::tibble(
      'ID' = 'ID_pred',
      'Input' = c(0,10),
      'Output' = c(0,0)
    )

    ## Draw random hyper-parameters if not provided
    if(hp %>% is.null()){
      hp = hp(kern = kern, noise = TRUE)
      cat(
        "The 'hp' argument has not been specified. Random hyper-parameters",
        "values have been drawn.\n \n"
      )
    }

    ## Create a boolean to remember no data were provided
    no_data = TRUE
  } else {
    no_data = FALSE
  }

  ## Remove the 'ID' column if present
  if ("Task_ID" %in% names(data)) {
    if (dplyr::n_distinct(data$Task_ID) > 1) {
      stop(
        "Problem in the 'Task_ID' column: different values are not allowed. ",
        "The prediction can only be performed for one task."
      )
    }
    data <- data %>%
      dplyr::select(-Task_ID)
  }

  ## Extract the list of different output IDs
  list_ID_output <- data$Output_ID %>% unique()

  data <- data %>%
    dplyr::group_by(Output_ID, Output, Input_ID) %>%
    # Add a unique number observation for the group
    dplyr::mutate(obs_num = row_number()) %>%
    dplyr::ungroup()

  ## To create the 'Reference' column as in the old MagmaClustR tibble format, we
  # need to pivot data to obtain one row per observation of (Task_ID, Output_ID).
  # In other words, inputs are no longer in "short" format; instead, we have one
  # column per input.
  data <- data %>%
    tidyr::pivot_wider(
      names_from = Input_ID,
      values_from = Input,
      names_prefix = "Input_"
    ) %>%
    # Keep 6 significant digits for Inputs to avoid numerical issues
    dplyr::mutate(across(starts_with("Input_"), ~ round(.x, 6))) %>%
    rowwise() %>%
    dplyr::mutate(
      Reference = paste(
        # Create output's prefix
        paste0("o", Output_ID),
        # Create the reference for each Output_ID
        paste(c_across(starts_with("Input_")), collapse = ":"),
        # Join output's prefix and reference
        sep = ";"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-obs_num)

  ## Get input column names
  if ("Reference" %in% names(data)) {
    names_col <- data %>%
      dplyr::select(- c(Output_ID, Output, Reference)) %>%
      names()
  } else {
    names_col <- data %>%
      dplyr::select(-c(Output_ID, Output)) %>%
      names()
  }

  ## Extract the observed Output (data points)
  data_obs <- data %>%
    dplyr::pull(Output)

  ## Name elements of data_obs with 'Reference' column
  names(data_obs) <- data$Reference

  ## Extract the observed inputs (reference Input + covariates)
  inputs_obs <- data %>%
    dplyr::select(-Output) %>%
    unique()

  ## Extract the observed (reference) Input
  input_obs <- inputs_obs %>% dplyr::pull(Reference)

  ## Define the target inputs to predict
  if (grid_inputs %>% is.null()) {
    set_grid <- function(data, size_grid) {
      seq(data %>% min(),
          data %>% max(),
          length.out = size_grid
      ) %>%
        return()
    }

    if (inputs_obs %>% names() %>% length() == 3) {
      size_grid <- 500
    } else if (inputs_obs %>% names() %>% length() > 3) {
      size_grid <- 1000^(1 / (ncol(inputs_obs) - 1)) %>% round()
      ## floor instead of round ?
    } else {
      stop(
        "The 'grid_inputs' argument should be a either a numerical vector ",
        "or a data frame depending on the context. Please read ?pred_magma()."
      )
    }

    # Create a unique grid of inputs (which will be replicated for each Output_ID)
    base_grid <- purrr::map(
      data %>% dplyr::select(tidyselect::all_of(names_col)),
      set_grid,
      size_grid
    ) %>%
      purrr::set_names(paste0("Input_", seq_along(.))) %>%
      expand.grid() %>%
      tibble::as_tibble()

    unique_outputs <- data %>%
      dplyr::distinct(Output_ID)

    # Cross base_grid with unique_outputs to replicate base_grid as many times as
    # unique_outputs and create 'Reference' column
    inputs_pred <- tidyr::crossing(unique_outputs, base_grid) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        Reference = paste(
          paste0("o", Output_ID),
          paste(c_across(starts_with("Input_")), collapse = ":"),
          sep = ";"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(Reference)

    input_pred <- inputs_pred %>% dplyr::pull(Reference)

  } else if (grid_inputs %>% is.vector()) {
    ## Test whether 'data' only provide the Input column and no covariates
    if (inputs_obs %>% names() %>% length() == 3) {
      input_temp <- grid_inputs %>%
        signif() %>%
        sort() %>%
        unique()

      inputs_pred <- tibble::tibble("Input_1" = rep(input_temp,
                                                    times = nrow(unique_outputs)),
                                    "Output_ID" = rep(unique_outputs %>%
                                                        purrr::as_vector(),
                                                      each = length(input_temp))) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          Reference = paste(
            paste0("o", Output_ID),
            paste(c_across(starts_with("Input_")), collapse = ":"),
            sep = ";"
          )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(Reference)

      input_pred <- inputs_pred %>% dplyr::pull(Reference)

    } else {
      stop(
        "The 'grid_inputs' argument should be a either a numerical vector ",
        "or a data frame depending on the context. Please read ?pred_magma()."
      )
    }
  } else if (grid_inputs %>% is.data.frame()) {
    grid_inputs <- grid_inputs %>%
      dplyr::group_by(Input_ID) %>%
      dplyr::mutate(id_ligne = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(
        names_from = Input_ID,
        values_from = Input,
        names_prefix = "Input_"
      ) %>%
      dplyr::select(-id_ligne) %>%
      # Keep 6 significant digits for Inputs to avoid numerical issues
      dplyr::mutate(across(starts_with("Input_"), ~ round(.x, 6))) %>%
      rowwise() %>%
      dplyr::mutate(
        Reference = paste(
          # Create output's prefix
          paste0("o", Output_ID),
          # Create the reference for each Output_ID
          paste(c_across(starts_with("Input_")), collapse = ":"),
          # Join output's prefix and reference
          sep = ";"
        )
      )

    ## Test whether 'data' has the same columns as grid_inputs
    if (names(inputs_obs) %>% setequal(names(grid_inputs))) {
      input_pred <- grid_inputs %>%
        # DO NOT arrange Reference because of the lexicographic order
        # (not armful but unnecessary in Magma case, tragic in MO case)
        dplyr::pull(Reference)

      inputs_pred <- grid_inputs %>%
        # DO NOT arrange Reference because of the lexicographic order
        # (not armful but unnecessary in Magma case, tragic in MO case)
        dplyr::select(names(inputs_obs))
    } else {
      stop(
        "The 'grid_inputs' argument should provide a column 'Input_ID', 'Input' ",
        " and 'Output_ID'."
      )
    }
  } else {
    stop(
      "The 'grid_inputs' argument should be a either a numerical vector ",
      "or a data frame depending on the context. Please read ?pred_gp()."
    )
  }

  if (mean %>% is.null()) {
    mean_obs <- rep(0, length(input_obs))
    mean_pred <- rep(0, length(input_pred))
    cat(
      "The 'prior_mean' argument has not been specified. The hyper_prior mean",
      "function is thus set to be 0 everywhere for all outputs.\n \n"
    )
  } else if (mean %>% is.vector()) {
    if (length(mean) == length(input_obs)) {
      mean_obs <- mean
      mean_pred <- mean
    } else if (length(mean) == length(data$Output_ID %>% unique())) {
      # Get the unique and sorted Output_IDs
      unique_outputs_sorted <- data$Output_ID %>% unique() %>% sort()

      # Create a lookup table: "o1" -> prior_mean[1], "o2" -> prior_mean[2], etc.
      # This assumes the prior_mean vector is provided in the sorted order of
      # Output_IDs.
      mean_map <- setNames(mean, paste0("o", unique_outputs_sorted))

      # Extract the prefix ("o1", "o2", etc.) from each element in all_input
      input_obs_prefixes <- stringr::str_extract(input_obs, "o[0-9]+")
      input_pred_prefixes <- stringr::str_extract(input_pred, "o[0-9]+")

      # Build m_0 using the lookup table; it will automatically repeat the correct
      # value for each prefix.
      mean_obs <- mean_map[input_obs_prefixes] %>% unname()
      mean_pred <- mean_map[input_pred_prefixes] %>% unname()

      cat(
        "A constant hyper_prior mean has been set for each output.\n \n"
      )
    } else {
      stop(
        "Incorrect format for the 'mean' argument. Please read ",
        "?pred_gp() for details."
      )
    }

  } else if (mean %>% is.function()) {
    mean <- mean(input_obs %>%
                   dplyr::select(-Reference)
    )
  } else {
    stop(
      "Incorrect format for the 'prior_mean' argument. Please read ",
      "?hyperposterior() for details."
    )
  }

  ## Learn the hyper-parameters if not provided
  if (hp %>% is.null()) {
    if (kern %>% is.function()) {
      stop(
        "When using a custom kernel function the 'hp' argument is ",
        "mandatory, in order to provide the name of the hyper-parameters. ",
        "You can use the function 'hp()' to easily generate a tibble of random",
        " hyper-parameters with the desired format, or use 'train_gp()' to ",
        "learn ML estimators for a better fit of data."
      )
    } else if (kern %>% is.character()) {
      hp <- quiet(
        train_gp(data,
                 prior_mean = mean_obs,
                 ini_hp = hp(kern, noise = T),
                 kern = kern,
                 hyperpost = NULL,
                 pen_diag = pen_diag
        )
      )

      cat(
        "The 'hp' argument has not been specified. The 'train_gp()' function",
        "(with random initialisation) has been used to learn ML estimators",
        "for the hyper-parameters associated with the 'kern' argument.\n \n"
      )
    } else {
      stop(
        "Incorrect format for the 'kern' argument. Please read ?pred_gp() for ",
        "details."
      )
    }
  }

  ## Remove the noise of the hp for evaluating some of the sub-matrix
  if ("noise" %in% names(hp)) {
    hp_rm_noi <- hp %>% dplyr::select(-noise)
    noise <- exp(hp[["noise"]])
  } else {
    hp_rm_noi <- hp
    noise <- 0
  }

  if(length(data$Output_ID %>% unique()) > 1){
    ## MO case
    ## Compute the required sub-matrix for prediction
    inv_obs <- kern_to_cov(inputs_obs, kern, hp) %>%
      chol_inv_jitter(pen_diag = pen_diag) %>%
      `rownames<-`(as.character(input_obs)) %>%
      `colnames<-`(as.character(input_obs))

    cov_pred <- kern_to_cov(inputs_pred, kern, hp_rm_noi)
    cov_crossed <- kern_to_cov(inputs_obs, kern, hp_rm_noi, input_2 = inputs_pred)
  } else {
    ## Single output case
    ## Compute the required sub-matrix for prediction
    inv_obs <- kern_to_cov(inputs_obs %>%
                             dplyr:: select(-Output_ID),
                           kern,
                           hp) %>%
      chol_inv_jitter(pen_diag = pen_diag) %>%
      `rownames<-`(as.character(input_obs)) %>%
      `colnames<-`(as.character(input_obs))

    cov_pred <- kern_to_cov(inputs_pred %>%
                              dplyr:: select(-Output_ID),
                            kern,
                            hp_rm_noi)

    cov_crossed <- kern_to_cov(inputs_obs %>%
                                 dplyr:: select(-Output_ID),
                               kern,
                               hp_rm_noi,
                               input_2 = inputs_pred %>%
                                 dplyr:: select(-Output_ID))
  }

  ## Order data_obs - mean_obs according to the column names of inv_obs
  obs <- data_obs - mean_obs
  obs <- obs[colnames(inv_obs)]

  ## Compute the posterior mean
  pred_mean <- (mean_pred +
                  t(cov_crossed) %*% inv_obs %*% obs) %>%
    as.vector()

  names(pred_mean) <- input_pred

  ## Compute the posterior covariance matrix
  pred_cov <- cov_pred - t(cov_crossed) %*% inv_obs %*% cov_crossed


  ## If no data were originally provided, return the GP prior as a prediction
  if(no_data){
    pred_mean = mean_pred
    pred_cov = cov_pred

    data = c()
  }

  if(length(data$Output_ID %>% unique()) > 1){
    ## MO case
    ## Get the pred_cov rownames and colnames to create a noise vector with correct
    ## dimensions
    point_names <- rownames(pred_cov)
    output_ids <- as.numeric(sub("^o(\\d+);.*", "\\1", point_names))
    full_noise_vector <- noise[output_ids]
  } else {
    ## Single output case
    full_noise_vector <- noise
  }


  ## Create a tibble of values and associated uncertainty from a GP prediction
  pred_gp <- tibble::tibble(
    "Mean" = pred_mean,
    "Var" = diag(pred_cov) + full_noise_vector
  ) %>%
    dplyr::mutate(inputs_pred) %>%
    dplyr::select(-Reference)

  ## Display the graph of the prediction if expected
  if (plot) {

    ## Display samples only in 1D and Credible Interval otherwise
    display_samples = dplyr::if_else(ncol(pred_gp) == 4, TRUE, FALSE)

    plot_gp(
      list('pred' = pred_gp, 'cov' = pred_cov),
      data = data,
      samples = display_samples
    ) %>% print()
  }

  ## Add the posterior covariance matrix in the results if expected
  if (get_full_cov) {
    pred_gp = list('pred' = pred_gp, 'cov' = pred_cov)
  }
  return(list("pred_gp" = pred_gp,
              "cov_crossed_times_inv_cov" = t(cov_crossed) %*% inv_obs
              )
  )
}



#' Compute the hyper-posterior distribution in Magma
#'
#' Compute the parameters of the hyper-posterior Gaussian distribution of the
#' mean process in Magma (similarly to the expectation step of the EM
#' algorithm used for learning). This hyper-posterior distribution, evaluated
#' on a grid of inputs provided through the \code{grid_inputs} argument, is a
#' key component for making prediction in Magma, and is required in the function
#' \code{\link{pred_magma}}.
#'
#' @param trained_model A list, containing  the information coming from a
#'    Magma model, previously trained using the \code{\link{train_magma}}
#'    function. If \code{trained_model} is not provided, the arguments
#'    \code{data}, \code{hp_0}, \code{hp_t}, \code{kern_0}, and \code{kern_t}
#'    are all required.
#' @param data A tibble or data frame. Required columns: 'Input',
#'    'Output'. Additional columns for covariates can be specified.
#'    The 'Input' column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    'Output' column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference 'Input'. Recovered from \code{trained_model} if not
#'    provided.
#' @param hp_0 A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_0}. Recovered from \code{trained_model} if not
#'    provided.
#' @param hp_t A tibble or data frame of hyper-parameters
#'    associated with \code{kern_t}. Recovered from \code{trained_model} if not
#'    provided.
#' @param kern_0 A kernel function, associated with the mean GP.
#'    Several popular kernels
#'    (see \href{https://www.cs.toronto.edu/~duvenaud/cookbook/}{The Kernel
#'    Cookbook}) are already implemented and can be selected within the
#'    following list:
#'    - "SE": (default value) the Squared Exponential Kernel (also called
#'        Radial Basis Function or Gaussian kernel),
#'    - "LIN": the Linear kernel,
#'    - "PERIO": the Periodic kernel,
#'    - "RQ": the Rational Quadratic kernel.
#'    - convolution_kernel : the Convolution kernel used to manage MO scenario.
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not). Recovered from
#'    \code{trained_model} if not provided.
#' @param kern_t A kernel function, associated with the task GPs. ("SE",
#'    "PERIO", "RQ" and convolution_kernel are aso available here). Recovered from
#'    \code{trained_model} if not provided.
#' @param prior_mean Hyper-prior mean parameter of the mean GP. This argument,
#'    can be specified under various formats, such as:
#'    - NULL (default). The hyper-prior mean would be set to 0 everywhere.
#'    - A number. The hyper-prior mean would be a constant function.
#'    - A vector of the same length as all the distinct Input values in the
#'     \code{data} argument. This vector would be considered as the evaluation
#'     of the hyper-prior mean function at the training Inputs.
#'    - A function. This function is defined as the hyper-prior mean.
#'    - A tibble or data frame. Required columns: Input, Output. The Input
#'     values should include at least the same values as in the \code{data}
#'     argument.
#' @param grid_inputs A vector or a data frame, indicating the grid of
#'    additional reference inputs on which the mean process' hyper-posterior
#'    should be evaluated. If grid_inputs is a data frame, it should contain 3
#'    mandatory columns : 'Input_ID', 'Input' and 'Output_ID'.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#'
#' @return A list gathering the parameters of the mean processes'
#'         hyper-posterior distributions, namely:
#'         \itemize{
#'           \item mean: A tibble, the hyper-posterior mean parameter
#'                 evaluated at each training \code{Input}.
#'           \item cov: A matrix, the covariance parameter for the
#'                 hyper-posterior distribution of the mean process.
#'           \item pred: A tibble, the predicted mean and variance at
#'                 \code{Input} for the mean process' hyper-posterior
#'                 distribution under a format that allows the direct
#'                 visualisation as a GP prediction.
#'          }
#'
#' @export
#'
#' @examples
#' TRUE
hyperposterior <- function(trained_model = NULL,
                           data = NULL,
                           hp_0 = NULL,
                           hp_t = NULL,
                           kern_0 = NULL,
                           kern_t = NULL,
                           prior_mean = NULL,
                           grid_inputs = NULL,
                           pen_diag = 1e-10) {

  ## Check whether a model trained by train_magma() is provided
  if(trained_model %>% is.null()){
    ## Check whether all mandatory arguments are present otherwise
    if(is.null(data)|is.null(hp_0)|is.null(hp_t)|
       is.null(kern_0)|is.null(kern_t)){
      stop(
        "If no 'trained_model' argument is provided, the arguments 'data', ",
        "'hp_0', 'hp_t' 'kern_0', and 'kern_t' are all required."
      )
    }
  } else {
    ## For each argument, retrieve the value from 'trained_model' if missing
    if(data %>% is.null()){data = trained_model$ini_args$data}
    if(hp_0 %>% is.null()){hp_0 = trained_model$hp_0}
    if(hp_t %>% is.null()){hp_t = trained_model$hp_t}
    if(kern_0 %>% is.null()){kern_0 = trained_model$ini_args$kern_0}
    if(kern_t %>% is.null()){kern_t = trained_model$ini_args$kern_t}
    if(prior_mean %>% is.null()){prior_mean = trained_model$ini_args$prior_mean}

  }

  data <- data %>%
    group_by(Task_ID, Output_ID, Output, Input_ID) %>%
    # Add a unique number observation for the group
    mutate(obs_num = row_number()) %>%
    ungroup()

  ## To create the 'Reference' column as in the old MagmaClustR tibble format, we
  # need to pivot data to obtain one row per observation of (Task_ID, Output_ID).
  # In other words, inputs are no longer in "short" format; instead, we have one
  # column per input.
  data <- data %>%
    tidyr::pivot_wider(
      names_from = Input_ID,
      values_from = Input,
      names_prefix = "Input_"
    ) %>%
    # Keep 6 significant digits for Inputs to avoid numerical issues
    dplyr::mutate(across(starts_with("Input_"), ~ round(.x, 6))) %>%
    rowwise() %>%
    dplyr::mutate(
      Reference = paste(
        # Create output's prefix
        paste0("o", Output_ID),
        # Create the reference for each Output_ID
        paste(c_across(starts_with("Input_")), collapse = ":"),
        # Join output's prefix and reference
        sep = ";"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-obs_num)

  ## Check that tasks do not have duplicate inputs for each output
  task_duplicates <- data %>%
    dplyr::count(Task_ID, Reference) %>%
    dplyr::filter(n > 1)

  if (nrow(task_duplicates) > 0) {
    stop("Error: At least one task has duplicates, i.e. several 'Output' values",
         " for the same 'Output_ID'.")
  }

  ## Get input column names
  if (!("Reference" %in% (names(data)))) {
    names_col <- data %>%
      dplyr::select(- c(Task_ID, Output_ID, Output)) %>%
      names()
  } else {
    names_col <- data %>%
      dplyr::select(- c(Task_ID,  Output_ID, Output, Reference)) %>%
      names()
  }

  if (grid_inputs %>% is.null()) {
    ## Extract the union of all reference inputs provided in the training data
    all_inputs <- data %>%
      dplyr::select(-c(Task_ID, Output_ID, Output)) %>%
      unique() %>%
      tidyr::separate(Reference,
                      into = c("Output_ID_temp", "Input_temp"),
                      sep = ";",
                      remove = FALSE) %>%
      dplyr::mutate(Input_temp_numeric = as.numeric(Input_temp)) %>%
      dplyr::arrange(Output_ID_temp, Input_temp_numeric) %>%
      dplyr::select(c(Input_1, Reference))

    all_input <- all_inputs %>%
      ## DO NOT arrange Reference because of the lexicographic order
      # (not armful but unnecessary in Magma case, tragic in MO case)
      dplyr::pull(Reference)
    cat(
      "The argument 'grid_inputs' is NULL, the hyper-posterior distribution",
      "will only be evaluated on observed Input from 'data'.\n \n"
    )

  } else {
    ## If 'grid_input' is a vector, convert it to the correct format. We suppose
    ## that the vector form is only used when the user works with single input
    ## single output
    if(grid_inputs %>% is.vector()){
      grid_inputs <- tibble::tibble('Input' = grid_inputs,
                                    'Input_ID' = rep("1", length(grid_inputs)),
                                    'Output_ID' = rep(as.factor("1"), length(grid_inputs)))
    }

    ## Define the union among all reference Inputs and a specified grid
    grid_inputs <- grid_inputs %>%
      dplyr::group_by(Input_ID) %>%
      dplyr::mutate(id_ligne = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(
        names_from = Input_ID,
        values_from = Input,
        names_prefix = "Input_"
      ) %>%
      dplyr::select(-id_ligne) %>%
      # Keep 6 significant digits for Inputs to avoid numerical issues
      dplyr::mutate(across(starts_with("Input_"), ~ round(.x, 6))) %>%
      rowwise() %>%
      dplyr::mutate(
        Reference = paste(
          # Create output's prefix
          paste0("o", Output_ID),
          # Create the reference for each Output_ID
          paste(c_across(starts_with("Input_")), collapse = ":"),
          # Join output's prefix and reference
          sep = ";"
        )
      )

    all_inputs <- data %>%
      dplyr::select(Reference, tidyselect::all_of(names_col), Output_ID) %>%
      dplyr::union(grid_inputs) %>%
      unique() %>%
      tidyr::separate(Reference,
                      into = c("Output_ID_temp", "Input_temp"),
                      sep = ";",
                      remove = FALSE) %>%
      dplyr::mutate(Input_temp_numeric = as.numeric(Input_temp)) %>%
      dplyr::arrange(Output_ID_temp, Input_temp_numeric) %>%
      dplyr::select(c(Input_1, Reference, Output_ID))
    # DO NOT arrange Reference because of the lexicographic order
    # (not armful but unnecessary in Magma case, tragic in MO case)

    all_input <- all_inputs %>% dplyr::pull(Reference)

    if(length(data$Output_ID %>% unique()) == 1){
      all_inputs <- all_inputs %>%
        dplyr::select(-Output_ID)
    }
  }

  ## Initialise m_0 according to the value provided by the user
  if (prior_mean %>% is.null()) {
    m_0 <- rep(0, length(all_input))
    cat(
      "The 'prior_mean' argument has not been specified. The hyper_prior mean",
      "function is thus set to be 0 everywhere for all outputs.\n \n"
    )
  } else if (prior_mean %>% is.vector()) {

    if (length(prior_mean) == length(all_input)) {
      m_0 <- prior_mean
    } else if (length(prior_mean) == length(data$Output_ID %>% unique())) {
      # Get the unique and sorted Output_IDs
      unique_outputs_sorted <- data$Output_ID %>% unique() %>% sort()

      # Create a lookup table: "o1" -> prior_mean[1], "o2" -> prior_mean[2], etc.
      # This assumes the prior_mean vector is provided in the sorted order of
      # Output_IDs.
      prior_mean_map <- setNames(prior_mean, paste0("o", unique_outputs_sorted))

      # Extract the prefix ("o1", "o2", etc.) from each element in all_input
      all_input_prefixes <- stringr::str_extract(all_input, "o[0-9]+")

      # Build m_0 using the lookup table; it will automatically repeat the correct
      # value for each prefix.
      m_0 <- prior_mean_map[all_input_prefixes] %>% unname()

      cat(
        "A constant hyper_prior mean has been set for each output.\n \n"
      )
    } else {
      stop(
        "The 'prior_mean' argument is of length ", length(prior_mean),
        ", whereas the grid of training inputs is of length ",
        length(all_input)
      )
    }

  } else if (prior_mean %>% is.function()) {
    m_0 <- prior_mean(all_inputs %>%
                        dplyr::select(-.data$Reference)
    )
  } else {
    stop(
      "Incorrect format for the 'prior_mean' argument. Please read ",
      "?hyperposterior() for details."
    )
  }


  ## Certify that IDs are of type 'character'
  data$Task_ID <- data$Task_ID %>% as.character()

  if('Output_ID' %in% names(all_inputs)){
    all_inputs$Output_ID <- all_inputs$Output_ID %>% as.factor()
  }

  if(length(data$Output_ID %>% unique()) > 1){
    # Compute the convolutional covariance matrix of the mean process
    cov_0 <- kern_to_cov(input = all_inputs,
                         kern = kern_0,
                         hp = hp_0)

    all_references <- rownames(cov_0)
    inv_0 <- cov_0 %>% chol_inv_jitter(pen_diag = pen_diag)
    # matrixcalc::is.positive.semi.definite(inv_0)
    # Re-apply the stored names to the inverted matrix
    dimnames(inv_0) <- list(all_references, all_references)
  } else {
    cov_0 <- kern_to_cov(input = all_inputs,
                         kern = kern_0,
                         hp = hp_0 %>%
                           dplyr::select(-Output_ID))

    all_references <- rownames(cov_0)
    inv_0 <- cov_0 %>% chol_inv_jitter(pen_diag = pen_diag)
    # Re-apply the stored names to the inverted matrix
    dimnames(inv_0) <- list(all_references, all_references)
  }


  list_inv_t <- list()
  list_ID_task <- unique(data$Task_ID)

  list_output_ID <-  data$Output_ID %>% unique()

  # For each task, compute its full multi-output inverse covariance matrix
  for (t in list_ID_task) {
    # Isolate the data and HPs for the current task
    db_t <- data %>% dplyr::filter(Task_ID == t) %>%
      dplyr::select(-c(Output, Task_ID))
    hp_t_indiv <- hp_t %>% dplyr::filter(Task_ID == t)

    if(length(list_output_ID) > 1){
      # Call kern_to_cov directly.
      # It will handle the multi-output structure and the noise addition internally.
      # 'kern_t' is expected to be the 'convolution_kernel' function.
      K_task_t <- kern_to_cov(
        input = db_t,
        kern = kern_t,
        hp = hp_t_indiv
      )
    } else{
      # Extract all_inputs to call kern_to_cov() on the single output case
      all_inputs_t <- data %>%
        dplyr::filter(Task_ID == t) %>%
        dplyr::select(-c(Task_ID, Output, Output_ID)) %>%
        unique()
        ## DO NOT arrange Reference because of the lexicographic order
        # (not armful but unnecessary in Magma case, tragic in MO case)

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

    # Re-apply the stored names to the inverted matrix
    dimnames(K_inv_t) <- list(task_references, task_references)

    # Add the inverted matrix to the list
    # The rownames are already correctly set by kern_to_cov
    list_inv_t[[t]] <- K_inv_t
  }

  ## Update the posterior distribution
  # Create a named list of output values, split by task
  list_output_t <- base::split(data$Output, list(data$Task_ID))

  ## Update Posterior Inverse Covariance
  post_inv <- inv_0
  for (inv_t in list_inv_t) {
    # Find the common input points between the mean process and the current task
    co_input <- intersect(row.names(inv_t), row.names(post_inv))

    # Add the task's contribution to the posterior inverse covariance
    post_inv[co_input, co_input] <- post_inv[co_input, co_input] +
      inv_t[co_input, co_input]
  }

  if(length(data$Output_ID %>% unique()) > 1){
    post_cov <- post_inv %>%
      chol_inv_jitter(pen_diag = pen_diag) %>%
      `rownames<-`(all_references) %>%
      `colnames<-`(all_references)
  } else {
    post_cov <- post_inv %>%
      chol_inv_jitter(pen_diag = pen_diag) %>%
      `rownames<-`(all_inputs$Reference)%>%
      `colnames<-`(all_inputs$Reference)
  }

  ## Update Posterior Mean
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

  ## Format the mean parameter of the hyper-posterior distribution
  mean <- tibble::tibble(all_inputs,
                         "Output" = post_mean
  )

  ## Format the GP prediction of the hyper-posterior mean (for direct plot)
  pred <- tibble::tibble(all_inputs,
                         "Mean" = post_mean,
                         "Var" = post_cov %>% diag() %>% as.vector()
  ) %>%
    dplyr::select(-Reference)

  list(
    "mean" = mean,
    "cov" = post_cov,
    "pred" = pred
  ) %>%
    return()
}


#' Magma prediction
#'
#' Compute the posterior predictive distribution in Magma. Providing data of any
#' new task, its trained hyper-parameters and a previously trained
#' Magma model, the predictive distribution is evaluated on any arbitrary inputs
#' that are specified through the 'grid_inputs' argument.
#'
#' @param data A tibble or data frame. Required columns: \code{Task_ID},
#'    \code{Input_ID}, \code{Input}, \code{Output_ID}, \code{Output}. Additional
#'    columns for covariates can be specified.
#'    The \code{Task_ID} column contains the unique names/codes used to identify each
#'    task (or batch of data).
#'    The \code{Input_ID} contains the unique names/codes used to identify each
#'    explanatory variable.
#'    The \code{Input} column should define the value of the explanatory
#'    variables that are used as reference for the observations (e.g. time for
#'    longitudinal data).
#'    The \code{Output_ID} contains the unique names/codes used to identify each
#'    output (response) variable.
#'    The \code{Output} column specifies the observed values (the response
#'    variables).
#' @param trained_model A list, containing  the information coming from a
#'    Magma model, previously trained using the \code{\link{train_magma}}
#'    function.
#' @param grid_inputs The grid of inputs (reference Inputs) values
#'    on which the GP should be evaluated. Ideally, this argument should be a
#'    tibble or a data frame, providing the same columns as \code{data}, except
#'    'Output' and 'Task_ID'. Nonetheless, in cases where \code{data} provides only one
#'    'Input' column, the \code{grid_inputs} argument can be NULL (default) or a
#'    vector. This vector would be used as reference input for prediction and if
#'    NULL, a vector of length 500 is defined, ranging between the min and max
#'    Input values of \code{data}. If \code{grid_inputs} is a vector, we replicate
#'    the vector as many times as the number of outputs.
#' @param hp A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern}. The columns/elements should be named
#'    according to the hyper-parameters that are used in \code{kern}. The
#'    function \code{\link{train_gp}} can be used to learn maximum-likelihood
#'    estimators of the hyper-parameters.
#' @param kern A kernel function, defining the covariance structure of the GP.
#'    Several popular kernels
#'    (see \href{https://www.cs.toronto.edu/~duvenaud/cookbook/}{The Kernel
#'    Cookbook}) are already implemented and can be selected within the
#'    following list:
#'    - "SE": (default value) the Squared Exponential Kernel (also called
#'        Radial Basis Function or Gaussian kernel),
#'    - "LIN": the Linear kernel,
#'    - "PERIO": the Periodic kernel,
#'    - "RQ": the Rational Quadratic kernel.
#'    - convolution_kernel : the Convolution kernel used to manage MO scenario.
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#' @param hyperpost A list, containing the elements 'mean' and 'cov', the
#'    parameters of the hyper-posterior distribution of the mean process.
#'    Typically, this argument should come from a previous learning using
#'    \code{\link{train_magma}}, or a previous prediction with
#'    \code{\link{pred_magma}}, with the argument \code{get_hyperpost} set to
#'    TRUE. The 'mean' element should be a data frame with two columns 'Input'
#'    and 'Output'. The 'cov' element should be a covariance matrix with
#'    colnames and rownames corresponding to the 'Input' in 'mean'. In all
#'    cases, the column 'Input' should contain all the values appearing both in
#'    the 'Input' column of \code{data} and in \code{grid_inputs}.
#' @param get_hyperpost A logical value, indicating whether the hyper-posterior
#'    distribution of the mean process should be returned. This can be useful
#'    when planning to perform several predictions on the same grid of inputs,
#'    since recomputation of the hyper-posterior can be prohibitive for high
#'    dimensional grids.
#' @param get_full_cov A logical value, indicating whether the full posterior
#'    covariance matrix should be returned.
#' @param plot A logical value, indicating whether a plot of the results is
#'    automatically displayed.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#'
#' @return A list, containing :
#'    - a tibble, representing Magma predictions as two column 'Mean' and
#'   'Var', evaluated on the \code{grid_inputs}. The column 'Input' and
#'   additional covariates columns are associated to each predicted values.
#'    If the \code{get_full_cov} or \code{get_hyperpost} arguments are TRUE,
#'    the function returns a list, in which the tibble described above is
#'    defined as 'pred_gp' and the full posterior covariance matrix is
#'    defined as 'cov', and the hyper-posterior distribution of the mean process
#'    is defined as 'hyperpost';
#'    - a ponderation matrix, corresponding to t(cov_crossed) %*% inv_obs.
#' @export
#'
#' @examples
#' TRUE
pred_magma <- function(data = NULL,
                       trained_model = NULL,
                       grid_inputs = NULL,
                       hp = NULL,
                       kern = "SE",
                       hyperpost = NULL,
                       get_hyperpost = FALSE,
                       get_full_cov = FALSE,
                       plot = TRUE,
                       pen_diag = 1e-10) {

  ## Return the mean process if no data is provided
  if(data %>% is.null()){
    ## Check whether trained_model is provided
    if(trained_model %>% is.null()){
      stop(
        "If 'data' is not provided, the 'trained_model' argument is needed ",
        "to provide the mean process as a generic prediction."
      )
    }
    else{
      ## Return the mean process as a generic prediction
      if(grid_inputs %>% is.null()){
        hyperpost = trained_model$hyperpost

      } else{
        ## Recompute the mean process at the required inputs if necessary
        hyperpost = hyperposterior(
          trained_model = trained_model,
          grid_inputs = grid_inputs,
          kern_0 = trained_model$ini_args$kern_0,
          kern_t = trained_model$ini_args$kern_t,
          hp_0 = trained_model$hp_0,
          hp_t = trained_model$hp_t,
        )

        ## Retain only grid_inputs for display purposes
        grid_inputs_wide <- grid_inputs %>%
          dplyr::group_by(Input_ID) %>%
          dplyr::mutate(id_ligne = dplyr::row_number()) %>%
          dplyr::ungroup() %>%
          tidyr::pivot_wider(
            names_from = Input_ID,
            values_from = Input,
            names_prefix = "Input_"
          )

        shared_columns <- intersect(
          names(grid_inputs_wide),
          names(hyperpost$pred)
        )

        hyperpost$pred <- grid_inputs_wide %>%
          dplyr::inner_join(hyperpost$pred, by = shared_columns) %>%
          dplyr::select(-id_ligne)
      }

      pred = hyperpost$pred

      ## Check whether the full posterior covariance should be returned
      if (get_full_cov) {
        pred <- list("pred" = pred, "cov" = hyperpost$cov)
      }
      return(pred)
    }
  }

  ## Keep a version of raw 'data-to-predict' (used only if plot == TRUE)
  raw_data <- data

  ## Remove the 'ID' column if present
  if ("Task_ID" %in% names(data)) {
    if (dplyr::n_distinct(data$Task_ID) > 1) {
      stop(
        "Problem in the 'Task_ID' column: different values are not allowed. ",
        "The prediction can only be performed for one task."
      )
    }
    data <- data %>% dplyr::select(-Task_ID)
  }

  data <- data %>%
    dplyr::group_by(Output_ID, Output, Input_ID) %>%
    # Add a unique number observation for the group
    dplyr::mutate(obs_num = row_number()) %>%
    dplyr::ungroup()

  ## To create the 'Reference' column as in the old MagmaClustR tibble format, we
  # need to pivot data to obtain one row per observation of (Task_ID, Output_ID).
  # In other words, inputs are no longer in "short" format; instead, we have one
  # column per input.
  data <- data %>%
    tidyr::pivot_wider(
      names_from = Input_ID,
      values_from = Input,
      names_prefix = "Input_"
    ) %>%
    # Keep 6 significant digits for Inputs to avoid numerical issues
    dplyr::mutate(across(starts_with("Input_"), ~ round(.x, 6))) %>%
    rowwise() %>%
    dplyr::mutate(
      Reference = paste(
        # Create output's prefix
        paste0("o", Output_ID),
        # Create the reference for each Output_ID
        paste(c_across(starts_with("Input_")), collapse = ":"),
        # Join output's prefix and reference
        sep = ";"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-obs_num)

  ## Get input column names
  if ("Reference" %in% names(data)) {
    names_col <- data %>%
      dplyr::select(- c(Output_ID, Output, Reference)) %>%
      names()
  } else {
    names_col <- data %>%
      dplyr::select(-c(Output_ID, Output)) %>%
      names()
  }

  ## Extract the observed Output (data points)
  data_obs <- data %>%
    dplyr::pull(Output)

  ## Name elements of data_obs with 'Reference' column
  names(data_obs) <- data$Reference

  ## Extract the observed (reference) Input
  input_obs <- data %>%
    dplyr::pull(Reference)

  ## Extract the observed inputs (reference Input + covariates)
  inputs_obs <- data %>%
    dplyr::select(-c(Output))

  ## Define the target inputs to predict
  if (grid_inputs %>% is.null()) {
    set_grid <- function(data, size_grid) {
      seq(data %>% min(),
          data %>% max(),
          length.out = size_grid
      ) %>%
        return()
    }

    if (inputs_obs %>% names() %>% length() == 3) {
      size_grid <- 500
    } else if (inputs_obs %>% names() %>% length() > 3) {
      size_grid <- 1000^(1 / (ncol(inputs_obs) - 1)) %>% round()
      ## floor instead of round ?
    } else {
      stop(
        "The 'grid_inputs' argument should be a either a numerical vector ",
        "or a data frame depending on the context. Please read ?pred_magma()."
      )
    }

    # Create a unique grid of inputs (which will be replicated for each Output_ID)
    base_grid <- purrr::map(
      data %>% dplyr::select(tidyselect::all_of(names_col)),
      set_grid,
      size_grid
    ) %>%
      purrr::set_names(paste0("Input_", seq_along(.))) %>%
      expand.grid() %>%
      tibble::as_tibble()

    unique_outputs <- data %>%
      dplyr::distinct(Output_ID)

    # Cross base_grid with unique_outputs to replicate base_grid as many times as
    # unique_outputs and create 'Reference' column
    inputs_pred <- tidyr::crossing(unique_outputs, base_grid) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        Reference = paste(
          paste0("o", Output_ID),
          paste(c_across(starts_with("Input_")), collapse = ":"),
          sep = ";"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(Reference)

    input_pred <- inputs_pred %>% dplyr::pull(Reference)

  } else if (grid_inputs %>% is.vector()) {
    ## Test whether 'data' only provide the Input column and no covariates
    if (inputs_obs %>% names() %>% length() == 3) {
      input_temp <- grid_inputs %>%
        signif() %>%
        sort() %>%
        unique()

      inputs_pred <- tibble::tibble("Input_1" = rep(input_temp,
                                                    times = nrow(unique_outputs)),
                                    "Output_ID" = rep(unique_outputs %>%
                                                        purrr::as_vector(),
                                                      each = length(input_temp))) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          Reference = paste(
            paste0("o", Output_ID),
            paste(c_across(starts_with("Input_")), collapse = ":"),
            sep = ";"
          )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(Reference)

      input_pred <- inputs_pred %>% dplyr::pull(Reference)

    } else {
      stop(
        "The 'grid_inputs' argument should be a either a numerical vector ",
        "or a data frame depending on the context. Please read ?pred_magma()."
      )
    }
  } else if (grid_inputs %>% is.data.frame()) {
    grid_inputs <- grid_inputs %>%
      dplyr::group_by(Input_ID) %>%
      dplyr::mutate(id_ligne = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(
        names_from = Input_ID,
        values_from = Input,
        names_prefix = "Input_"
      ) %>%
      dplyr::select(-id_ligne) %>%
      # Keep 6 significant digits for Inputs to avoid numerical issues
      dplyr::mutate(across(starts_with("Input_"), ~ round(.x, 6))) %>%
      rowwise() %>%
      dplyr::mutate(
        Reference = paste(
          # Create output's prefix
          paste0("o", Output_ID),
          # Create the reference for each Output_ID
          paste(c_across(starts_with("Input_")), collapse = ":"),
          # Join output's prefix and reference
          sep = ";"
        )
      )

    ## Test whether 'data' has the same columns as grid_inputs
    if (names(inputs_obs) %>% setequal(names(grid_inputs))) {
      input_pred <- grid_inputs %>%
        # DO NOT arrange Reference because of the lexicographic order
        # (not armful but unnecessary in Magma case, tragic in MO case)
        dplyr::pull(Reference)

      inputs_pred <- grid_inputs %>%
        # DO NOT arrange Reference because of the lexicographic order
        # (not armful but unnecessary in Magma case, tragic in MO case)
        dplyr::select(names(inputs_obs))
    } else {
      stop(
        "The 'grid_inputs' argument should provide a column 'Input_ID', 'Input' ",
        " and 'Output_ID'."
      )
    }
  } else {
    stop(
      "The 'grid_inputs' argument should be a either a numerical vector ",
      "or a data frame depending on the context. Please read ?pred_gp()."
    )
  }


  ## Define the union of all distinct reference Input
  all_inputs <- dplyr::union(inputs_obs, inputs_pred)
  # DO NOT arrange Reference because of the lexicographic order
  # (not armful but unnecessary in Magma case, tragic in MO case)
  all_input <- all_inputs$Reference


  ## Check whether the hyper-posterior is provided and recompute if necessary
  if (hyperpost %>% is.null()) {
    if (trained_model %>% is.null()) {
      stop(
        "If the 'hyperpost' argument is NULL, the 'trained_model' ",
        "should be provided, in order to extract or recompute mean process' ",
        "hyper-posterior distribution evaluated on the correct inputs."
      )
    } else {
      ## Get the hyper-posterior distribution from the trained model
      hyperpost <- trained_model$hyperpost
    }
  }

  ## Check hyperpost format
  if ((hyperpost %>% is.list()) &
      (!is.null(hyperpost$mean)) &
      (!is.null(hyperpost$cov))
  ) {
    ## Check hyperpost format (in particular presence of all reference Input)
    if (!all(all_input %in% hyperpost$mean$Reference)) {
      if (trained_model %>% is.null()) {
        stop(
          "hyperpost is not evaluated at the correct inputs, please use the ",
          "'trained_model' argument instead."
        )
      }
      cat(
        "The hyper-posterior distribution of the mean process provided in",
        "'hyperpost' argument isn't evaluated on the expected inputs.\n \n",
        "Start evaluating the hyper-posterior on the correct inputs...\n \n"
      )

      # Reformat all_inputs into a hyperposterior() friendly format
      all_inputs_long <- all_inputs %>%
        dplyr::select(-Reference) %>%
        tidyr::pivot_longer(
          cols = starts_with("Input_"),
          names_to = "Input_ID",
          values_to = "Input"
        ) %>%
        dplyr::mutate(Input_ID = stringr::str_remove(Input_ID, "Input_"))

      hyperpost <- hyperposterior(
        data = trained_model$ini_args$data,
        kern_0 = trained_model$ini_args$kern_0,
        kern_t = trained_model$ini_args$kern_t,
        hp_0 = trained_model$hp_0,
        hp_t = trained_model$hp_t,
        prior_mean = trained_model$ini_args$prior_mean,
        grid_inputs = all_inputs_long,
        pen_diag = pen_diag
      )
      cat("Done!\n \n")
    }
  } else {
    stop(
      "The format of the 'hyperpost' argument is not as expected. Please ",
      "read ?pred_magma() for details."
    )
  }

  ## Extract the mean parameter from the hyper-posterior

  mean_obs <- hyperpost$mean %>%
    dplyr::filter(Reference %in% input_obs) %>%
    # DO NOT arrange Reference because of the lexicographic order
    # (not armful but unnecessary in Magma case, tragic in MO case)
    dplyr::pull(Output)

  names(mean_obs) <- (hyperpost$mean %>%
                        dplyr::filter(Reference %in% input_obs))$Reference

  mean_pred <- hyperpost$mean %>%
    dplyr::filter(Reference %in% input_pred) %>%
    # DO NOT arrange Reference because of the lexicographic order
    # (not armful but unnecessary in Magma case, tragic in MO case)
    dplyr::pull(Output)

  names(mean_pred) <- (hyperpost$mean %>%
                         dplyr::filter(Reference %in% input_pred))$Reference

  ## Extract the covariance sub-matrices from the hyper-posterior
  post_cov_obs <- hyperpost$cov[
    as.character(input_obs),
    as.character(input_obs)
  ]
  post_cov_pred <- hyperpost$cov[
    as.character(input_pred),
    as.character(input_pred)
  ]
  post_cov_crossed <- hyperpost$cov[
    as.character(input_obs),
    as.character(input_pred)
  ]

  ## Extract or learn the hyper-parameters if not provided
  if (hp %>% is.null()) {
    if (!is.null(trained_model)) {
      ## Check whether hyper-parameters are shared if we have 'trained_model'
      if (
        tryCatch(trained_model$ini_args$shared_hp_tasks, error = function(e) FALSE)
      ) {
        ## Extract the hyper-parameters common to all 't'
        hp <- trained_model$hp_t %>%
          dplyr::slice(1:length(all_inputs$Output_ID %>% unique())) %>%
          dplyr::select(-Task_ID)
      } else if (kern %>% is.function()) {
        stop(
          "When using a custom kernel function the 'hp' argument is ",
          "mandatory, in order to provide the name of the hyper-parameters. ",
          "You can use the function 'hp()' to easily generate a tibble of ",
          "random hyper-parameters with the desired format, or use ",
          "'train_gp()' to learn ML estimators for a better fit."
        )
      } else if (kern %>% is.character()) {
        hp <- quiet(
          train_gp(data,
                   prior_mean = NULL,
                   ini_hp = hp(kern, noise = T),
                   kern = kern,
                   hyperpost = hyperpost,
                   pen_diag = pen_diag
          )
        )
        cat(
          "The 'hp' argument has not been specified. The 'train_gp()' function",
          "(with random initialisation) has been used to learn ML estimators",
          "for the hyper-parameters associated with the 'kern' argument.\n \n"
        )
      }
    } else if (kern %>% is.function()) {
      stop(
        "When using a custom kernel function the 'hp' argument is ",
        "mandatory, in order to provide the name of the hyper-parameters. ",
        "You can use the function 'hp()' to easily generate a tibble of random",
        " hyper-parameters with the desired format, or use 'train_gp()' to ",
        "learn ML estimators for a better fit of data."
      )
    } else if (kern %>% is.character()) {
      hp <- quiet(
        train_gp(data,
                 prior_mean = NULL,
                 ini_hp = hp(kern, noise = T),
                 kern = kern,
                 hyperpost = hyperpost,
                 pen_diag = pen_diag
        )
      )
      cat(
        "The 'hp' argument has not been specified. The 'train_gp()' function",
        "(with random initialisation) has been used to learn ML estimators",
        "for the hyper-parameters associated with the 'kern' argument.\n \n"
      )
    } else {
      stop(
        "Incorrect format for the 'kern' argument. Please read ?pred_gp() for ",
        "details."
      )
    }
  }

  if(length(inputs_obs$Output_ID %>% unique()) == 1){
    inputs_obs <- inputs_obs %>%
      dplyr::select(-Output_ID)
    # DO NOT arrange Reference because of the lexicographic order
    # (not armful but unnecessary in Magma case, tragic in MO case)

    inputs_pred <- inputs_pred %>%
      dplyr::select(-Output_ID)
    # DO NOT arrange Reference because of the lexicographic order
    # (not armful but unnecessary in Magma case, tragic in MO case)
  }

  ## Sum the covariance matrices on observed inputs and compute the inverse
  cov_obs <- kern_to_cov(inputs_obs, kern, hp) + post_cov_obs

  inv_obs <- cov_obs %>%
    chol_inv_jitter(pen_diag = pen_diag) %>%
    `rownames<-`(as.character(input_obs)) %>%
    `colnames<-`(as.character(input_obs))

  ## Remove the noise of the hp for evaluating some of the sub-matrix
  if ("noise" %in% names(hp)) {
    hp_rm_noi <- hp %>% dplyr::select(-noise)
    noise <- exp(hp[["noise"]])
  } else {
    hp_rm_noi <- hp
    noise <- 0
  }

  ## Compute the required sub-matrix for prediction
  cov_pred <- kern_to_cov(inputs_pred, kern, hp_rm_noi) + post_cov_pred
  cov_crossed <- kern_to_cov(inputs_obs, kern, hp_rm_noi,
                             input_2 = inputs_pred
  ) +
    post_cov_crossed

  ## Compute the posterior mean of a GP
  pred_mean <- (mean_pred +
                  t(cov_crossed) %*% inv_obs %*% (data_obs - mean_obs)) %>%
    as.vector()

  ## Compute the posterior covariance matrix of a GP
  pred_cov <- cov_pred - t(cov_crossed) %*% inv_obs %*% cov_crossed

  if(length(all_inputs$Output_ID %>% unique()) > 1){
    ## Get the pred_cov rownames and colnames to create a noise vector with
    ## correct dimensions
    point_names <- rownames(pred_cov)
    output_ids <- as.numeric(sub("^o(\\d+);.*", "\\1", point_names))
    full_noise_vector <- noise[output_ids]
  } else {
    full_noise_vector <- noise
  }


  ## Create a tibble of values and associated uncertainty from a GP prediction
  pred_gp <- tibble::tibble(
    "Mean" = pred_mean,
    "Var" = diag(pred_cov) + full_noise_vector
  ) %>%
    dplyr::mutate(inputs_pred) %>%
    dplyr::select(-Reference)

  # Compute ponderation matrix to visualize the contribution of each observation
  # to each predicted point
  ponderation_matrix <- t(cov_crossed) %*% inv_obs

  ## Display the graph of the prediction if expected
  if (plot) {
    ## Check whether training data are available
    if (trained_model %>% is.null()) {
      data_train <- NULL
    } else {
      data_train <- trained_model$ini_args$data
    }

    ## Add 'cov' to display samples
    res <- list("pred" = pred_gp)
    res[["cov"]] <- pred_cov

    ## Display samples only in 1D and Credible Interval otherwise
    display_samples = dplyr::if_else(ncol(pred_gp) == 4, TRUE, FALSE)

    if(length(data_train$Output_ID %>% unique()) == 1){
      res$pred$Output_ID <- as.factor("1")
    }

    ## Plot results
    plot_gp(res,
            data = raw_data,
            data_train = data_train %>%
              group_by(Task_ID, Output_ID, Output, Input_ID) %>%
              dplyr::mutate(obs_num = row_number()) %>%
              dplyr::ungroup(),
            prior_mean = hyperpost$mean %>%
              dplyr::select(-Reference),
            samples = display_samples
            ) %>%
      print()
  }

  res <- pred_gp
  ## Check whether posterior covariance or hyper-posterior should be returned
  if (get_full_cov | get_hyperpost) {
    res <- list("pred" = pred_gp)
    if (get_full_cov) {
      res[["cov"]] <- pred_cov
    }
    if (get_hyperpost) {
      res[["hyperpost"]] <- hyperpost
    }
  }

  return(
    list("pred_gp" = res,
         "ponderation_matrix" = ponderation_matrix)
    )
}

#' Magma prediction for ploting GIFs
#'
#' Generate a Magma or classic GP prediction under a format that is compatible
#' with a further GIF visualisation of the results. For a Magma prediction,
#' either the \code{trained_model} or \code{hyperpost} argument is required.
#' Otherwise, a classic GP prediction is applied and the prior mean can be
#' specified through the \code{mean} argument.
#'
#' @param data  A tibble or data frame. Required columns: 'Input',
#'    'Output'. Additional columns for covariates can be specified.
#'    The 'Input' column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    'Output' column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference 'Input'.
#' @param trained_model A list, containing  the information coming from a
#'    Magma model, previously trained using the \code{\link{train_magma}}
#'    function.
#' @param grid_inputs The grid of inputs (reference Input and covariates) values
#'    on which the GP should be evaluated. Ideally, this argument should be a
#'    tibble or a data frame, providing the same columns as \code{data}, except
#'    'Output'. Nonetheless, in cases where \code{data} provides only one
#'    'Input' column, the \code{grid_inputs} argument can be NULL (default) or a
#'    vector. This vector would be used as reference input for prediction and if
#'    NULL, a vector of length 500 is defined, ranging between the min and max
#'    Input values of \code{data}.
#' @param mean Mean parameter of the GP. This argument can be specified under
#'    various formats, such as:
#'    - NULL (default). The mean would be set to 0 everywhere.
#'    - A number. The mean would be a constant function.
#'    - A function. This function is defined as the mean.
#'    - A tibble or data frame. Required columns: Input, Output. The Input
#'     values should include at least the same values as in the \code{data}
#'     argument.
#' @param hp A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern}. The columns/elements should be named
#'    according to the hyper-parameters that are used in \code{kern}. The
#'    function \code{\link{train_gp}} can be used to learn maximum-likelihood
#'    estimators of the hyper-parameters,
#' @param kern A kernel function, defining the covariance structure of the GP.
#'    Several popular kernels
#'    (see \href{https://www.cs.toronto.edu/~duvenaud/cookbook/}{The Kernel
#'    Cookbook}) are already implemented and can be selected within the
#'    following list:
#'    - "SE": (default value) the Squared Exponential Kernel (also called
#'        Radial Basis Function or Gaussian kernel),
#'    - "LIN": the Linear kernel,
#'    - "PERIO": the Periodic kernel,
#'    - "RQ": the Rational Quadratic kernel.
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#' @param hyperpost A list, containing the elements 'mean' and 'cov', the
#'    parameters of the hyper-posterior distribution of the mean process.
#'    Typically, this argument should from a previous learning using
#'    \code{\link{train_magma}}, or a previous prediction with
#'    \code{\link{pred_magma}}, with the argument \code{get_hyperpost} set to
#'    TRUE. The 'mean' element should be a data frame with two columns 'Input'
#'    and 'Output'. The 'cov' element should be a covariance matrix with
#'    colnames and rownames corresponding to the 'Input' in 'mean'. In all
#'    cases, the column 'Input' should contain all the values appearing both in
#'    the 'Input' column of \code{data} and in \code{grid_inputs}.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#'
#' @return A tibble, representing Magma or GP predictions as two column 'Mean'
#'    and 'Var', evaluated on the \code{grid_inputs}. The column 'Input' and
#'    additional covariates columns are associated to each predicted values. An
#'    additional 'Index' column is created for the sake of GIF creation using
#'    the function \code{\link{plot_gif}}
#' @export
#'
#' @examples
#' TRUE
pred_gif <- function(data,
                     trained_model = NULL,
                     grid_inputs = NULL,
                     hyperpost = NULL,
                     mean = NULL,
                     hp = NULL,
                     kern = "SE",
                     pen_diag = 1e-10) {
  ## Remove possible missing data
  data <- data %>% tidyr::drop_na()

  ## Extract the inputs (reference Input + covariates)
  inputs <- data %>% dplyr::select(-.data$Output)
  ## Remove the 'ID' column if present
  if ("ID" %in% names(data)) {
    inputs <- inputs %>% dplyr::select(-.data$ID)
  }
  min_data <- min(data$Input)
  max_data <- max(data$Input)

  ## Define the target inputs to predict
  if (grid_inputs %>% is.null()) {
    ## Test whether 'data' only provide the Input column and no covariates
    if (inputs %>% names() %>% length() == 1) {
      grid_inputs <- tibble::tibble(
        "Input" = seq(min_data, max_data, length.out = 500)
      )
    } else if (inputs %>% names() %>% length() == 2) {
      ## Define a default grid for 'Input'
      grid_inputs <- tibble::tibble(
        "Input" = rep(seq(min_data, max_data, length.out = 20), each = 20)
      )
      ## Add a grid for the covariate
      name_cova <- inputs %>%
        dplyr::select(-.data$Input) %>%
        names()
      cova <- inputs[name_cova]
      grid_inputs[name_cova] <- rep(
        seq(min(cova), max(cova), length.out = 20),
        times = 20
      )
    } else {
      stop(
        "The 'grid_inputs' argument should be a either a numerical vector ",
        "or a data frame depending on the context. Please read ?pred_gp()."
      )
    }
  }

  all_pred <- tibble::tibble()
  for (j in 1:nrow(data)) {
    cat(" =>", j)
    ## Extract the sample of the 'j' first data points
    data_j <- data %>% dplyr::slice(1:j)
    if (!is.null(trained_model) | !is.null(hyperpost)) {
      ## Compute Magma prediction for this sub-dataset
      all_pred <- quiet(pred_magma(data_j,
        trained_model = trained_model,
        hp = hp,
        kern = kern,
        grid_inputs = grid_inputs,
        hyperpost = hyperpost,
        get_hyperpost = FALSE,
        get_full_cov = FALSE,
        plot = FALSE,
        pen_diag = pen_diag
      )) %>%
        dplyr::mutate(Index = j) %>%
        dplyr::bind_rows(all_pred)
    } else {
      all_pred <- quiet(pred_gp(data_j,
        mean = mean,
        hp = hp,
        kern = kern,
        grid_inputs = grid_inputs,
        get_full_cov = FALSE,
        plot = FALSE,
        pen_diag = pen_diag
      )) %>%
        dplyr::mutate(Index = j) %>%
        dplyr::bind_rows(all_pred)
    }
  }
  return(all_pred)
}

#' Compute the hyper-posterior distribution for each cluster in MagmaClust
#'
#' Recompute the E-step of the VEM algorithm in MagmaClust for a new set of
#' reference \code{Input}. Once training is completed, it can be necessary to
#' evaluate the hyper-posterior distributions of the mean processes at specific
#' locations, for which we want to make predictions. This process is directly
#' implemented in the \code{\link{pred_magmaclust}} function but the user
#' might want to use \code{hyperpost_clust} for a tailored control of
#' the prediction procedure.
#'
#' @param trained_model A list, containing  the information coming from a
#'    Magma model, previously trained using the \code{\link{train_magma}}
#'    function. If \code{trained_model} is not provided, the arguments
#'    \code{data}, \code{mixture}, \code{hp_k}, \code{hp_t}, \code{kern_k}, and
#'    \code{kern_t} are all required.
#' @param data A tibble or data frame. Required columns: \code{Task_ID},
#'    \code{Input_ID}, \code{Input}, \code{Output_ID}, \code{Output}. Additional
#'    columns for covariates can be specified.
#'    The \code{Task_ID} column contains the unique names/codes used to identify
#'    each task (or batch of data).
#'    The \code{Input_ID} contains the unique names/codes used to identify each
#'    explanatory variable.
#'    The \code{Input} column should define the value of the explanatory
#'    variables that are used as reference for the observations (e.g. time for
#'    longitudinal data).
#'    The \code{Output_ID} contains the unique names/codes used to identify each
#'    output (response) variable.
#'    The \code{Output} column specifies the observed values (the response
#'    variables).
#' @param mixture A tibble or data frame, indicating the mixture probabilities
#'     of each cluster for each task. Required column: \code{ID}.
#'     Recovered from \code{trained_model} if not
#'     provided.
#' @param hp_k A tibble or data frame of hyper-parameters
#'    associated with \code{kern_k}. Recovered from \code{trained_model} if not
#'    provided.
#' @param hp_t A tibble or data frame of hyper-parameters
#'    associated with \code{kern_t}. Recovered from \code{trained_model} if not
#'    provided.
#' @param kern_k A kernel function, associated with the mean GPs.
#'    Several popular kernels
#'    (see \href{https://www.cs.toronto.edu/~duvenaud/cookbook/}{The Kernel
#'    Cookbook}) are already implemented and can be selected within the
#'    following list:
#'    - "SE": (default value) the Squared Exponential Kernel (also called
#'        Radial Basis Function or Gaussian kernel),
#'    - "LIN": the Linear kernel,
#'    - "PERIO": the Periodic kernel,
#'    - "RQ": the Rational Quadratic kernel.
#'    - convolution_kernel: the Convolution kernel used to manage MO scenario.
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not). Recovered from
#'    \code{trained_model} if not provided.
#' @param kern_t A kernel function, associated with the task GPs. ("SE",
#'    "LIN", PERIO", "RQ" and convolution_kernel are also available here).
#'    Recovered from \code{trained_model} if not provided.
#' @param prior_mean_k The set of hyper-prior mean parameters (m_k) for the K
#'    mean GPs, one value for each cluster.
#'    cluster. This argument can be specified under various formats, such as:
#'    - NULL (default). All hyper-prior means would be set to 0 everywhere.
#'    - A numerical vector of the same length as the number of clusters.
#'    Each number is associated with one cluster, and considered
#'    to be the hyper-prior mean parameter of the cluster (i.e. a constant
#'    function at all \code{Input}).
#'    - A list of functions. Each function is associated with one cluster. These
#'    functions are all evaluated at all \code{Input} values, to provide
#'    specific hyper-prior mean vectors for each cluster.
#' @param grid_inputs A vector or a data frame, indicating the grid of
#'    additional reference inputs on which the mean process' hyper-posterior
#'    should be evaluated.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#'
#' @return A list containing the parameters of the mean processes'
#'         hyper-posterior distribution, namely:
#'         \itemize{
#'            \item mean: A list of tibbles containing, for each cluster, the
#'                  hyper-posterior mean parameters evaluated at each
#'                  \code{Input}.
#'            \item cov: A list of matrices containing, for each cluster, the
#'                  hyper-posterior covariance parameter of the mean process.
#'            \item mixture: A tibble, indicating the mixture probabilities in
#'                  each cluster for each task.
#'          }
#'
#' @export
#'
#' @examples
#' TRUE
#' Compute the hyper-posterior distribution for each cluster in MagmaClust
#'
#' @export
hyperposterior_clust <- function(trained_model = NULL,
                                 data = NULL,
                                 mixture = NULL,
                                 hp_k = NULL,
                                 hp_t = NULL,
                                 kern_k = NULL,
                                 kern_t = NULL,
                                 prior_mean_k = NULL,
                                 grid_inputs = NULL,
                                 pen_diag = 1e-10) {
  ## Check whether a model trained by train_magma() is provided
  if(trained_model %>% is.null()){
    ## Check whether all mandatory arguments are present otherwise
    if(is.null(data)|is.null(prior_mean_k)|is.null(hp_k)|is.null(hp_t)|is.null(mixture)|
       is.null(kern_k)|is.null(kern_t)){
      stop(
        "If no 'trained_model' argument is provided, the arguments 'data', ",
        "'prior_mean_k', mixture', 'hp_k', 'hp_t' 'kern_k', and 'kern_t', ",
        "are all required."
      )
    }
  } else {
    ## For each argument, retrieve the value from 'trained_model' if missing
    if(data %>% is.null()){data = trained_model$ini_args$data}
    if(mixture %>% is.null()){mixture = trained_model$hyperpost$mixture}
    if(hp_k %>% is.null()){hp_k = trained_model$hp_k}
    if(hp_t %>% is.null()){hp_t = trained_model$hp_t}
    if(kern_k %>% is.null()){kern_k = trained_model$ini_args$kern_k}
    if(kern_t %>% is.null()){kern_t = trained_model$ini_args$kern_t}
    if(prior_mean_k %>% is.null()){
      prior_mean_k = trained_model$m_k
    }
  }

  data <- data %>%
    dplyr::group_by(Task_ID, Output_ID, Output, Input_ID) %>%
    # Add a unique number observation for the group
    dplyr::mutate(obs_num = row_number()) %>%
    dplyr::ungroup()

  ## To create the 'Reference' column as in the old MagmaClustR tibble format
  data <- data %>%
    tidyr::pivot_wider(
      names_from = Input_ID,
      values_from = Input,
      names_prefix = "Input_"
    ) %>%
    # Keep 6 significant digits for Inputs to avoid numerical issues
    dplyr::mutate(across(starts_with("Input_"), ~ round(.x, 6))) %>%
    rowwise() %>%
    dplyr::mutate(
      Reference = paste(
        # Create output's prefix
        paste0("o", Output_ID),
        # Create the reference for each Output_ID
        paste(c_across(starts_with("Input_")), collapse = ":"),
        # Join output's prefix and reference
        sep = ";"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-obs_num)

  ## Check that tasks do not have duplicate inputs for each output
  task_duplicates <- data %>%
    dplyr::count(Task_ID, Reference) %>%
    dplyr::filter(n > 1)

  if (nrow(task_duplicates) > 0) {
    stop("Error: At least one task has duplicates, i.e. several 'Output' values",
         " for the same 'Output_ID'.")
  }

  ## Get input column names
  if (!("Reference" %in% (names(data)))) {
    names_col <- data %>%
      dplyr::select(- c(Task_ID, Output_ID, Output)) %>%
      names()
  } else {
    names_col <- data %>%
      dplyr::select(- c(Task_ID,  Output_ID, Output, Reference)) %>%
      names()
  }

  ## Get the number of clusters
  nb_cluster <- length(hp_k$Cluster_ID %>% unique)
  ## Get the name of clusters
  ID_k <- hp_k$Cluster_ID %>%
    unique() %>%
    as.character()

  ## Certify that IDs are of type 'character'
  data$Task_ID <- data$Task_ID %>% as.character()

  ## Get the list of output IDs
  list_ID_output <- data$Output_ID %>% unique()

  ## Get the list of task IDs
  list_ID_task <- data$Task_ID %>% unique()

  if (grid_inputs %>% is.null()) {
    ## Extract the union of all reference inputs provided in the training data
    all_inputs <- data %>%
      dplyr::select(-c(Task_ID, Output)) %>%
      unique() %>%
      tidyr::separate(Reference,
                      into = c("Output_ID_temp", "Input_temp"),
                      sep = ";",
                      remove = FALSE) %>%
      dplyr::mutate(Input_temp_numeric = as.numeric(Input_temp)) %>%
      dplyr::arrange(Output_ID_temp, Input_temp_numeric) %>%
      dplyr::select(c(Input_1, Reference, Output_ID))

    all_input <- all_inputs %>%
      dplyr::pull(Reference)
    cat(
      "The argument 'grid_inputs' is NULL, the hyper-posterior distribution",
      "will only be evaluated on observed Input from 'data'.\n \n"
    )

  } else {
    ## Handle grid_inputs formatting (vector vs tibble)
    if(grid_inputs %>% is.vector()){
      grid_inputs <- tibble::tibble('Input' = grid_inputs,
                                    'Input_ID' = rep("1", length(grid_inputs)),
                                    'Output_ID' = rep(as.factor("1"), length(grid_inputs)))
    }

    if(grid_inputs %>% tibble::is_tibble() & !("Reference" %in% names(grid_inputs))){
      grid_inputs <- grid_inputs %>%
        dplyr::group_by(Input_ID) %>%
        dplyr::mutate(id_ligne = dplyr::row_number()) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider(
          names_from = Input_ID,
          values_from = Input,
          names_prefix = "Input_"
        ) %>%
        dplyr::select(-id_ligne) %>%
        dplyr::mutate(across(starts_with("Input_"), ~ round(.x, 6))) %>%
        rowwise() %>%
        dplyr::mutate(
          Reference = paste(
            paste0("o", Output_ID),
            paste(c_across(starts_with("Input_")), collapse = ":"),
            sep = ";"
          )
        )
    }

    all_inputs <- data %>%
      dplyr::select(Reference, tidyselect::all_of(names_col), Output_ID) %>%
      dplyr::union(grid_inputs) %>%
      unique() %>%
      tidyr::separate(Reference,
                      into = c("Output_ID_temp", "Input_temp"),
                      sep = ";",
                      remove = FALSE) %>%
      dplyr::mutate(Input_temp_numeric = as.numeric(Input_temp)) %>%
      dplyr::arrange(Output_ID_temp, Input_temp_numeric) %>%
      dplyr::select(c(Input_1, Reference, Output_ID))

    all_input <- all_inputs %>% dplyr::pull(Reference)

    if(length(data$Output_ID %>% unique()) == 1){
      all_inputs <- all_inputs %>%
        dplyr::select(-Output_ID)
    }
  }

  ## Initialise m_k according to the value provided by the user
  m_k <- list()

  # Define clusters' names
  names_k <- paste0("K", 1:nb_cluster)

  # if (prior_mean_k %>% is.null()) {
  #   for (k in 1:nb_cluster) {
  #     m_k[[names_k[k]]] <- rep(0, length(all_input))
  #   }
  #   cat(
  #     "The 'prior_mean' argument has not been specified. The hyper_prior mean",
  #     "function is thus set to be 0 everywhere.\n \n"
  #   )
  # } else if (prior_mean_k[[1]] %>% is.function()) {
  #   for (k in 1:nb_cluster) {
  #     all_inputs %>% dplyr::select(- c(Output_ID, Reference))
  #     m_k[[names_k[k]]] <- prior_mean_k[[k]](all_inputs)
  #   }
  # } else if (prior_mean_k %>% is.vector()) {
  #
  #   # Récupération des Outputs triés
  #   unique_outputs_sorted <- list_ID_output %>% unlist() %>% unique() %>% sort()
  #   num_outputs <- length(unique_outputs_sorted) # Nombre d'outputs
  #
  #   if (length(prior_mean_k) == nb_cluster * num_outputs) {
  #
  #     # Extract the prefix of each point of the grid all_input (ex: "o1", "o2")
  #     all_input_prefixes <- stringr::str_extract(all_input, "o[0-9]+")
  #
  #     ## Create a list named by cluster
  #     for (k in 1:nb_cluster) {
  #
  #       all_inputs_k <- all_inputs
  #       start_index <- (k - 1) * num_outputs + 1
  #       end_index   <- k * num_outputs
  #       vals_cluster_k <- prior_mean_k[start_index:end_index]
  #
  #       # Create a corresponding table mapping:
  #       # "o1" -> mean_val_1, "o2" -> mean_val_2
  #       prior_mean_map <- setNames(vals_cluster_k,
  #                                  paste0("o", unique_outputs_sorted))
  #
  #       # Assign it to the correctly named element of the list
  #       # Map values to the full grid based on prefixes
  #       m_k[[names_k[k]]] <- prior_mean_map[all_input_prefixes] %>% unname()
  #       names(m_k[[names_k[k]]]) <- all_inputs$Reference
  #     }
  #
  #   } else if (length(prior_mean_k) == num_outputs) {
  #     # One mean per Output (all clusters share the same)
  #     prior_mean_map <- setNames(prior_mean_k,
  #                                paste0("o", unique_outputs_sorted))
  #     all_input_prefixes <- stringr::str_extract(all_input, "o[0-9]+")
  #
  #     for (k in 1:nb_cluster) {
  #       m_k[[names_k[k]]] <- prior_mean_map[all_input_prefixes] %>% unname()
  #     }
  #   } else {
  #     stop(sprintf("Incorrect length for prior_mean_k. Expected %d or %d, got %d.",
  #                  nb_cluster * num_outputs, num_outputs, length(prior_mean_k)))
  #   }
  # }

  # # Label switching check
  # if (!is.null(prior_mean_k) && !is.null(mixture)) {
  #   # cat("Checking for label switching across multiple outputs (Hyperposterior step)...\n")
  #
  #   # Helper: Get profile from Prior (m_k)
  #   get_prior_profile <- function(cluster_name) {
  #     browser()
  #     vals <- m_k[[cluster_name]]
  #     # We rely on all_input vector for names implicitly if not named,
  #     # but m_k lists are usually just vectors.
  #     # Let's use all_input to map back to Output_IDs
  #     if(is.null(names(vals))) refs <- all_input else refs <- names(vals)
  #
  #     oids <- stringr::str_extract(refs, "^o[0-9]+")
  #
  #     tibble::tibble(Output_ID = oids, value = vals) %>%
  #       dplyr::group_by(Output_ID) %>%
  #       dplyr::summarise(mean_val = mean(value), .groups = "drop") %>%
  #       dplyr::arrange(Output_ID) %>%
  #       dplyr::pull(mean_val)
  #   }
  #
  #   prior_signatures <- do.call(rbind, lapply(names_k, get_prior_profile))
  #   rownames(prior_signatures) <- names_k
  #
  #   # Helper: Get profile from Data + Mixture
  #   # Note: 'data' has columns: Task_ID, Output_ID, Output, Reference...
  #   df_check <- data %>%
  #     dplyr::mutate(OID_label = paste0("o", Output_ID)) %>%
  #     dplyr::select(Task_ID, OID_label, Output) %>%
  #     dplyr::left_join(mixture, by = "Task_ID")
  #
  #   get_empiric_profile <- function(cluster_col) {
  #     # browser() # À retirer une fois le debug fini
  #     col_name_str <- as.character(cluster_col)
  #     # Calcul de la somme des poids pour ce cluster
  #     sum_weights <- sum(df_check[[col_name_str]], na.rm = TRUE)
  #     # Si le cluster est vide (poids quasi nuls), on retourne NA directement
  #     if (sum_weights < 1e-9) {
  #       # On retourne un vecteur de NA de la taille du nombre d'outputs uniques
  #       # (Supposons que OID_label est trié et complet)
  #       n_outputs <- dplyr::n_distinct(df_check$OID_label)
  #       return(rep(NA_real_, n_outputs))
  #     }
  #
  #     df_check %>%
  #       dplyr::group_by(OID_label) %>%
  #       dplyr::summarise(
  #         w_mean = stats::weighted.mean(Output, w = .data[[col_name_str]], na.rm = TRUE),
  #         .groups = "drop"
  #       ) %>%
  #       dplyr::arrange(OID_label) %>%
  #       dplyr::pull(w_mean)
  #   }
  #
  #   # Ensure we look at the clusters present in the mixture
  #   # The mixture columns usually match ID_k
  #   empiric_signatures <- do.call(rbind, lapply(ID_k, get_empiric_profile))
  #   rownames(empiric_signatures) <- ID_k
  #
  #   # Compute distance matrix between Prior (rows) and Empiric (cols of loop, rows of matrix)
  #   dist_mat <- matrix(NA, nrow = nb_cluster, ncol = nb_cluster,
  #                      dimnames = list(ID_k, names_k))
  #
  #   for(i in 1:nb_cluster) { # Loop over Empiric (ID_k)
  #     for(j in 1:nb_cluster) { # Loop over Prior (names_k)
  #       # Check dimensions to avoid crash if priors and data don't span same outputs
  #       if(length(empiric_signatures[i, ]) == length(prior_signatures[j, ])) {
  #         dist_mat[i, j] <- sum((empiric_signatures[i, ] - prior_signatures[j, ])^2)
  #       } else {
  #         dist_mat[i, j] <- Inf # Safety
  #       }
  #     }
  #   }
  #
  #   # 1. On cherche pour CHAQUE PRIOR (colonne) quel est son cluster empirique le plus proche
  #   # Cela nous donne un vecteur de taille nb_cluster (ex: c(1, 1, 2) -> Priors 1 et 2 vont vers Empiric 1)
  #   # Note : which.min ignore les Inf, donc les clusters empiriques vides ne seront pas choisis ici.
  #   assigned_empiric_indices <- apply(dist_mat, 2, which.min)
  #
  #   # Création de la matrice pour stocker les nouvelles moyennes de priors réordonnées
  #   # Elle doit avoir la même structure que prior_signatures
  #   new_prior_means <- matrix(NA, nrow = nrow(prior_signatures), ncol = ncol(prior_signatures))
  #   rownames(new_prior_means) <- ID_k # On garde les IDs des clusters empiriques (K1, K2...)
  #
  #   # Liste pour suivre les priors qui ont été utilisés
  #   used_priors <- numeric()
  #
  #   # 2. Boucle sur les clusters EMPIRIQUES (i) pour construire leur nouveau prior
  #   for(i in 1:nb_cluster) {
  #
  #     # Quels priors ont "choisi" ce cluster empirique i ?
  #     priors_merged_here <- which(assigned_empiric_indices == i)
  #
  #     if(length(priors_merged_here) > 1) {
  #       # --- CAS DE FUSION ---
  #       # Plusieurs priors pointent vers ce cluster empirique.
  #       # On prend la moyenne arithmétique de leurs signatures (Moyenne des vecteurs de moyennes)
  #
  #       # colMeans gère correctement la moyenne par colonne (par Output)
  #       new_prior_means[i, ] <- colMeans(prior_signatures[priors_merged_here, , drop = FALSE])
  #
  #       # On note ces priors comme utilisés
  #       used_priors <- c(used_priors, priors_merged_here)
  #
  #     } else if(length(priors_merged_here) == 1) {
  #       # --- CAS 1-pour-1 ---
  #       # Un seul prior pointe ici, on le copie simplement
  #       new_prior_means[i, ] <- prior_signatures[priors_merged_here, ]
  #       used_priors <- c(used_priors, priors_merged_here)
  #
  #     } else {
  #       # --- CAS CLUSTER EMPIRIQUE VIDE ---
  #       # Aucun prior n'a choisi ce cluster (souvent car sa ligne est Inf ou très loin)
  #       # On laisse NA pour l'instant, on remplira avec les priors orphelins juste après
  #     }
  #   }
  #
  #   # 3. Gestion des Orphelins (Priors non assignés) et des Trous (Clusters vides)
  #   # Si des priors ont fusionné ailleurs, cela laisse des clusters empiriques vides.
  #   # Il faut réassigner les priors "perdus" (ceux qui n'ont pas été choisis du tout)
  #   # aux clusters vides pour ne pas perdre de diversité (mécanisme de sauvetage).
  #
  #   all_priors <- 1:nb_cluster
  #   orphan_priors <- setdiff(all_priors, used_priors)
  #   empty_empiric_slots <- which(is.na(new_prior_means[, 1])) # On regarde la 1ère colonne pour voir les NA
  #
  #   # Si on a des trous et des orphelins, on remplit
  #   if(length(empty_empiric_slots) > 0 && length(orphan_priors) > 0) {
  #
  #     # On peut faire un matching simple ou aléatoire ici.
  #     # Prenons les dans l'ordre pour simplifier.
  #     n_fill <- min(length(empty_empiric_slots), length(orphan_priors))
  #
  #     for(k in 1:n_fill) {
  #       slot_idx <- empty_empiric_slots[k]
  #       prior_idx <- orphan_priors[k]
  #       new_prior_means[slot_idx, ] <- prior_signatures[prior_idx, ]
  #     }
  #   }
  #
  #   # S'il reste encore des NA (ex: plus de clusters vides que de priors dispo - rare),
  #   # on peut réinitialiser ou dupliquer un existant.
  #   # Ici, sécurité simple : remplacer les NA restants par 0 ou la moyenne globale.
  #   new_prior_means[is.na(new_prior_means)] <- 0
  #
  #   # Mise à jour des noms pour la suite du code (si besoin)
  #   # Note : new_prior_means contient maintenant les valeurs numériques directement
  #   # Si votre code attend une liste de noms de clusters, la logique change légèrement,
  #   # mais pour mettre à jour les hyper-paramètres, c'est cette matrice `new_prior_means` qu'il faut utiliser.
  # }

  ## Create a list named by cluster with evaluation of the prior mean (m_k) at all Input locations
  for (k in 1:nb_cluster) {
    # Extrait le numéro après le 'o' (ex: "o1" -> 1, "o2" -> 2)
    out_indices <- as.numeric(gsub("^o([0-9]+);.*", "\\1", all_input))

    # Assigne directement en utilisant ces indices pour piocher dans prior_mean_k
    m_k[[ID_k[k]]] <- prior_mean_k[[ID_k[k]]][out_indices]
  }

  ## Format a sequence of inputs for all clusters
  t_clust <- tidyr::expand_grid("Cluster_ID" = names(m_k),
                                all_inputs
  )

  list_inv_k <- list()
  list_inv_t <- list()

  # Loop over clusters
  for(k in t_clust$Cluster_ID %>% unique){
    t_clust_k <- t_clust %>%
      dplyr::filter(Cluster_ID == k) %>%
      dplyr::select(-Cluster_ID)

    hp_k_subset <- hp_k %>%
      dplyr::filter(Cluster_ID == k) %>%
      dplyr::select(-prop_mixture)

    if(length(list_ID_output) > 1 && !(kern_t %>% is.character())){
      cov_k <- kern_to_cov(input = t_clust_k,
                           kern = kern_k,
                           hp = hp_k_subset %>% dplyr::select(-Cluster_ID))
    } else {
      cov_k <- kern_to_cov(input = t_clust_k,
                           kern = kern_k,
                           hp = hp_k_subset %>%
                             dplyr::select(-c(Cluster_ID, Output_ID))
      )
    }

    references <- rownames(cov_k)
    inv_k <- cov_k %>% chol_inv_jitter(pen_diag = pen_diag)
    dimnames(inv_k) <- list(references, references)
    list_inv_k[[k]] <- inv_k
  }

  # Loop over tasks
  for (t in list_ID_task) {
    data_t <- data %>% dplyr::filter(Task_ID == t) %>%
      dplyr::select(-c(Output, Task_ID))
    hp_t_indiv <- hp_t %>% dplyr::filter(Task_ID == t)

    if(length(list_ID_output) > 1 && !(kern_t %>% is.character())){
      K_task_t <- kern_to_cov(input = data_t, kern = kern_t, hp = hp_t_indiv)
    } else{
      all_inputs_t <- data %>%
        dplyr::filter(Task_ID == t) %>%
        dplyr::select(-c(Task_ID, Output, Output_ID)) %>%
        unique()
      K_task_t <- kern_to_cov(input = all_inputs_t, kern = kern_t, hp = hp_t_indiv)
    }

    task_references <- rownames(K_task_t)
    K_inv_t <- K_task_t %>% chol_inv_jitter(pen_diag = pen_diag)
    dimnames(K_inv_t) <- list(task_references, task_references)
    list_inv_t[[t]] <- K_inv_t
  }

  ## Create a named list of Output values for all tasks
  list_output_t <- base::split(data$Output, list(data$Task_ID))

  ## Update each mu_k parameters for each cluster
  floop <- function(k) {
    post_inv <- list_inv_k[[k]]
    tau_k <- mixture %>% dplyr::select(Task_ID, k)

    for (t in list_inv_t %>% names()) {
      tau_t_k <- tau_k %>% dplyr::filter(Task_ID == t) %>% dplyr::pull(k)
      inv_t <- list_inv_t[[t]]
      co_input <- intersect(row.names(inv_t), row.names(post_inv))
      post_inv[co_input, co_input] <- post_inv[co_input, co_input] +
        tau_t_k * inv_t[co_input, co_input]
    }
    post_inv %>%
      chol_inv_jitter(pen_diag = pen_diag) %>%
      `rownames<-`(all_input) %>%
      `colnames<-`(all_input) %>%
      return()
  }
  cov_k <- sapply(tidyselect::all_of(names(m_k)),
                  floop,
                  simplify = FALSE,
                  USE.NAMES = TRUE)

  ## Update the posterior mean for each cluster
  floop2 <- function(k) {
    weighted_k <- list_inv_k[[k]] %*% m_k[[k]]

    for (t in list_inv_t %>% names()) {
      tau_t_k <- mixture %>% dplyr::filter(Task_ID == t) %>% dplyr::pull(k)
      weighted_t <- tau_t_k * list_inv_t[[t]] %*% list_output_t[[t]]
      co_input <- intersect(row.names(weighted_t), row.names(weighted_k))
      weighted_k[co_input, ] <- weighted_k[co_input, ] +
        weighted_t[co_input, ]
    }

    post_mean <- cov_k[[k]] %*% weighted_k %>% as.vector()
    tibble::tibble(all_inputs, "Output" = post_mean) %>% return()
  }
  mean_k <- sapply(names(m_k), floop2, simplify = FALSE, USE.NAMES = TRUE)

  ## Format the GP prediction of the hyper-posterior mean (for direct plot)
  floop_pred <- function(k) {
    tibble::tibble(mean_k[[k]],
                   "Var" = cov_k[[k]] %>% diag() %>% as.vector()
    ) %>%
      dplyr::rename("Mean" = Output) %>%
      dplyr::select(-Reference) %>%
      return()
  }
  pred <- sapply(ID_k, floop_pred, simplify = FALSE, USE.NAMES = TRUE)

  list("mean" = mean_k, "cov" = cov_k, "mixture" = mixture, "pred" = pred) %>%
    return()
}

#' MagmaClust prediction
#'
#' Compute the posterior predictive distribution in MagmaClust.
#' Providing data from any new task, its trained hyper-parameters
#' and a previously trained MagmaClust model, the multi-task posterior
#' distribution is evaluated on any arbitrary inputs that are specified through
#' the 'grid_inputs' argument. Due to the nature of the model, the prediction is
#' defined as a mixture of Gaussian distributions. Therefore the present
#' function computes the parameters of the predictive distribution
#' associated with each cluster, as well as the posterior mixture probabilities
#' for this new task.
#'
#' @param data  A tibble or data frame. Required columns: \code{Input},
#'    \code{Output}. Additional columns for covariates can be specified.
#'    The \code{Input} column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    \code{Output} column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference 'Input'. If NULL, the mixture of mean processes from
#'    \code{trained_model} is returned as a generic prediction.
#' @param trained_model A list, containing  the information coming from a
#'    MagmaClust model, previously trained using the
#'    \code{\link{train_magmaclust}} function. If \code{trained_model} is set to
#'    NULL, the \code{hyperpost} and \code{prop_mixture} arguments are mandatory
#'    to perform required re-computations for the prediction to succeed.
#' @param grid_inputs The grid of inputs (reference Input and covariates) values
#'    on which the GP should be evaluated. Ideally, this argument should be a
#'    tibble or a data frame, providing the same columns as \code{data}, except
#'    'Output'. Nonetheless, in cases where \code{data} provides only one
#'    'Input' column, the \code{grid_inputs} argument can be NULL (default) or a
#'    vector. This vector would be used as reference input for prediction and if
#'    NULL, a vector of length 500 is defined, ranging between the min and max
#'    Input values of \code{data}.
#' @param mixture A tibble or data frame, indicating the mixture probabilities
#'     of each cluster for the new task.
#'     If NULL, the \code{\link{train_gp_clust}} function is used to compute
#'     these posterior probabilities according to \code{data}.
#' @param hp A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern}. The columns/elements should be named
#'    according to the hyper-parameters that are used in \code{kern}. The
#'    \code{\link{train_gp_clust}} function can be used to learn
#'    maximum-likelihood estimators of the hyper-parameters.
#' @param kern A kernel function, defining the covariance structure of the GP.
#'    Several popular kernels
#'    (see \href{https://www.cs.toronto.edu/~duvenaud/cookbook/}{The Kernel
#'    Cookbook}) are already implemented and can be selected within the
#'    following list:
#'    - "SE": (default value) the Squared Exponential Kernel (also called
#'        Radial Basis Function or Gaussian kernel),
#'    - "LIN": the Linear kernel,
#'    - "PERIO": the Periodic kernel,
#'    - "RQ": the Rational Quadratic kernel,
#'    - convolution_kernel: the Convolution kernel used to manage MO scenario.
#'    Compound kernels can be created as sums or products of the above kernels.
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#' @param hyperpost A list, containing the elements \code{mean}, \code{cov} and
#'   \code{mixture} the parameters of the hyper-posterior distributions of the
#'    mean processes. Typically, this argument should come from a previous
#'    learning using \code{\link{train_magmaclust}}, or a previous prediction
#'    with \code{\link{pred_magmaclust}}, with the argument \code{get_hyperpost}
#'    set to TRUE.
#' @param prop_mixture A tibble or a named vector of the mixture proportions.
#'    Each name of column or element should refer to a cluster. The value
#'    associated with each cluster is a number between 0 and 1. If both
#'    \code{mixture} and \code{trained_model} are set to NULL, this argument
#'    allows to recompute mixture probabilities, thanks to the \code{hyperpost}
#'    argument and the \code{\link{train_gp_clust}} function.
#' @param get_hyperpost A logical value, indicating whether the hyper-posterior
#'    distributions of the mean processes should be returned. This can be useful
#'    when planning to perform several predictions on the same grid of inputs,
#'    since recomputation of the hyper-posterior can be prohibitive for high
#'    dimensional grids.
#' @param get_full_cov A logical value, indicating whether the full posterior
#'    covariance matrices should be returned.
#' @param plot A logical value, indicating whether a plot of the results is
#'    automatically displayed.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#'
#' @return A list of GP prediction results composed of:
#'         \itemize{
#'           \item pred: As sub-list containing, for each cluster:
#'             \itemize{
#'              \item pred_gp: A tibble, representing the GP predictions as two
#'                 column \code{Mean} and \code{Var}, evaluated on the
#'                 \code{grid_inputs}. The column \code{Input} and additional
#'                 covariates columns are associated with each predicted values.
#'              \item proba: A number, the posterior probability associated with
#'                 this cluster.
#'              \item cov (if \code{get_full_cov} = TRUE): A matrix, the full
#'                 posterior covariance matrix associated with this cluster.
#'              }
#'           \item mixture: A tibble, indicating the mixture probabilities
#'              of each cluster for the predicted task.
#'           \item hyperpost (if \code{get_hyperpost} = TRUE): A list,
#'              containing the hyper-posterior distributions information useful
#'              for visualisation purposes.
#'          }
#'
#' @export
#'
#' @examples
#' TRUE
pred_magmaclust <- function(data = NULL,
                            trained_model = NULL,
                            grid_inputs = NULL,
                            mixture = NULL,
                            hp = NULL,
                            kern = "SE",
                            hyperpost = NULL,
                            prop_mixture = NULL,
                            get_hyperpost = FALSE,
                            get_full_cov = TRUE,
                            plot = TRUE,
                            pen_diag = 1e-10) {
  ## Return the mean process if no data is provided
  if(data %>% is.null()){
    ## Check whether trained_model is provided
    if(trained_model %>% is.null()){
      stop(
        "If 'data' is not provided, the 'trained_model' argument is needed ",
        "to provide the mixture of mean processes as a generic prediction."
      )
    }
    else{
      ## Return the mean process as a generic prediction
      if(grid_inputs %>% is.null()){
        hyperpost = trained_model$hyperpost

      } else{
        ## Recompute the mean process at the required inputs if necessary
        hyperpost = hyperposterior_clust(
          trained_model = trained_model,
          grid_inputs = grid_inputs,
          kern_k = trained_model$ini_args$kern_k,
          kern_t = trained_model$ini_args$kern_t,
          hp_k = trained_model$hp_k,
          hp_t = trained_model$hp_t,
        )

        ## Retain only grid_inputs for display purposes
        for(k in names(hyperpost$pred))
        {
          ## Retain only grid_inputs for display purposes
          grid_inputs_wide <- grid_inputs %>%
            dplyr::group_by(Input_ID) %>%
            dplyr::mutate(id_ligne = dplyr::row_number()) %>%
            dplyr::ungroup() %>%
            tidyr::pivot_wider(
              names_from = Input_ID,
              values_from = Input,
              names_prefix = "Input_"
            )

          shared_columns <- intersect(
            names(grid_inputs_wide),
            names(hyperpost$pred[[k]])
          )

          hyperpost$pred[[k]] <- grid_inputs_wide %>%
            dplyr::inner_join(hyperpost$pred[[k]], by = shared_columns) %>%
            dplyr::select(-id_ligne)
        }
      }
    }
      names_k = hyperpost$pred %>% names()

      ## Compute the generic mixture weights
      mixture = hyperpost$mixture %>%
        dplyr::select(-Task_ID) %>%
        dplyr::summarise(dplyr::across(tidyselect::everything(), mean)) %>%
        dplyr::mutate('Task_ID' = 'Task_ID_pred', .before = 1)

      ## Add the ID and Proba columns to 'pred'
      floop_k = function(k){
        hyperpost$pred[[k]] %>%
          dplyr::mutate(
            'Task_ID' = 'Task_ID_pred',
            'Proba' = mixture[[k]],
            .before = 1
          ) %>%
          return()
      }
      pred_k = sapply(names_k, floop_k, simplify = FALSE, USE.NAMES = TRUE)

      ## Compute the mixture mean and variance of predictions
      mixture_mean <- 0
      for (k in names_k)
      {
        proba_k <- mixture %>% dplyr::pull(k)
        mixture_mean = mixture_mean + proba_k * pred_k[[k]]$Mean
      }

      mixture_var <- 0
      for (k in names_k)
      {
        proba_k <- mixture %>% dplyr::pull(k)
        ## Cov(mixture) = Sum_k{ tau_k * (C_k + (m_k - m)(m_k - m)T) }
        mixture_var <- mixture_var +
          proba_k * (pred_k[[k]]$Var + (pred_k[[k]]$Mean - mixture_mean)^2)
      }

      ## Create the mixture of GPs prediction
      mixture_pred <- pred_k[[names_k[1]]] %>%
        dplyr::mutate(
        "Mean" = mixture_mean,
        "Var" = mixture_var %>% as.vector()
        ) %>%
        dplyr::select(-Proba)

      ## Define the list to return with the correct format
      pred <- list(
        "pred" = pred_k,
        'mixture' = mixture,
        'mixture_pred' = mixture_pred)

      ## Display the graph of the prediction if expected
      if (plot) {
        data_train <- trained_model$ini_args$data

        ## Add 'cov' to display samples
        pred[["cov"]] <- hyperpost$cov

        ## Display samples only in 1D and Credible Interval otherwise
        if(ncol(mixture_pred) == 4){
          ## Plot 1D prediction
          plot_magmaclust(
            pred,
            data = data,
            data_train = data_train,
            prior_mean = hyperpost$mean,
            samples = TRUE
          ) %>%
            print()
        } else {
          ## Plot 2D prediction

          plot_magmaclust(
            pred,
            samples = FALSE
          ) %>%
            print()
        }

      ## Check whether posterior covariance should be returned
      if (!get_full_cov) {
        pred[["cov"]] <- NULL
      }
    }
    return(pred)
  }

  ## Add an 'ID' column if present
  if ("Task_ID" %in% names(data)) {
    if (dplyr::n_distinct(data$Task_ID) > 1) {
      stop(
        "Problem in the 'Task_ID' column: different values are not allowed. ",
        "The prediction can only be performed for one task."
      )
    }

    ## Get 'Task_ID' of the task to predict
    ID_task_pred <- unique(data$Task_ID)

  } else {
    ## Set 'Task_ID' of the task to predict to 'Task_ID_pred' if not provided
    ID_task_pred <- "Task_ID_pred"
    data <- data %>% dplyr::mutate("Task_ID" = "Task_ID_pred", .before = 1)
  }

  if("Reference" %in% names(data)){
    ## Get input column names
    names_col <- data %>%
      dplyr::select(- c(Output, Output_ID, Reference)) %>%
      names()
  } else{
    names_col <- data %>%
      dplyr::select(- c(Output, Output_ID)) %>%
      names()
  }

  data <- data %>%
    dplyr::group_by(Output_ID, Output, Input_ID) %>%
    # Add a unique number observation for the group
    dplyr::mutate(obs_num = row_number()) %>%
    dplyr::ungroup()

  ## To create the 'Reference' column as in the old MagmaClustR tibble format, we
  # need to pivot data to obtain one row per observation of (Task_ID, Output_ID).
  # In other words, inputs are no longer in "short" format; instead, we have one
  # column per input.
  data <- data %>%
    tidyr::pivot_wider(
      names_from = Input_ID,
      values_from = Input,
      names_prefix = "Input_"
    ) %>%
    # Keep 6 significant digits for Inputs to avoid numerical issues
    dplyr::mutate(across(starts_with("Input_"), ~ round(.x, 6))) %>%
    rowwise() %>%
    dplyr::mutate(
      Reference = paste(
        # Create output's prefix
        paste0("o", Output_ID),
        # Create the reference for each Output_ID
        paste(c_across(starts_with("Input_")), collapse = ":"),
        # Join output's prefix and reference
        sep = ";"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-obs_num)

  ## Extract the observed Output (data points)
  data_obs <- data %>%
    dplyr::pull(Output)

  ## Extract the Reference
  input_obs <- data %>%
    dplyr::pull(Reference)

  ## Extract the observed inputs (reference Input + covariates)
  inputs_obs <- data %>%
    dplyr::select(-c(Output, Task_ID))

  ## Define the target inputs to predict
  if (grid_inputs %>% is.null()) {
    set_grid <- function(data, size_grid) {
      seq(data %>% min(),
          data %>% max(),
          length.out = size_grid
      ) %>%
        return()
    }

    if (inputs_obs %>% names() %>% length() == 3) {
      size_grid <- 500
    } else if (inputs_obs %>% names() %>% length() > 3) {
      size_grid <- 1000^(1 / (ncol(inputs_obs) - 1)) %>% round()
      ## floor instead of round ?
    } else {
      stop(
        "The 'grid_inputs' argument should be a either a numerical vector ",
        "or a data frame depending on the context. Please read ?pred_magma()."
      )
    }

    # Create a unique grid of inputs (which will be replicated for each Output_ID)
    base_grid <- purrr::map(
      data %>% dplyr::select(tidyselect::all_of(names_col)),
      set_grid,
      size_grid
    ) %>%
      purrr::set_names(paste0("Input_", seq_along(.))) %>%
      expand.grid() %>%
      tibble::as_tibble()

    unique_outputs <- data %>%
      dplyr::distinct(Output_ID)

    # Cross base_grid with unique_outputs to replicate base_grid as many times as
    # unique_outputs and create 'Reference' column
    inputs_pred <- tidyr::crossing(unique_outputs, base_grid) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        Reference = paste(
          paste0("o", Output_ID),
          paste(c_across(starts_with("Input_")), collapse = ":"),
          sep = ";"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(Reference)

    input_pred <- inputs_pred %>% dplyr::pull(Reference)

  } else if (grid_inputs %>% is.vector()) {
    ## Test whether 'data' only provide the Input column and no covariates
    if (inputs_obs %>% names() %>% length() == 3) {
      input_temp <- grid_inputs %>%
        signif() %>%
        sort() %>%
        unique()

      inputs_pred <- tibble::tibble("Input_1" = rep(input_temp,
                                                    times = nrow(unique_outputs)),
                                    "Output_ID" = rep(unique_outputs %>%
                                                        purrr::as_vector(),
                                                      each = length(input_temp))) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          Reference = paste(
            paste0("o", Output_ID),
            paste(c_across(starts_with("Input_")), collapse = ":"),
            sep = ";"
          )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(Reference)

      input_pred <- inputs_pred %>% dplyr::pull(Reference)

    } else {
      stop(
        "The 'grid_inputs' argument should be a either a numerical vector ",
        "or a data frame depending on the context. Please read ?pred_magma()."
      )
    }
  } else if (grid_inputs %>% is.data.frame()) {
    grid_inputs <- grid_inputs %>%
      dplyr::group_by(Input_ID) %>%
      dplyr::mutate(id_ligne = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(
        names_from = Input_ID,
        values_from = Input,
        names_prefix = "Input_"
      ) %>%
      dplyr::select(-id_ligne) %>%
      # Keep 6 significant digits for Inputs to avoid numerical issues
      dplyr::mutate(across(starts_with("Input_"), ~ round(.x, 6))) %>%
      rowwise() %>%
      dplyr::mutate(
        Reference = paste(
          # Create output's prefix
          paste0("o", Output_ID),
          # Create the reference for each Output_ID
          paste(c_across(starts_with("Input_")), collapse = ":"),
          # Join output's prefix and reference
          sep = ";"
        )
      )

    ## Test whether 'data' has the same columns as grid_inputs
    if (names(inputs_obs) %>% setequal(names(grid_inputs))) {
      input_pred <- grid_inputs %>%
        # DO NOT arrange Reference because of the lexicographic order
        # (not armful but unnecessary in Magma case, tragic in MO case)
        dplyr::pull(Reference)

      inputs_pred <- grid_inputs %>%
        # DO NOT arrange Reference because of the lexicographic order
        # (not armful but unnecessary in Magma case, tragic in MO case)
        dplyr::select(names(inputs_obs))
    } else {
      stop(
        "The 'grid_inputs' argument should provide a column 'Input_ID', 'Input' ",
        " and 'Output_ID'."
      )
    }
  } else {
    stop(
      "The 'grid_inputs' argument should be a either a numerical vector ",
      "or a data frame depending on the context. Please read ?pred_gp()."
    )
  }

  ## Define the union of all distinct reference Input
  all_inputs <- dplyr::union(inputs_obs, inputs_pred)
  all_inputs <- all_inputs %>%
    arrange(Output_ID, Input_1)

  # DO NOT arrange Reference because of the lexicographic order
  # (not armful but unnecessary in Magma case, tragic in MO case)
  all_input <- all_inputs$Reference

  ## Check whether the hyper-posterior is provided and recompute if necessary
  if (hyperpost %>% is.null()) {
    if (trained_model %>% is.null()) {
      stop(
        "If the 'hyperpost' argument is NULL, the 'trained_model' ",
        "should be provided, in order to extract or recompute mean processes' ",
        "hyper-posterior distributions evaluated at the correct inputs."
      )
    } else {
      ## Get the hyper-posterior distribution from the trained model
      hyperpost <- trained_model$hyperpost

      ## Check hyperpost format (in particular presence of all reference Input)
      if (!all(all_input %in% hyperpost$mean[[1]]$Reference)) {
        if (trained_model %>% is.null()) {
          stop(
            "hyperpost is not evaluated at the correct inputs, please use the ",
            "'trained_model' argument instead."
          )
        }
        cat(
          "The hyper-posterior distribution of the mean process provided in",
          "'hyperpost' argument isn't evaluated on the expected inputs.",
          "Start evaluating the hyper-posterior on the correct inputs...\n \n"
        )

        hyperpost <- hyperposterior_clust(
          data = trained_model$ini_args$data,
          mixture = trained_model$hyperpost$mixture,
          hp_k = trained_model$hp_k,
          hp_t = trained_model$hp_t,
          kern_k = trained_model$ini_args$kern_k,
          kern_t = trained_model$ini_args$kern_t,
          prior_mean_k = trained_model$m_k,
          # prior_mean_k = trained_model$ini_args$prior_mean_k,
          grid_inputs = all_inputs,
          pen_diag = pen_diag
        )
        cat("Done!\n \n")
      }
    }
  } else if (hyperpost %>% is.list()) {
    ## Check hyperpost format
    if (!is.null(hyperpost$mean)) {
      ## Check hyperpost format (in particular presence of all reference Input)
      if (!all(all_input %in% hyperpost$mean[[1]]$Reference)) {
        stop(
          "The hyper-posterior distribution of the mean processes provided ",
          "in the 'hyperpost' argument isn't evaluated on expected inputs. ",
          "Please provide a 'trained_model' argument for re-computation. "
        )
      }
    } else {
      stop(
        "The format of the 'hyperpost' argument is not as expected. Please ",
        "read ?pred_magmaclust() for details."
      )
    }
  } else {
    stop(
      "The format of the 'hyperpost' argument is not as expected. Please ",
      "read ?pred_magmaclust() for details."
    )
  }

  ## Get clusters' names
  ID_k <- hyperpost$mean %>% names()

  ## Extract or learn the hyper-parameters if not provided
  if (hp %>% is.null()) {
    if (!is.null(trained_model)) {
      ## Check whether hyper-parameters are shared between tasks if we have
      ## 'trained_model'
      if (tryCatch(trained_model$ini_args$shared_hp_tasks,
                   error = function(e) FALSE
      )) {
        ## Extract the hyper-parameters shared netween all 't'
        hp <- trained_model$hp_t %>%
          dplyr::slice(1:length(data$Output_ID %>% unique())) %>%
          dplyr::mutate("Task_ID" = ID_task_pred)
        ## Extract and format the mixture proportions
        prop_mixture <- trained_model$hp_k %>%
          dplyr::select(Cluster_ID, prop_mixture) %>%
          dplyr::distinct() %>%
          dplyr::pull(prop_mixture, name = Cluster_ID)

        mixture <- update_mixture(
          data,
          hyperpost$mean,
          hyperpost$cov,
          hp,
          trained_model$ini_args$kern_t,
          prop_mixture,
          pen_diag
        )
      } else if (kern %>% is.function()) {
        ## Extract the mixture proportions
        prop_mixture <- trained_model$hp_k %>%
          dplyr::distinct(Cluster_ID, prop_mixture) %>%
          dplyr::pull(prop_mixture, name = Cluster_ID)

        cat(
          "The 'hp' argument has not been specified. The 'train_gp_clust()'",
          "function (with random initialisation) will be used to learn ML",
          "estimators for hyper-parameters and mixture probabilities... \n \n"
        )

        ini_hp <- hp(kern = convolution_kernel,
                    noise = T,
                    list_task_ID = ID_task_pred,
                    list_output_ID = data$Output_ID %>% unique(),
                    shared_hp_outputs = FALSE,
                    shared_hp_tasks = FALSE)

        hp_mix <- train_gp_clust(
          data,
          prop_mixture = prop_mixture,
          ini_hp = ini_hp,
          kern = kern,
          hyperpost = hyperpost,
          pen_diag = pen_diag
        )

        ## Extract values of hyper-parameters and mixture probabilities
        hp <- hp_mix$hp
        mixture <- hp_mix$mixture

      } else if (kern %>% is.character()) {
        ## Extract the mixture proportions
        prop_mixture <- trained_model$hp_k %>%
          dplyr::pull(prop_mixture, name = Cluster_ID)

        cat(
          "The 'hp' argument has not been specified. The 'train_gp_clust()'",
          "function (with random initialisation) will be used to learn ML",
          "estimators for hyper-parameters and mixture probabilities... \n \n"
        )

        ini_hp <- hp(kern,
                      noise = T,
                      list_task_ID = ID_task_pred,
                      list_output_ID = data$Output_ID %>% unique(),
                      shared_hp_outputs = FALSE,
                      shared_hp_tasks = FALSE) %>%
          dplyr::filter(Task_ID == ID_task_pred)

        hp_mix <- train_gp_clust(
          data,
          prop_mixture = prop_mixture,
          ini_hp = ini_hp,
          kern = kern,
          hyperpost = hyperpost,
          pen_diag = pen_diag
        )

        ## Extract values of hyper-parameters and mixture probabilities
        hp <- hp_mix$hp
        mixture <- hp_mix$mixture
      }
    } else if (kern %>% is.function()) {
      stop(
        "When using a custom kernel function the 'hp' argument is ",
        "mandatory, in order to provide the name of the hyper-parameters. ",
        "You can use the function 'hp()' to easily generate a tibble of random",
        " hyper-parameters with the desired format, or use 'train_gp_clust()' ",
        "to learn ML estimators for a better fit."
      )
    } else if (kern %>% is.character()) {
      cat(
        "The 'hp' argument has not been specified. The 'train_gp_clust()'",
        "function (with random initialisation) will be used to learn ML",
        "estimators for hyper-parameters and mixture probabilities... \n \n"
      )
      ## If 'prop_mixture' has no name, use the ID_k to retrieve them
      if (!("ID" %in% prop_mixture) & !is.null(prop_mixture)) {
        names(prop_mixture) <- ID_k
      }

      ini_hp <- hp(kern,
                   noise = T,
                   list_task_ID = ID_task_pred,
                   list_output_ID = data$Output_ID %>% unique(),
                   shared_hp_outputs = FALSE,
                   shared_hp_tasks = FALSE) %>%
        filter(Task_ID == ID_task_pred)

      hp_mix <- train_gp_clust(
        data,
        prop_mixture = prop_mixture,
        ini_hp = ini_hp,
        kern = kern,
        hyperpost = hyperpost,
        pen_diag = pen_diag
      )

      ## Extract values of hyper-parameters and mixture probabilities
      hp <- hp_mix$hp
      mixture <- hp_mix$mixture
    } else {
      stop(
        "Incorrect format for the 'kern' argument. Please read ?pred_gp() for ",
        "details."
      )
    }
  }

  ## Remove the noise of the hp for evaluating some of the sub-matrix
  if ("noise" %in% names(hp)) {
    hp_rm_noi <- hp %>% dplyr::select(-noise)
    noise <- exp(hp[["noise"]])
  } else {
    hp_rm_noi <- hp
    noise <- 0
  }

  ## Initialisation if we want to recover full_cov
  full_cov <- list()
  ## Initialisation of the mixture prediction
  mixture_mean <- 0

  ponderation_matrix <- NULL
  floop <- function(k) {
    ## Extract the mean parameter from the hyper-posterior
    mean_obs <- hyperpost$mean[[k]] %>%
      dplyr::filter(Reference %in% input_obs) %>%
      # DO NOT arrange Reference because of the lexicographic order
      # (not armful but unnecessary in Magma case, tragic in MO case)
      dplyr::pull(Output)

    mean_pred <- hyperpost$mean[[k]] %>%
      dplyr::filter(Reference %in% input_pred) %>%
      # DO NOT arrange Reference because of the lexicographic order
      # (not armful but unnecessary in Magma case, tragic in MO case)
      dplyr::pull(Output)

    ## Extract the covariance sub-matrices from the hyper-posterior
    post_cov_obs <- hyperpost$cov[[k]][
      as.character(input_obs),
      as.character(input_obs)
    ]
    post_cov_pred <- hyperpost$cov[[k]][
      as.character(input_pred),
      as.character(input_pred)
    ]
    post_cov_crossed <- hyperpost$cov[[k]][
      as.character(input_obs),
      as.character(input_pred)
    ]

    if(length(inputs_obs$Output_ID %>% unique()) == 1){
      inputs_obs <- inputs_obs %>%
        dplyr::select(-Output_ID)
      # DO NOT arrange Reference because of the lexicographic order
      # (not armful but unnecessary in Magma case, tragic in MO case)

      inputs_pred <- inputs_pred %>%
        dplyr::select(-Output_ID)
      # DO NOT arrange Reference because of the lexicographic order
      # (not armful but unnecessary in Magma case, tragic in MO case)
    }

    ## Sum the covariance matrices on observed inputs and compute the inverse
    cov_obs <- kern_to_cov(inputs_obs, kern, hp) + post_cov_obs

    inv_obs <- cov_obs %>%
      chol_inv_jitter(pen_diag = pen_diag) %>%
      `rownames<-`(as.character(input_obs)) %>%
      `colnames<-`(as.character(input_obs))

    ## Compute the other required sum of sub-matrices for prediction
    cov_pred <- kern_to_cov(inputs_pred, kern, hp_rm_noi) + post_cov_pred
    cov_crossed <- kern_to_cov(
      inputs_obs,
      kern,
      hp_rm_noi,
      input_2 = inputs_pred
    ) + post_cov_crossed

    ## Compute the posterior mean of a GP
    pred_mean <- (mean_pred +
                    t(cov_crossed) %*% inv_obs %*% (data_obs - mean_obs)) %>%
      as.vector()

    ## Compute the posterior covariance matrix of a GP
    pred_cov <- cov_pred - t(cov_crossed) %*% inv_obs %*% cov_crossed

    if(mixture[k] == 1){
      ponderation_matrix <<- t(cov_crossed) %*% inv_obs
    }

    if(length(all_inputs$Output_ID %>% unique()) > 1){
      ## Get the pred_cov rownames and colnames to create a noise vector with correct
      ## dimensions
      point_names <- rownames(pred_cov)
      output_ids <- as.numeric(sub("^o(\\d+);.*", "\\1", point_names))
      full_noise_vector <- noise[output_ids]
    } else {
      full_noise_vector <- noise
    }

    # browser()

    ## Keep track of the full predicted covariances
    full_cov[[k]] <<- pred_cov + full_noise_vector

    ## Select the adequate task if necessary
    proba <- mixture %>%
      dplyr::filter(Task_ID == ID_task_pred) %>%
      dplyr::pull(k)

    ## Combine cluster-specific predictions into a mixture prediction
    mixture_mean <<- mixture_mean + proba * pred_mean

    ## Create a tibble of values and associated uncertainty from a GP prediction
    if(length(grid_inputs$Output_ID %>% unique()) == 1){
      inputs_pred$Output_ID <- as.factor("1")
    }

    tibble::tibble(
      "Task_ID" = ID_task_pred,
      "Proba" = proba,
      "Mean" = pred_mean,
      "Var" = (diag(pred_cov) + full_noise_vector) %>% as.vector()
    ) %>%
      dplyr::mutate(inputs_pred) %>%
      dplyr::select(-Reference) %>%
      return()
  }
  pred <- sapply(ID_k, floop, simplify = FALSE, USE.NAMES = TRUE)

  ## Compute the mixture variance of predictions
  mixture_var <- 0
  for (k in ID_k)
  {
    proba_k <- mixture %>%
      dplyr::filter(Task_ID == ID_task_pred) %>%
      dplyr::pull(k)
    ## Cov(mixture) = Sum_k{ tau_k * (C_k + (m_k - m)(m_k - m)T) }
    mixture_var <- mixture_var +
      proba_k * (pred[[k]]$Var + (pred[[k]]$Mean - mixture_mean)^2)
  }

  ## Create a tibble of values for the mixture prediction
  mixture_pred <- tibble::tibble(
    "Task_ID" = ID_task_pred,
    "Mean" = mixture_mean,
    "Var" = mixture_var %>% as.vector()
  ) %>%
    dplyr::mutate(inputs_pred) %>%
    dplyr::select(-Reference)

  res <- list("pred" = pred, "mixture" = mixture, "mixture_pred" = mixture_pred,
              "ponderation_matrix" = ponderation_matrix)

  ## Check whether hyper-posterior should be returned
  if (get_hyperpost) {
    res[["hyperpost"]] <- hyperpost
  }

  ## Add 'cov' to display samples
  res[["cov"]] <- full_cov

  ## Display the graph of the prediction if expected
  if (plot) {
    ## Check whether training data are available
    if (trained_model %>% is.null()) {
      data_train <- NULL
    } else {
      data_train <- trained_model$ini_args$data
    }

    ## Display samples only in 1D and Credible Interval otherwise
    if(ncol(mixture_pred) == 4){
      ## Plot 1D prediction
      plot_magmaclust(
        res,
        data = data,
        data_train = data_train,
        prior_mean = hyperpost$mean,
        samples = TRUE
      ) %>%
        print()
    } else {
      ## Plot 2D prediction
      plot_magmaclust(
        res,
        data = data,
        samples = FALSE
      ) %>%
        print()
    }

  }

  ## Check whether posterior covariance should be returned
  if (!get_full_cov) {
    res[["cov"]] <- NULL
  }

  return(res)
}

#' Indicates the most probable cluster
#'
#' @param mixture A tibble or data frame containing mixture probabilities.
#'
#' @return A tibble, retaining only the most probable cluster. The column
#'    \code{Cluster} indicates the the cluster's name whereas \code{Proba}
#'    refers to its associated probability. If \code{Task_ID} is initially
#'    a column of \code{mixture} (optional), the function returns the most
#'    probable cluster for all the different \code{Task_ID} values.
#'
#' @export
#'
#' @examples
#' TRUE
proba_max_cluster <- function(mixture) {
  if ("Task_ID" %in% names(mixture)) {
    mixture %>%
      tidyr::pivot_longer(-Task_ID) %>%
      dplyr::group_by(Task_ID) %>%
      dplyr::filter(value == max(value)) %>%
      dplyr::rename("Cluster" = name, "Proba" = value)
  } else {
    mixture %>%
      tidyr::pivot_longer(tidyselect::everything()) %>%
      dplyr::filter(value == max(value)) %>%
      dplyr::rename("Cluster" = name, "Proba" = value)
  }
}

#' Allocate training data into the most probable cluster
#'
#' @param trained_model A list, containing  the information coming from a
#'    MagmaClust model, previously trained using the
#'    \code{\link{train_magmaclust}} function.
#'
#' @return The original dataset used to train the MagmaClust model, with
#'    additional 'Cluster' and associated 'Proba' columns, indicating the most
#'    probable cluster for each task at the end of the training
#'    procedure.
#'
#' @export
#'
#' @examples
#' TRUE
data_allocate_cluster <- function(trained_model) {
  ## Check format and extract useful information
  if (!is.list(trained_model)) {
    stop("Please provide a trained MagmaClust model.")
  } else {
    db <- trained_model$ini_args$data
    mixture <- trained_model$hyperpost$mixture
  }

  ## Remove old 'Cluster' column if necessary
  if ("Cluster" %in% names(db)) {
    db <- db %>% dplyr::select(-.data$Cluster)
  }

  max_clust <- proba_max_cluster(mixture)
  db %>%
    dplyr::left_join(max_clust, by = "ID") %>%
    return()
}

