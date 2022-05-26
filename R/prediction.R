#' Gaussian Process prediction
#'
#' Compute the posterior distribution of a simple GP, using the formalism of
#' Magma. By providing observed data, the prior mean and covariance
#' matrix (by defining a kernel and its associated hyper-parameters), the mean
#' and covariance parameters of the posterior distribution are computed on the
#' grid of inputs that has been specified. This predictive distribution can be
#' evaluated on any arbitrary inputs since a GP is an infinite-dimensional
#' object.
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
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#' @param grid_inputs The grid of inputs (reference Input and covariates) values
#'    on which the GP should be evaluated. Ideally, this argument should be a
#'    tibble or a data frame, providing the same columns as \code{data}, except
#'    'Output'. Nonetheless, in cases where \code{data} provides only one
#'    'Input' column, the \code{grid_inputs} argument can be NULL (default) or a
#'    vector. This vector would be used as reference input for prediction and if
#'    NULL, a vector of length 500 is defined, ranging between the min and max
#'    Input values of \code{data}.
#' @param get_full_cov A logical value, indicating whether the full posterior
#'    covariance matrix should be returned.
#' @param plot A logical value, indicating whether a plot of the results is
#'    automatically displayed.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#'
#' @return A tibble, representing the GP predictions as two column 'Mean' and
#'   'Var', evaluated on the \code{grid_inputs}. The column 'Input' and
#'   additional covariates columns are associated to each predicted values.
#'    If the \code{get_full_cov} argument is TRUE, the function returns a list,
#'    in which the tibble described above is defined as 'pred' and the full
#'    posterior covariance matrix is defined as 'cov'.
#' @export
#'
#' @examples
#' db <- simu_db(M = 1, N = 10)
#' grid_inputs <- tibble::tibble(Input = 0:10, Covariate = 1:11)
#' hp <- tibble::tibble("variance" = 2, "lengthscale" = 1)
#'
#' pred_gp(db, grid_inputs = grid_inputs)
pred_gp <- function(data,
                    mean = NULL,
                    hp = NULL,
                    kern = "SE",
                    grid_inputs = NULL,
                    get_full_cov = FALSE,
                    plot = TRUE,
                    pen_diag = 1e-8) {
  ## Extract the observed Output (data points)
  data_obs <- data %>%
    dplyr::arrange(.data$Input) %>%
    dplyr::pull(.data$Output)

  ## Extract the observed (reference) Input
  input_obs <- data %>%
    dplyr::arrange(.data$Input) %>%
    dplyr::pull(.data$Input)
  ## Extract the observed inputs (reference Input + covariates)
  inputs_obs <- data %>%
    dplyr::arrange(.data$Input) %>%
    dplyr::select(-.data$Output)

  ## Remove the 'ID' column if present
  if ("ID" %in% names(data)) {
    inputs_obs <- inputs_obs %>% dplyr::select(-.data$ID)
    if (dplyr::n_distinct(data$ID) > 1) {
      stop(
        "Problem in the 'ID' column: different values are not allowed. ",
        "The prediction can only be performed for one individual/task."
      )
    }
  }

  ## Define the target inputs to predict
  if (grid_inputs %>% is.null()) {
    ## Test whether 'data' only provide the Input column and no covariates
    if (inputs_obs %>% names() %>% length() == 1) {
      input_pred <- seq(min(data$Input), max(data$Input), length.out = 500)
      inputs_pred <- tibble::tibble("Input" = input_pred)
    } else if (inputs_obs %>% names() %>% length() == 2) {
      ## Define a default grid for 'Input'
      input_pred <- rep(
        seq(min(data$Input), max(data$Input), length.out = 20),
        each = 20
      )
      inputs_pred <- tibble::tibble("Input" = input_pred)
      ## Add a grid for the covariate
      name_cova <- inputs_obs %>%
        dplyr::select(-.data$Input) %>%
        names()
      cova <- inputs_obs[name_cova]
      inputs_pred[name_cova] <- rep(
        seq(min(cova), max(cova), length.out = 20),
        times = 20
      )
    } else {
      stop(
        "The 'grid_inputs' argument should be a either a numerical vector ",
        "or a data frame depending on the context. Please read ?pred_gp()."
      )
    }
  } else if (grid_inputs %>% is.vector()) {
    ## Test whether 'data' only provide the Input column and no covariates
    if (inputs_obs %>% names() %>% length() == 1) {
      input_pred <- grid_inputs %>% sort()
      inputs_pred <- tibble::tibble("Input" = input_pred)
    } else {
      stop(
        "The 'grid_inputs' argument should be a either a numerical vector ",
        "or a data frame depending on the context. Please read ?pred_gp()."
      )
    }
  } else if (grid_inputs %>% is.data.frame()) {
    ## Test whether 'data' has the same columns as grid_inputs
    if (all(names(inputs_obs) %in% names(grid_inputs))) {
      input_pred <- grid_inputs %>%
        dplyr::arrange(.data$Input) %>%
        dplyr::pull(.data$Input)

      inputs_pred <- grid_inputs %>%
        dplyr::arrange(.data$Input) %>%
        dplyr::select(names(inputs_obs))
    } else {
      stop(
        "The 'grid_inputs' argument should provide a column 'Input', and ",
        "the same additional covariate columns as contained in 'data'."
      )
    }
  } else {
    stop(
      "The 'grid_inputs' argument should be a either a numerical vector ",
      "or a data frame depending on the context. Please read ?pred_gp()."
    )
  }

  ## Define of extract the mean values under an adequate format
  if (mean %>% is.null()) {
    mean_obs <- rep(0, length(input_obs))
    mean_pred <- rep(0, length(input_pred))
    cat(
      "The 'mean' argument has not been specified. The ",
      "mean function is thus set to be 0 everywhere.\n \n"
    )
  } else if (mean %>% is.vector()) {
    if (length(mean) == 1) {
      mean_obs <- rep(mean, length(input_obs))
      mean_pred <- rep(mean, length(input_pred))
      cat(
        "The mean function has been set to be", mean, "everywhere.\n \n"
      )
    } else {
      stop(
        "Incorrect format for the 'mean' argument. Please read ",
        "?pred_gp() for details."
      )
    }
  } else if (mean %>% is.function()) {
    mean_obs <- mean(input_obs)
    mean_pred <- mean(input_pred)
  } else if (mean %>% is.data.frame()) {
    if (all(c("Output", "Input") %in% names(mean))) {
      mean_obs <- mean %>%
        dplyr::filter(.data$Input %in% input_obs) %>%
        dplyr::arrange(.data$Input) %>%
        dplyr::pull(.data$Output)

      mean_pred <- mean %>%
        dplyr::filter(.data$Input %in% input_pred) %>%
        dplyr::arrange(.data$Input) %>%
        dplyr::pull(.data$Output)

      if ((length(mean_obs) != length(input_obs)) |
        (length(mean_pred) != length(input_pred))) {
        stop(
          "Problem in the length of the mean parameter. The ",
          "'mean' argument should provide an Output value for each Input ",
          "value appearing in the data and for each 'grid_inputs' values."
        )
      }
    } else {
      stop(
        "If the 'mean' argument is provided as a data frame, it ",
        "should contain the mandatory column names: 'Output', 'Input'"
      )
    }
  } else {
    stop(
      "Incorrect format for the 'mean' argument. Please read ",
      "?pred_gp() for details."
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
    hp_rm_noi <- hp %>% dplyr::select(-.data$noise)
    noise <- exp(hp[["noise"]])
  } else {
    hp_rm_noi <- hp
    noise <- 0
  }

  ## Compute the required sub-matrix for prediction
  inv_obs <- kern_to_inv(inputs_obs, kern, hp, pen_diag)
  cov_pred <- kern_to_cov(inputs_pred, kern, hp_rm_noi)
  cov_crossed <- kern_to_cov(inputs_obs, kern, hp_rm_noi, input_2 = inputs_pred)

  ## Compute the posterior mean
  pred_mean <- (mean_pred +
    t(cov_crossed) %*% inv_obs %*% (data_obs - mean_obs)) %>%
    as.vector()

  ## Compute the posterior covariance matrix
  pred_cov <- cov_pred - t(cov_crossed) %*% inv_obs %*% cov_crossed

  ## Create a tibble of values and associated uncertainty from a GP prediction
  pred_gp <- tibble::tibble(
    "Mean" = pred_mean,
    "Var" = diag(pred_cov) + noise
  ) %>%
    dplyr::mutate(inputs_pred)

  ## Display the graph of the prediction if expected
  if (plot) {
    plot_gp(pred_gp, data = data) %>% print()
  }

  ## Add the posterior covariance matrix in the results if expected
  if (get_full_cov) {
    list("pred" = pred_gp, "cov" = pred_cov) %>% return()
  } else {
    pred_gp %>% return()
  }
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
#' @param data A tibble or data frame. Required columns: 'Input',
#'    'Output'. Additional columns for covariates can be specified.
#'    The 'Input' column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    'Output' column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference 'Input'.
#' @param hp_0 A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_0}.
#' @param hp_i A tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}.
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
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#' @param kern_i A kernel function, associated with the individual GPs. ("SE",
#'    "PERIO" and "RQ" are aso available here)
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
#' @param grid_inputs A vector, indicating the grid of additional reference
#'    inputs on which the mean process' hyper-posterior should be evaluated.
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
#' db <- simu_db(N = 10, common_input = TRUE)
#' hp_0 <- hp()
#' hp_i <- hp("SE", list_ID = unique(db$ID))
#' grid_inputs <- seq(0, 10, 0.1)
#' hyperposterior(db, hp_0, hp_i, "SE", "SE", grid_inputs = grid_inputs)
hyperposterior <- function(data,
                           hp_0,
                           hp_i,
                           kern_0,
                           kern_i,
                           prior_mean = NULL,
                           grid_inputs = NULL,
                           pen_diag = 1e-8) {
  if (grid_inputs %>% is.null()) {
    ## Define the union of all reference Inputs in the dataset
    all_input <- unique(data$Input) %>% sort()
    cat(
      "The argument 'grid_inputs' is NULL, the hyper-posterior distribution",
      "will only be evaluated on observed Input from 'data'.\n \n"
    )
  } else {
    ## Define the union among all reference Inputs and a specified grid
    all_input <- unique(data$Input) %>%
      union(grid_inputs) %>%
      sort()
  }

  ## Initialise m_0 according to the value provided by the user
  if (prior_mean %>% is.null()) {
    m_0 <- rep(0, length(all_input))
    cat(
      "The 'prior_mean' argument has not been specified. The hyper_prior mean",
      "function is thus set to be 0 everywhere.\n \n"
    )
  } else if (prior_mean %>% is.vector()) {
    if (length(prior_mean) == length(all_input)) {
      m_0 <- prior_mean
    } else if (length(prior_mean) == 1) {
      m_0 <- rep(prior_mean, length(all_input))
      cat(
        "The provided 'prior_mean' argument is of length 1. Thus, the",
        "hyper-prior mean function has set to be constant everywhere.\n \n"
      )
    } else {
      stop(
        "The 'prior_mean' argument is of length ", length(prior_mean),
        ", whereas the grid of training inputs is of length ",
        length(all_input)
      )
    }
  } else if (prior_mean %>% is.function()) {
    m_0 <- prior_mean(all_input)
  } else if (prior_mean %>% is.data.frame()) {
    if (all(c("Output", "Input") %in% names(prior_mean))) {
      m_0 <- prior_mean %>%
        dplyr::filter(.data$Input %in% all_input) %>%
        dplyr::arrange(.data$Input) %>%
        dplyr::pull(.data$Output)

      if (length(m_0) != length(all_input)) {
        stop(
          "Problem in the length of the hyper_prior mean parameter. The ",
          "'pior_mean' argument should provide an Output value for each Input ",
          "value appearing in the training data."
        )
      }
    } else {
      stop(
        "If the 'prior_mean' argument is provided as a data frame, it ",
        "should contain the mandatory column names: 'Output', 'Input'"
      )
    }
  } else {
    stop(
      "Incorrect format for the 'prior_mean' argument. Please read ",
      "?hyperposterior() for details."
    )
  }

  ## Compute all the inverse covariance matrices
  inv_0 <- kern_to_inv(all_input, kern_0, hp_0, pen_diag)
  list_inv_i <- list_kern_to_inv(data, kern_i, hp_i, pen_diag)
  ## Create a named list of Output values for all individuals
  list_output_i <- base::split(data$Output, list(data$ID))

  ## Update the posterior inverse covariance ##
  post_inv <- inv_0
  for (inv_i in list_inv_i)
  {
    ## Collect the input's common indices between mean and individual processes
    co_input <- intersect(row.names(inv_i), row.names(post_inv))
    ## Sum the common inverse covariance's terms
    post_inv[co_input, co_input] <- post_inv[co_input, co_input] +
      inv_i[co_input, co_input]
  }
  ##############################################

  ## Update the posterior mean ##
  weighted_0 <- inv_0 %*% m_0

  for (i in names(list_inv_i))
  {
    ## Compute the weighted mean for the i-th individual
    weighted_i <- list_inv_i[[i]] %*% list_output_i[[i]]
    ## Collect the input's common indices between mean and individual processes
    co_input <- intersect(row.names(weighted_i), row.names(weighted_0))
    ## Sum the common weighted mean terms
    weighted_0[co_input, ] <- weighted_0[co_input, ] +
      weighted_i[co_input, ]
  }

  ## Fast or slow matrix inversion if nearly singular
  post_cov <- tryCatch(post_inv %>% chol() %>% chol2inv(), error = function(e) {
    MASS::ginv(post_inv)
  }) %>%
    `rownames<-`(all_input) %>%
    `colnames<-`(all_input)
  ## Compute the updated mean parameter
  post_mean <- post_cov %*% weighted_0 %>% as.vector()
  ##############################################

  ## Format the mean parameter of the hyper-posterior distribution
  mean <- tibble::tibble(
    "Input" = all_input,
    "Output" = post_mean
  )

  ## Format the GP prediction of the hyper-posterior mean (for direct plot)
  pred <- tibble::tibble(
    "Input" = all_input,
    "Mean" = post_mean,
    "Var" = post_cov %>% diag() %>% as.vector()
  )
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
#' new individual/task, its trained hyper-parameters and a previously trained
#' Magma model, the predictive distribution is evaluated on any arbitrary inputs
#' that are specified through the 'grid_inputs' argument.
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
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#' @param grid_inputs The grid of inputs (reference Input and covariates) values
#'    on which the GP should be evaluated. Ideally, this argument should be a
#'    tibble or a data frame, providing the same columns as \code{data}, except
#'    'Output'. Nonetheless, in cases where \code{data} provides only one
#'    'Input' column, the \code{grid_inputs} argument can be NULL (default) or a
#'    vector. This vector would be used as reference input for prediction and if
#'    NULL, a vector of length 500 is defined, ranging between the min and max
#'    Input values of \code{data}.
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
#' @return A tibble, representing Magma predictions as two column 'Mean' and
#'   'Var', evaluated on the \code{grid_inputs}. The column 'Input' and
#'   additional covariates columns are associated to each predicted values.
#'    If the \code{get_full_cov} or \code{get_hyperpost} arguments are TRUE,
#'    the function returns a list, in which the tibble described above is
#'    defined as 'pred_gp' and the full posterior covariance matrix is
#'    defined as 'cov', and the hyper-posterior distribution of the mean process
#'    is defined as 'hyperpost'.
#' @export
#'
#' @examples
#' db <- simu_db(M = 1, N = 10, covariate = FALSE)
#' grid_inputs <- seq(0, 10, 0.1)
#' all_input <- union(db$Input, grid_inputs) %>% sort()
#' hyperpost <- list(
#'   "mean" = tibble::tibble(Input = all_input, Output = 0),
#'   "cov" = kern_to_cov(all_input, "SE", hp("SE"))
#' )
#'
#' pred_magma(db, grid_inputs = grid_inputs, hyperpost = hyperpost)
pred_magma <- function(data,
                       trained_model = NULL,
                       hp = NULL,
                       kern = "SE",
                       grid_inputs = NULL,
                       hyperpost = NULL,
                       get_hyperpost = FALSE,
                       get_full_cov = FALSE,
                       plot = TRUE,
                       pen_diag = 1e-8) {
  ## Extract the observed Output (data points)
  data_obs <- data %>%
    dplyr::arrange(.data$Input) %>%
    dplyr::pull(.data$Output)

  ## Extract the observed (reference) Input
  input_obs <- data %>%
    dplyr::arrange(.data$Input) %>%
    dplyr::pull(.data$Input)

  ## Extract the observed inputs (reference Input + covariates)
  inputs_obs <- data %>%
    dplyr::arrange(.data$Input) %>%
    dplyr::select(-.data$Output)

  ## Remove the 'ID' column if present
  if ("ID" %in% names(data)) {
    inputs_obs <- inputs_obs %>% dplyr::select(-.data$ID)
    if (dplyr::n_distinct(data$ID) > 1) {
      stop(
        "Problem in the 'ID' column: different values are not allowed. ",
        "The prediction can only be performed for one individual/task."
      )
    }
  }

  ## Define the target inputs to predict
  if (grid_inputs %>% is.null()) {
    ## Test whether 'data' only provide the Input column and no covariates
    if (inputs_obs %>% names() %>% length() == 1) {
      input_pred <- seq(min(data$Input), max(data$Input), length.out = 500)
      inputs_pred <- tibble::tibble("Input" = input_pred)
    } else if (inputs_obs %>% names() %>% length() == 2) {
      ## Define a default grid for 'Input'
      input_pred <- rep(
        seq(min(data$Input), max(data$Input), length.out = 20),
        each = 20
      )
      inputs_pred <- tibble::tibble("Input" = input_pred)
      ## Add a grid for the covariate
      name_cova <- inputs_obs %>%
        dplyr::select(-.data$Input) %>%
        names()
      cova <- inputs_obs[name_cova]
      inputs_pred[name_cova] <- rep(
        seq(min(cova), max(cova), length.out = 20),
        times = 20
      )
    } else {
      stop(
        "The 'grid_inputs' argument should be a either a numerical vector ",
        "or a data frame depending on the context. Please read ?pred_magma()."
      )
    }
  } else if (grid_inputs %>% is.vector()) {
    ## Test whether 'data' only provide the Input column and no covariates
    if (inputs_obs %>% names() %>% length() == 1) {
      input_pred <- grid_inputs %>% sort()
      inputs_pred <- tibble::tibble("Input" = input_pred)
    } else {
      stop(
        "The 'grid_inputs' argument should be a either a numerical vector ",
        "or a data frame depending on the context. Please read ?pred_magma()."
      )
    }
  } else if (grid_inputs %>% is.data.frame()) {
    ## Test whether 'data' has the same columns as grid_inputs
    if (names(inputs_obs) %>% setequal(names(grid_inputs))) {
      input_pred <- grid_inputs %>%
        dplyr::arrange(.data$Input) %>%
        dplyr::pull(.data$Input)

      inputs_pred <- grid_inputs %>%
        dplyr::arrange(.data$Input) %>%
        dplyr::select(names(inputs_obs))
    } else {
      stop(
        "The 'grid_inputs' argument should provide a column 'Input', and ",
        "the same additional covariate columns as contained in 'data'."
      )
    }
  } else {
    stop(
      "The 'grid_inputs' argument should be a either a numerical vector ",
      "or a data frame depending on the context. Please read ?pred_gp()."
    )
  }

  ## Define the union of all distinct reference Input
  all_input <- union(input_obs, input_pred) %>% sort()

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
    if (!all(all_input %in% hyperpost$mean$Input)) {
      cat(
        "The hyper-posterior distribution of the mean process provided in",
        "'hyperpost' argument isn't evaluated on the expected inputs.\n \n",
        "Start evaluating the hyper-posterior on the correct inputs...\n \n"
      )
      hyperpost <- hyperposterior(
        data = trained_model$ini_args$data,
        kern_0 = trained_model$ini_args$kern_0,
        kern_i = trained_model$ini_args$kern_i,
        hp_0 = trained_model$hp_0,
        hp_i = trained_model$hp_i,
        prior_mean = trained_model$ini_args$prior_mean,
        grid_inputs = all_input,
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
    dplyr::right_join(inputs_obs, by = "Input") %>%
    dplyr::arrange(.data$Input) %>%
    dplyr::pull(.data$Output)

  mean_pred <- hyperpost$mean %>%
    dplyr::right_join(inputs_pred, by = "Input") %>%
    dplyr::arrange(.data$Input) %>%
    dplyr::pull(.data$Output)

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
      ## Check whether hyper-parameters are common if we have 'trained_model'
      if (
        tryCatch(trained_model$ini_args$common_hp, error = function(e) FALSE)
        ) {
        ## Extract the hyper-parameters common to all 'i'
        hp <- trained_model$hp_i %>%
          dplyr::slice(1) %>%
          dplyr::select(-.data$ID)
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

  ## Sum the covariance matrices on observed inputs and compute the inverse
  cov_obs <- kern_to_cov(inputs_obs, kern, hp) + post_cov_obs
  diag <- diag(x = pen_diag, ncol = ncol(cov_obs), nrow = nrow(cov_obs))

  inv_obs <- tryCatch((cov_obs + diag) %>% chol() %>% chol2inv(),
    error = function(e) {
      MASS::ginv(cov_obs + diag)
    }
  ) %>%
    `rownames<-`(as.character(input_obs)) %>%
    `colnames<-`(as.character(input_obs))

  ## Remove the noise of the hp for evaluating some of the sub-matrix
  if ("noise" %in% names(hp)) {
    hp_rm_noi <- hp %>% dplyr::select(-.data$noise)
    noise <- exp(hp[["noise"]])
  } else {
    hp_rm_noi <- hp
    noise <- 0
  }

  ## Compute the required sub-matrix for prediction
  cov_pred <- kern_to_cov(inputs_pred, kern, hp_rm_noi) + post_cov_pred
  cov_crossed <- kern_to_cov(inputs_obs, kern, hp_rm_noi,
                             input_2 = inputs_pred) +
    post_cov_crossed

  ## Compute the posterior mean of a GP
  pred_mean <- (mean_pred +
    t(cov_crossed) %*% inv_obs %*% (data_obs - mean_obs)) %>%
    as.vector()

  ## Compute the posterior covariance matrix of a GP
  pred_cov <- cov_pred - t(cov_crossed) %*% inv_obs %*% cov_crossed

  ## Create a tibble of values and associated uncertainty from a GP prediction
  pred_gp <- tibble::tibble(
    "Mean" = pred_mean,
    "Var" = diag(pred_cov) + noise
  ) %>%
    dplyr::mutate(inputs_pred)

  ## Display the graph of the prediction if expected
  if (plot) {
    plot_gp(pred_gp, data = data, prior_mean = hyperpost$mean) %>%
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

  return(res)
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
#' @param mean Mean parameter of the GP. This argument can be specified under
#'    various formats, such as:
#'    - NULL (default). The mean would be set to 0 everywhere.
#'    - A number. The mean would be a constant function.
#'    - A function. This function is defined as the mean.
#'    - A tibble or data frame. Required columns: Input, Output. The Input
#'     values should include at least the same values as in the \code{data}
#'     argument.
#' @param trained_model A list, containing  the information coming from a
#'    Magma model, previously trained using the \code{\link{train_magma}}
#'    function.
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
#' @param grid_inputs The grid of inputs (reference Input and covariates) values
#'    on which the GP should be evaluated. Ideally, this argument should be a
#'    tibble or a data frame, providing the same columns as \code{data}, except
#'    'Output'. Nonetheless, in cases where \code{data} provides only one
#'    'Input' column, the \code{grid_inputs} argument can be NULL (default) or a
#'    vector. This vector would be used as reference input for prediction and if
#'    NULL, a vector of length 500 is defined, ranging between the min and max
#'    Input values of \code{data}.
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
#' \donttest{
#' db <- simu_db(M = 1)
#' grid_inputs <- tibble::tibble(
#'   Input = seq(0, 10, 0.1),
#'   Covariate = seq(-5, 5, 0.1)
#' )
#' pred_gif(db, grid_inputs = grid_inputs)
#' }
pred_gif <- function(data,
                     trained_model = NULL,
                     hyperpost = NULL,
                     mean = NULL,
                     hp = NULL,
                     kern = "SE",
                     grid_inputs = NULL,
                     pen_diag = 1e-8) {
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
