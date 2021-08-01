#' Gaussian Process prediction
#'
#' Compute the posterior distribution of a simple GP, using the formalism of
#' \code{magma}. By providing observed data, the prior mean and covariance
#' matrix (by defining a kernel and its associated hyper-parameters), the mean
#' and covariance parameters of the posterior distribution are computed on the
#' grid of inputs that has been specified. This predictive distribution can be
#' evaluated on any arbitrary inputs since a GP is an infinite-dimensional
#' object.
#'
#' @param data  A tibble or data frame. Columns required: \code{Input},
#'    \code{Output}. Additional columns for covariates can be specified.
#'    The \code{Input} column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    \code{Output} column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference \code{Input}.
#' @param hp A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern}. The columns/elements should be named
#'    according to the hyper-parameters that are used in \code{kern}. The
#'    function \code{\link{train_gp}} can be used to learn maximum-likelihood
#'    estimators of the hyper-parameters,
#' @param kern A kernel function, defining the covariance structure of the GP.
#' @param mean Mean parameter of the GP. This argument can be specified under
#'    various formats, such as:
#'    - NULL (default). The mean would be set to 0 everywhere.
#'    - A number. The mean would be a constant function.
#'    - A function. This function is defined as the mean.
#'    - A tibble or data frame. Columns required: Input, Output. The Input
#'     values should include at least the same values as in the \code{data}
#'     argument.
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
#'    in which the tibble described above is defined as 'pred_gp' and the full
#'    posterior covariance matrix is defined as 'cov'.
#' @export
#'
#' @examples
#' db <- simu_db(M = 1, N = 10)
#' grid_inputs <- tibble::tibble(Input = 1:10, Covariate = 2:11)
#' hp <- tibble::tibble("variance" = 2, "lengthscale" = 1)
#'
#' pred_gp(db, hp, grid_inputs = grid_inputs)
pred_gp <- function(data,
                    hp,
                    kern = "SE",
                    mean = NULL,
                    grid_inputs = NULL,
                    get_full_cov = FALSE,
                    plot = T,
                    pen_diag = 0.01) {
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
  }

  ## Define the target inputs to predict
  if (grid_inputs %>% is.null()) {
    ## Test whether 'data' only provide the Input column and no covariates
    if (inputs_obs %>% names() %>% length() == 1) {
      input_pred <- seq(min(data$Input), max(data$Input), length.out = 500)
      inputs_pred <- tibble::tibble("Input" = input_pred)
    }
    else {
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
    }
    else {
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
  }
  else if (mean %>% is.vector()) {
    if (length(mean) == 1) {
      mean_obs <- rep(mean, length(input_obs))
      mean_pred <- rep(mean, length(input_pred))
      cat(
        "The mean function has been set to be", mean, "everywhere.\n \n"
      )
    }
    else {
      stop(
        "Incorrect format for the 'mean' argument. Please read ",
        "?pred_gp() for details."
      )
    }
  }
  else if (mean %>% is.function()) {
    mean_obs <- mean(input_obs)
    mean_pred <- mean(input_pred)
  }
  else if (mean %>% is.data.frame()) {
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
  }
  else {
    stop(
      "Incorrect format for the 'mean' argument. Please read ",
      "?pred_gp() for details."
    )
  }

  ## Compute the required sub-matrix for prediction
  inv_obs <- kern_to_inv(inputs_obs, kern, hp, pen_diag)
  cov_pred <- kern_to_cov(inputs_pred, kern, hp)

  ## Check the nature of the 'kern' argument
  if (is.character(kern)) {
    if (kern == "SE") {
      kernel <- se_kernel
    }
    else if (kern == "PERIO") {
      kernel <- perio_kernel
    }
    else if (kern == "RQ") {
      kernel <- rq_kernel
    }
  }
  else if (is.function(kern)) {
    kernel <- kern
  }
  ## Transform the batches of input into lists
  l_inputs_obs <- split(
    t(inputs_obs),
    rep(1:nrow(inputs_obs),
      each = ncol(inputs_obs)
    )
  )
  l_inputs_pred <- split(
    t(inputs_pred),
    rep(1:nrow(inputs_pred),
      each = ncol(inputs_pred)
    )
  )
  ## Compute the covariance matrix between observed and predicted inputs
  cov_crossed <- outer(
    l_inputs_obs, l_inputs_pred,
    Vectorize(function(x, y) kernel(x, y, hp))
  ) %>%
    `rownames<-`(as.character(input_obs)) %>%
    `colnames<-`(as.character(input_pred))

  ## Compute the posterior mean
  pred_mean <- (mean_pred +
    t(cov_crossed) %*% inv_obs %*% (data_obs - mean_obs)) %>%
    as.vector()

  ## Compute the posterior covariance matrix
  pred_cov <- cov_pred - t(cov_crossed) %*% inv_obs %*% cov_crossed

  ## Create a tibble of values and associated uncertainty from a GP prediction
  pred_gp <- tibble::tibble(
    "Mean" = pred_mean,
    "Var" = pred_cov %>% diag()
  ) %>%
    dplyr::mutate(inputs_pred)

  ## Display the graph of the prediction if expected
  if (plot) {
    plot_gp(pred_gp, data) %>% print()
  }

  ## Add the posterior covariance matrix in the results if expected
  if (get_full_cov) {
    list("pred_gp" = pred_gp, "cov" = pred_cov) %>% return()
  }
  else {
    pred_gp %>% return()
  }
}

#' Magma prediction
#'
#' @param db tibble of data columns required ('input', 'Output')
#' @param grid_inputs grid_inputs on which we want a prediction
#' @param mean_mu mean value of mean GP at grid_inputs (obs + pred) (matrix dim: grid_inputs x 1, with Input rownames)
#' @param cov_mu covariance of mean GP at grid_inputs (obs + pred) (square matrix, with Input row/colnames)
#' @param hp list of hyperparameters for the kernel of the GP
#' @param kern kernel associated to the covariance function of the GP
#'
#' @return pamameters of the gaussian density predicted at grid_inputs
#' @export
#'
#' @examples
#' db <- simu_db(M = 1, N = 10)
#' grid_inputs <- tibble::tibble(Input = 1:10, Covariate = 2:11)
#' hp <- tibble::tibble("variance" = 2, "lengthscale" = 1)
#' all_input = union(db$Input, grid_inputs$Input) %>% sort()
#' hyperpost <- list(
#'   "mean" = tibble::tibble(Input = all_input, Output = 0),
#'   "cov" = kern_to_cov(1:10, "SE", hp)
#' )
#'
#' pred_magma(db, hp, grid_inputs = grid_inputs, trained_model = bla)
pred_magma <- function(data,
                       hp,
                       kern = "SE",
                       trained_model = NULL,
                       grid_inputs = NULL,
                       hyperpost = NULL,
                       plot = TRUE,
                       get_hyperpost = FALSE,
                       get_full_cov = FALSE,
                       pen_diag = 0.01) {
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
  }
browser()
  ## Define the target inputs to predict
  if (grid_inputs %>% is.null()) {
    ## Test whether 'data' only provide the Input column and no covariates
    if (inputs_obs %>% names() %>% length() == 1) {
      input_pred <- seq(min(data$Input), max(data$Input), length.out = 500)
      inputs_pred <- tibble::tibble("Input" = input_pred)
    }
    else {
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
    }
    else {
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
  all_input <- union(input_obs, input_pred)
  ## Check whether the hyper-posterior is provided and recompute if necessary
  if (is.null(hyperpost)) {
    if (trained_model %>% is.null()) {
      stop(
        "If the 'hyperpost' argument is NULL, the 'trained_model' ",
        "should be provided, in order to extract or recompute mean process' ",
        "hyper-posterior distribution evaluated on the correct inputs."
      )
    } else if (!all(all_input %in% trained_model$pred_post$Input)) {
      hyperpost <- e_step(
        db = trained_model$fct_args$data,
        m_0 = union(trained_model$fct_args$prior_mean, ,
        kern_0 = trained_model$fct_args$kern_0,
        kern_i = trained_model$fct_args$kern_i,
        hp_0 = trained_model$hp_0,
        hp_i = trained_model$hp_i,
        pen_diag = pen_diag,
        grid_inputs = all_input
      )
    }
  } else if (is.vector(hyperpost$mean$Input) &
    is.matrix(hyperpost$cov)) {
    ## Check whether the inputs in 'hyperpost' are correct
    if (!all(all_input %in% hyperpost$mean$Input) |
      !all(as.character(all_input) %in% colnames(hyperpost$cov))) {
      if (trained_model %>% is.null()) {
        stop(
          "The hyper-posterior distribution of the mean process provided in ",
          "'hyperpost' argument isn't evaluated on the expected inputs. The ",
          "'trained_model' argument is needed for an adequate re-conmputing."
        )
      }
      else {
        cat(
          "The hyper-posterior distribution of the mean process provided in the",
          "'hyperpost' argument isn't evaluated on the expected inputs. \n \n",
          "Start evaluating the hyper-posterior on the correct inputs..."
        )

        hyperpost <- e_step(
          db = trained_model$fct_args$data,
          m_0 = trained_model$fct_args$prior_mean,
          kern_0 = trained_model$fct_args$kern_0,
          kern_i = trained_model$fct_args$kern_i,
          hp_0 = trained_model$hp_0,
          hp_i = trained_model$hp_i,
          pen_diag = pen_diag,
          grid_inputs = all_input
        )

        cat("Done!\n \n")
      }
    }
  } else {
    stop(
      "The format of the 'hyperpost' argument is not as exepected. Please ",
      "read ?pred_magma() for details."
    )
  }
  ## Extract the mean parameter from the hyper-posterior
  mean_obs <- hyperpost$mean %>%
    dplyr::filter(.data$Input %in% input_obs) %>%
    dplyr::arrange(.data$Input) %>%
    dplyr::pull(.data$Output)

  mean_pred <- hyperpost$mean %>%
    dplyr::filter(.data$Input %in% input_pred) %>%
    dplyr::arrange(.data$Input) %>%
    dplyr::pull(.data$Output)

  ## Extract the covariance sub-matrices from the hyper-posterior
  post_cov_obs <- hyperpost$cov[as.character(input_obs), as.character(input_obs)]
  post_cov_pred <- hyperpost$cov[
    as.character(input_pred),
    as.character(input_pred)
  ]
  post_cov_crossed <- hyperpost$cov[
    as.character(input_obs),
    as.character(input_pred)
  ]

  ## Sum the covariance matrices on oberved inputs and compute the inverse
  cov_obs <- kern_to_cov(inputs_obs, kern, hp) + post_cov_obs
  diag <- diag(x = pen_diag, ncol = ncol(cov_obs), nrow = nrow(cov_obs))

  inv_obs <- tryCatch((cov_obs + diag) %>% chol() %>% chol2inv(),
    error = function(e) {
      MASS::ginv(cov_obs + diag)
    }
  ) %>%
    `rownames<-`(as.character(input_obs)) %>%
    `colnames<-`(as.character(input_obs))

  ## Sum the covariance matrices on prediction inputs
  cov_pred <- kern_to_cov(inputs_pred, kern, hp) + post_cov_pred

  ## Check the nature of the 'kern' argument
  if (is.character(kern)) {
    if (kern == "SE") {
      kernel <- se_kernel
    }
    else if (kern == "PERIO") {
      kernel <- perio_kernel
    }
    else if (kern == "RQ") {
      kernel <- rq_kernel
    }
  }
  else if (is.function(kern)) {
    kernel <- kern
  }
  ## Transform the batches of input into lists
  l_inputs_obs <- split(
    t(inputs_obs),
    rep(1:nrow(inputs_obs),
      each = ncol(inputs_obs)
    )
  )
  l_inputs_pred <- split(
    t(inputs_pred),
    rep(1:nrow(inputs_pred),
      each = ncol(inputs_pred)
    )
  )
  ## Compute the covariance matrix between observed and predicted inputs
  cov_crossed <- outer(
    l_inputs_obs, l_inputs_pred,
    Vectorize(function(x, y) kernel(x, y, hp))
  ) %>%
    `+`(post_cov_crossed) %>%
    `rownames<-`(as.character(input_obs)) %>%
    `colnames<-`(as.character(input_pred))

  ## Compute the posterior mean of a GP
  pred_mean <- (mean_pred +
    t(cov_crossed) %*% inv_obs %*% (data_obs - mean_obs)) %>%
    as.vector()

  ## Compute the posterior covariance matrix of a GP
  pred_cov <- cov_pred - t(cov_crossed) %*% inv_obs %*% cov_crossed

  ## Create a tibble of values and associated uncertainty from a GP prediction
  pred_gp <- tibble::tibble(
    "Mean" = pred_mean,
    "Var" = pred_cov %>% diag()
  ) %>%
    dplyr::mutate(inputs_pred)

  ## Display the graph of the prediction if expected
  if (plot) {
    plot_gp(pred_gp, data) %>% print()
  }

  res <- pred_gp
  ## Check whether posterior covariance or hyper-posterior should be returned
  if (get_full_cov | get_hyperpost) {
    res <- list("pred_gp" = pred_gp)
    if (get_full_cov) {
      res[["cov"]] <- pred_cov
    }
    if (get_hyperpost) {
      res[["hyperpost"]] <- hyperpost
    }
  }

  return(res)
}

#' Prediction Gaussian Process Animate
#'
#' Function used to generate data compatible with a GIF ploting of the results
#'
#' @param db Your database, the same as for a classic GP prediction
#' @param timestamps dl
#' @param mean_mu mean value of mean GP at timestamps (obs + pred)
#' @param cov_mu covariance of mean GP at timestamps (obs + pred)
#' @param kern kernel of your choice use with your hyperparameters
#' @param hp list of hyperparameters of your model
#'
#' @return tibble of classic GP predictions but with an inscreasing number of data points considered as 'observed'
#' @export
#'
#' @examples
#' TRUE
pred_magma_animate <- function(db, timestamps = NULL, mean_mu = 0, cov_mu = NULL,
                               kern, hp) {
  db %>% dplyr::arrange(db$input)
  all_pred <- tibble::tibble()

  if (is.null(grid_inputs)) {
    grid_inputs <- seq(min(db$input), max(db$input), length.out = 500)
  }

  for (j in 1:nrow(db))
  {
    pred_j <- pred_magma(db[1:j, ], grid_inputs, mean_mu, cov_mu, kern, hp) %>% dplyr::mutate(Nb_data = j)
    all_pred <- all_pred %>% rbind(pred_j)
  }
  return(all_pred)
}
