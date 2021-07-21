#' Training Magma with an EM algorithm
#'
#' The hyper-parameters and the hyper-posterior distribution involved in Magma
#' can be learned thanks to an EM algorithm implemented in \code{train_magma}.
#' By providing a dataset, the model hypotheses (prior mean parameter and
#' covariance kernels) and initialisation values for the hyper-parameters, the
#' function compute maximum likelihood estimates of the HPs as well as the
#' mean and covariance parameters of the Gaussian hyper-posterior distribution
#' of the mean process.
#'
#' @param data A tibble or data frame. Columns required: \code{ID}, \code{Input}
#'    , \code{Output}.
#'    Additional columns for covariates can be specified.
#'    The \code{ID} column contains the unique names/codes used to identify each
#'    individual/task (or batch of data).
#'    The \code{Input} column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    \code{Output} column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at each
#'    reference \code{Input}.
#' @param prior_mean Prior mean parameter (m_0) of the mean GP. This argument
#'    can be specified under various formats, such as:
#'    - NULL (default). The prior mean would be set to 0 everywhere.
#'    - A number. The prior mean would be a constant function.
#'    - A vector of the same length as all the distinct Input values in the
#'     \code{data} argument. This vector would be considered as the evaluation
#'     of the prior mean function at the training Inputs.
#'    - A function. This function is defined as the prior mean.
#'    - A tibble or data frame. Columns required: Input, Output. The Input
#'     values should include at least the same values as in the \code{data}
#'     argument.
#'
#' @param ini_hp_0 A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_0}, the mean process' kernel. The
#'    columns/elements should be named according to the hyper-parameters
#'    that are used in \code{kern_0}.
#' @param ini_hp_i A tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}, the individual processes' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each individual/task. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_i}.
#' @param kern_0 A kernel function, associated with the mean GP.
#' @param kern_i A kernel function, associated with the individual GPs.
#' @param common_hp A logical value, indicating whether the set of
#'    hyper-parameters is assumed to be common to all indiviuals.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#' @param grid_inputs A vector, indicating the grid of additional reference
#'    inputs on which the mean process' hyper-posterior should be evaluated.
#' @param n_iter_max A number, indicating the maximum number of iterations of
#'    the EM algorithm to proceed while not reaching convergence.
#' @param cv_threshold A number, indicating the threshold of the likelihood gain
#'    under which the EM algorithm will stop. The convergence condition is
#'    defined as the difference of likelihoods between two consecutive steps,
#'    divided by the absolute value of the last one
#'    ( (LL_n - LL_n-1) / |LL_n| ).
#'
#' @details The user can specify custom kernel functions for the argument
#'    \code{kern_0} and \code{kern_i}. The hyper-parameters used in the kernel
#'    should have explicit names, and be contained within the \code{hp}
#'    argument. \code{hp} should typically be defined as a named vector or a
#'    data frame. Although it is not mandatory for the \code{train_magma}
#'    function to run, gradients can provided within kernel function definition.
#'    See for example \code{\link{se_kernel}} to create a custom kernel
#'    function displaying an adequate format to be used in Magma.
#'
#' @return A list, containing the results of the EM algorithm used for training
#'    in Magma. The elements of the list are:
#'    - hp_0: A tibble containing the trained hyper-parameters for the mean
#'    process' kernel.
#'    - hp_i: A tibble containing all the trained hyper-parameters for the
#'    individual processes' kernels.
#'    - post_mean: A tibble containing the values of hyper-posterior's mean
#'    parameter (\code{Output}) evaluated at each training reference
#'    \code{Input}.
#'    - post_cov: A matrix, covariance parameter of the hyper-posterior
#'    distribution of the mean process.
#'    - pred_post: A tibble, gathering mean and covariance parameters of the
#'    mean process' hyper-posterior distribution under a format that allows
#'    direct visualisation as a GP prediction.
#'    - Converged: A logical value indicated whether the EM algorithm converged
#'    or not.
#'    - Training_time: Total running time of the complete training.
#'
#' @export
#'
#' @examples
#' db = simu_db()
#' train_magma(db)
train_magma <- function(data,
                        prior_mean = NULL,
                        ini_hp_0 = NULL,
                        ini_hp_i = NULL,
                        kern_0 = "SE",
                        kern_i = "SE",
                        common_hp = T,
                        grid_inputs = NULL,
                        pen_diag = 0.01,
                        n_iter_max = 25,
                        cv_threshold = 1e-3) {

  ## Check for the correct format of the training data
  if (data %>% is.data.frame()) {
    if (!all(c("ID", "Output", "Input") %in% names(data))) {
      stop(
        "The 'data' argument should be a tibble or a data frame containing ",
        "at least the mandatory column names: 'ID', 'Output' and 'Input'"
      )
    }
  } else {
    stop(
      "The 'data' argument should be a tibble or a data frame containing ",
      "at least the mandatory column names: 'ID', 'Output' and 'Input'"
    )
  }

  ## Certify that IDs are of type 'character'
  data$ID <- data$ID %>% as.character()
  ## Extract the list of differnt IDs
  list_ID <- data$ID %>% unique()
  ## Extract the union of all reference inputs provided in the training data
  all_inputs <- data$Input %>%
    unique() %>%
    sort()
  ## Define the union of training inputs and the provided grid
  if(!is.null(grid_inputs)){
    if (grid_inputs %>% is.vector()) {
      all_inputs_grid = all_inputs %>%
        union(grid_inputs) %>%
        sort()
    }
    else {
      grid_inputs = NULL
      warning("The argument 'grid_inputs' should be a vector. The ",
            "hyper-posterior will only be evaluated on training data inputs.",
            "\n \n")
    }
  }
  ## Set m_0_grid to NULL in case no correct value is provided
  m_0_grid = NULL

  ## Initialise m_0 according to the value provided by the user
  if (prior_mean %>% is.null()) {
    m_0 <- rep(0, length(all_inputs))
    if(!is.null(grid_inputs)){m_0_grid <- rep(0, length(all_inputs_grid))}
    cat(
      "The 'prior_mean' argument has not been specified. The prior mean",
      "function is thus set to be 0 everywhere.\n \n"
    )
  }
  else if (prior_mean %>% is.vector()) {
    if (length(prior_mean) == length(all_inputs)) {
      m_0 <- prior_mean
    } else if (length(prior_mean) == 1) {
      m_0 <- rep(prior_mean, length(all_inputs))
      if(!is.null(grid_inputs)){
        m_0_grid <- rep(prior_mean, length(all_inputs_grid))
      }
      cat(
        "The provided 'prior_mean' argument is of length 1. Thus, the prior",
        "mean function has set to be constant everywhere.\n \n"
      )
    }
    else {
      stop(
        "The 'prior_mean' argument is of length ", length(prior_mean),
        ", whereas the grid of training inputs is of length ",
        length(all_inputs)
      )
    }
  }
  else if (prior_mean %>% is.function()) {
    m_0 <- prior_mean(all_inputs)
    if(!is.null(grid_inputs)){m_0_grid <- prior_mean(all_inputs_grid)}
  }
  else if (prior_mean %>% is.data.frame()) {
    if (all(c("Output", "Input") %in% names(prior_mean))) {
      m_0 <- prior_mean %>%
        dplyr::filter(.data$Input %in% all_inputs) %>%
        dplyr::pull(.data$Output)
      if(!is.null(grid_inputs)){
      m_0_grid <- prior_mean %>%
          dplyr::filter(.data$Input %in% all_inputs_grid) %>%
          dplyr::pull(.data$Output)
      }
      if (length(m_0) != length(all_inputs)) {
        stop(
          "Problem in the length of the prior mean parameter. The ",
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
  }
  else {
    stop(
      "Incorrect format for the 'prior_mean' argument. Please read ",
      "?train_magma() for details."
    )
  }

  ## Track the total training time
  t_1 <- Sys.time()

  ## Initialise the mean process' hp according to user's values
  if (kern_0 %>% is.function()) {
    if (ini_hp_0 %>% is.null()) {
      stop(
        "When using a custom kernel function the 'ini_hp_0' argument is ",
        "mandatory, in order to provide the name of the hyper-parameters. ",
        "You can use the function 'hp()' to easily generate a tibble of random",
        " hyper-parameters with the desired format for initialisation."
      )
    }
  } else {
    if (ini_hp_0 %>% is.null()) {
      hp_0 <- hp(kern_0)
      cat(
        "The 'ini_hp_0' argument has not been specified. Random values of",
        "hyper-parameters for the mean process are used as initialisation.\n \n"
      )
    } else {
      hp_0 <- ini_hp_0
    }
  }

  ## Initialise the individual process' hp according to user's values
  if (kern_i %>% is.function()) {
    if (ini_hp_i %>% is.null()) {
      stop(
        "When using a custom kernel function the 'ini_hp_i' argument is ",
        "mandatory, in order to provide the name of the hyper-parameters. ",
        "You can use the function 'hp()' to easily generate a tibble of random",
        " hyper-parameters with the desired format for initialisation."
      )
    }
  } else {
    if (ini_hp_i %>% is.null()) {
      hp_i <- hp(kern_i, list_ID, common_hp = common_hp)
      cat(
        "The 'ini_hp_i' argument has not been specified. Random values of",
        "hyper-parameters for the individal processes are used as",
        "initialisation.\n \n"
      )
    } else if (!("ID" %in% names(ini_hp_i))) {
      ## Create a full tibble of common HPs if the column ID is not specified
      hp_i <- tibble::tibble(ID = list_ID, dplyr::bind_rows(ini_hp_i))
    } else if (!(all(as.character(ini_hp_i$ID) %in% as.character(list_ID)) &
      all(as.character(list_ID) %in% as.character(ini_hp_i$ID)))) {
      stop(
        "The 'ID' column in 'ini_hp_i' is different from the 'ID' of the ",
        "'data'."
      )
    }
    else {
      hp_i <- ini_hp_i
    }
  }

  ## Initialise the monitoring information
  cv <- FALSE
  logL_monitoring <- -Inf

  ## Iterate E-step and M-step until convergence
  for (i in 1:n_iter_max)
  {
    ## Track the running time for each iteration of the EM algorithm
    t_i_1 <- Sys.time()

    ## E-Step of Magma
    post <- e_step(
      db = data,
      m_0 = m_0,
      kern_0 = kern_0,
      kern_i = kern_i,
      hp_0 = hp_0,
      hp_i = hp_i,
      pen_diag = pen_diag,
    )

    ## M-Step of Magma
    new_hp <- m_step(
      db = data,
      m_0 = m_0,
      kern_0 = kern_0,
      kern_i = kern_i,
      old_hp_0 = hp_0,
      old_hp_i = hp_i,
      post_mean = post$mean,
      post_cov = post$cov,
      common_hp = common_hp,
      pen_diag = pen_diag
    )
    new_hp_0 <- new_hp$hp_0
    new_hp_i <- new_hp$hp_i

    ## In case something went wrong during the optimization
    if (any(is.na(new_hp_0)) | any(is.na(new_hp_i)) ) {
      warning(paste0("The M-step encountered an error at iteration : ", i))
      warning(
        "Training has stopped and hyper-parameters values from the ",
        "last valid iteration are returned."
      )
      break
    }

    ## Monitoring of the full log-likelihood
    new_logL_monitoring <- logL_monitoring(
      hp_0 = new_hp_0,
      hp_i = new_hp_i,
      db = data,
      m_0 = m_0,
      kern_0 = kern_0,
      kern_i = kern_i,
      post_mean = post$mean,
      post_cov = post$cov,
      pen_diag = pen_diag
    ) + 0.5 * log(det(post$cov))

    diff_logL <- new_logL_monitoring - logL_monitoring
    if (diff_logL < - 0.1) {
      warning("The likelihood descreased. Possible numerical issues.")
    }

    ## Update HPs values and the log-likelihood monitoring
    hp_0 <- new_hp_0
    hp_i <- new_hp_i
    logL_monitoring <- new_logL_monitoring

    ## Compute the convergence ratio
    eps <- diff_logL / abs(logL_monitoring)

    ## Provide monitoring information
    t_i_2 <- Sys.time()
    paste0(
      "EM algorithm, step ", i, ": ",
      difftime(t_i_2, t_i_1, units = "secs") %>% round(2),
      " seconds \n \n"
    ) %>%
      cat()

    paste0("Value of the likelihood: ",
           logL_monitoring,
           " --- Convergence ratio = ",
           eps,
           "\n \n") %>%
      cat()

    ## Check the convergence condition
    if (abs(eps) < cv_threshold) {
      cat("The EM algorithm successfully converged, training is completed.",
          "\n \n")
      cv <- TRUE
      break
    }
  }
  if (!cv & (i == n_iter_max)) {
    warning("The EM algorithm has reached the maximum number of iterations ",
            "before convergence, training might be sub-optimal \n \n")
  }

  ## Evaluate the hyper-posterior on the grid of inputs if provided
  if (!is.null(grid_inputs)) {
    if (m_0_grid %>% is.null()){
      warning("The format of 'prior_mean' doesn't allow the evaluation on the ",
      "provided grid of inputs. Thus, the hyper-posterior has only been ",
      "evaluated on training data inputs.\n \n")
    }
    else {
      cat("Start evaluating hyper-posterior distribution of the mean process",
           "on the provided grid of inputs... \n \n")

      post <- e_step(
        db = data,
        m_0 = m_0_grid,
        kern_0 = kern_0,
        kern_i = kern_i,
        hp_0 = hp_0,
        hp_i = hp_i,
        pen_diag = pen_diag,
        grid_inputs = grid_inputs
      )
      cat("Done!\n \n")
    }
  }

  ## Create a variable for directly plotting the mean process' hyper-posterior
  pred_post <- tibble::tibble(
    "Input" = post$mean %>% dplyr::pull(.data$Input),
    "Mean" = post$mean %>% dplyr::pull(.data$Output),
    "Var" = post$cov %>% diag()
  )

  t_2 <- Sys.time()
  ## Create and return the list of elements from the trained model
  list(
    "hp_0" = hp_0,
    "hp_i" = hp_i,
    "post_mean" = post$mean,
    "post_cov" = post$cov,
    "pred_post" = pred_post,
    "Converged" = cv,
    "Training_time" = difftime(t_2, t_1, units = "secs")
  ) %>%
    return()
}

#' Training new gaussian process
#'
#' @param data Database with all individuals in training set
#' @param mean_mu mean value of mean GP at timestamps (obs + pred)
#' @param cov_mu covariance value of mean GP at timestamps (obs + pred)
#' @param ini_hp_i Initial values of each parameters of the HP to start the
#'    training.
#' @param kern_i Kernel associated to individual GPs.
#'
#' @return list of trained HP
#' @export
#'
#' @examples
#' TRUE
train_new_gp <- function(data, mean_mu, cov_mu, ini_hp_i, kern_i) {
  if (is.vector(mean_mu)) {
    mean <- mean_mu
  }
  else {
    mean <- mean_mu %>%
      dplyr::filter(data$Timestamp %in% data$Timestamp) %>%
      dplyr::pull(data$Output) %>%
      as.vector()
  }
  if (length(mean) == 1) {
    mean <- rep(mean, length(data$Timestamp))
  }

  if (is.matrix(cov_mu)) {
    post_cov <- cov_mu[paste0("X", data$Timestamp), paste0("X", data$Timestamp)]
  }
  else {
    post_cov <- 0
  }

  new_hp <- optimr::opm(ini_hp_i,
    fn = logL_GP, gr = gr_GP, db = data,
    mean = mean, kern = kern_i, post_cov = post_cov,
    method = "L-BFGS-B", control = list(kkt = FALSE)
  )

  ## If something went wrong during the optimization
  if (new_hp[1, ] %>% anyNA()) {
    warning("Training encountered an error and the function returns initial
    values of the hyperparameters")
    new_hp <- ini_hp_i
  }

  tibble::as_tibble(new_hp) %>% return()
}
