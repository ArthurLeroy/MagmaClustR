#' Training Magma with an EM algorithm
#'
#' The hyper-parameters and the hyper-posterior distribution involved in Magma
#' can be learned thanks to an EM algorithm implemented in \code{train_magma}.
#' By providing a dataset, the model hypotheses (hyper-prior mean parameter and
#' covariance kernels) and initialisation values for the hyper-parameters, the
#' function computes maximum likelihood estimates of the HPs as well as the
#' mean and covariance parameters of the Gaussian hyper-posterior distribution
#' of the mean process.
#'
#' @param data A tibble or data frame. Required columns: \code{ID}, \code{Input}
#'    , \code{Output}. Additional columns for covariates can be specified.
#'    The \code{ID} column contains the unique names/codes used to identify each
#'    individual/task (or batch of data).
#'    The \code{Input} column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    \code{Output} column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference \code{Input}.
#' @param prior_mean Hyper-prior mean parameter (m_0) of the mean GP. This
#'    argument can be specified under various formats, such as:
#'    - NULL (default). The hyper-prior mean would be set to 0 everywhere.
#'    - A number. The hyper-prior mean would be a constant function.
#'    - A vector of the same length as all the distinct Input values in the
#'     \code{data} argument. This vector would be considered as the evaluation
#'     of the hyper-prior mean function at the training Inputs.
#'    - A function. This function is defined as the hyper_prior mean.
#'    - A tibble or data frame. Required columns: Input, Output. The Input
#'     values should include at least the same values as in the \code{data}
#'     argument.
#' @param ini_hp_0 A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_0}, the mean process' kernel. The
#'    columns/elements should be named according to the hyper-parameters
#'    that are used in \code{kern_0}. If NULL (default), random values are used
#'    as initialisation.
#' @param ini_hp_i A tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}, the individual processes' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each individual/task. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_i}. Compared to \code{ini_hp_0} should contain an additional
#'    'noise' column to initialise the noise hyper-parameter of the model. If
#'     NULL (default), random values are used as initialisation.
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
#'    "PERIO" and "RQ" are also available here).
#' @param common_hp A logical value, indicating whether the set of
#'    hyper-parameters is assumed to be common to all individuals.
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
#'    ( \eqn{(LL_n - LL_n-1) / |LL_n|} ).
#' @param fast_approx A boolean, indicating whether the EM algorithm should
#'    stop after only one iteration of the E-step. This advanced feature is
#'    mainly used to provide a faster approximation of the model selection
#'    procedure, by preventing any optimisation over the hyper-parameters.
#'
#'
#' @details The user can specify custom kernel functions for the argument
#'    \code{kern_0} and \code{kern_i}. The hyper-parameters used in the kernel
#'    should have explicit names, and be contained within the \code{hp}
#'    argument. \code{hp} should typically be defined as a named vector or a
#'    data frame. Although it is not mandatory for the \code{train_magma}
#'    function to run, gradients can be provided within kernel function
#'    definition. See for example \code{\link{se_kernel}} to create a custom kernel
#'    function displaying an adequate format to be used in Magma.
#'
#' @return A list, gathering the results of the EM algorithm used for training
#'    in Magma. The elements of the list are:
#'    - hp_0: A tibble of the trained hyper-parameters for the mean
#'    process' kernel.
#'    - hp_i: A tibble of all the trained hyper-parameters for the
#'    individual processes' kernels.
#'    - hyperpost: A sub-list gathering the parameters of the mean processes'
#'    hyper-posterior distributions, namely:
#'      \itemize{
#'        \item mean: A tibble, the hyper-posterior mean parameter
#'           (\code{Output}) evaluated at each training reference \code{Input}.
#'        \item cov: A matrix, the covariance parameter for the hyper-posterior
#'           distribution of the mean process.
#'        \item pred: A tibble, the predicted mean and variance at \code{Input}
#'           for the mean process' hyper-posterior distribution under a format
#'           that allows the direct visualisation as a GP prediction.
#'      }
#'    - ini_args: A list containing the initial function arguments and values
#'    for the hyper-prior mean, the hyper-parameters. In particular, if
#'    those arguments were set to NULL, \code{ini_args} allows us to retrieve
#'    the (randomly chosen) initialisations used during training.
#'    - seq_loglikelihood: A vector, containing the sequence of log-likelihood
#'    values associated with each iteration.
#'    - converged: A logical value indicated whether the EM algorithm converged
#'    or not.
#'    - training_time: Total running time of the complete training.
#'
#' @export
#'
#' @examples
#' TRUE
train_magma <- function(data,
                        prior_mean = NULL,
                        ini_hp_0 = NULL,
                        ini_hp_i = NULL,
                        kern_0 = "SE",
                        kern_i = "SE",
                        common_hp = T,
                        grid_inputs = NULL,
                        pen_diag = 1e-8,
                        n_iter_max = 25,
                        cv_threshold = 1e-3,
                        fast_approx = FALSE) {

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
  all_input <- data$Input %>%
    unique() %>%
    sort()

  ## Initialise m_0 according to the value provided by the user
  if (prior_mean %>% is.null()) {
    m_0 <- rep(0, length(all_input))
    cat(
      "The 'prior_mean' argument has not been specified. The hyper_prior mean",
      "function is thus set to be 0 everywhere.\n \n"
    )
  }
  else if (prior_mean %>% is.vector()) {
    if (length(prior_mean) == length(all_input)) {
      m_0 <- prior_mean
    } else if (length(prior_mean) == 1) {
      m_0 <- rep(prior_mean, length(all_input))
      cat(
        "The provided 'prior_mean' argument is of length 1. Thus, the",
        "hyper_prior mean function has set to be constant everywhere.\n \n"
      )
    }
    else {
      stop(
        "The 'prior_mean' argument is of length ", length(prior_mean),
        ", whereas the grid of training inputs is of length ",
        length(all_input)
      )
    }
  }
  else if (prior_mean %>% is.function()) {
    m_0 <- prior_mean(all_input)
  }
  else if (prior_mean %>% is.data.frame()) {
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
  }
  else {
    stop(
      "Incorrect format for the 'prior_mean' argument. Please read ",
      "?hyperposterior() for details."
    )
  }

  ## Track the total training time
  t_1 <- Sys.time()

  ## Initialise the mean process' HPs according to user's values
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
      hp_i <- hp(kern_i, list_ID = list_ID, common_hp = common_hp, noise = TRUE)
      cat(
        "The 'ini_hp_i' argument has not been specified. Random values of",
        "hyper-parameters for the individal processes are used as",
        "initialisation.\n \n"
      )
    } else if (!("ID" %in% names(ini_hp_i))) {
      ## Create a full tibble of common HPs if the column ID is not specified
      hp_i <- tibble::tibble(
        ID = list_ID,
        dplyr::bind_rows(ini_hp_i)
        )
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

  ## Add a 'noise' hyper-parameter if absent
  if(!('noise' %in% names(hp_i))){
    if(common_hp){
      hp_i = hp_i %>% dplyr::mutate(hp(NULL, noise = T))
    } else{
      hp_i = hp_i %>%
        dplyr::left_join(hp(NULL, list_ID = hp_i$ID, noise = T), by = 'ID')
    }
  }

  ## Keep an history of the (possibly random) initial values of hyper-parameters
  hp_i_ini = hp_i
  hp_0_ini = hp_0
  ## Initialise the monitoring information
  cv <- FALSE
  logL_monitoring <- -Inf
  seq_loglikelihood = c()

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
      pen_diag = pen_diag
    )

    ## Break after E-step if we can to compute the fast approximation
    if(fast_approx){
      ## Track the log-likelihood values
      seq_loglikelihood <- logL_monitoring(
        hp_0 = hp_0,
        hp_i = hp_i,
        db = data,
        m_0 = m_0,
        kern_0 = kern_0,
        kern_i = kern_i,
        post_mean = post$mean,
        post_cov = post$cov,
        pen_diag = pen_diag)

      cv = FALSE
      break
    }

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

    ## Monitoring the complete log-likelihood
    new_logL_monitoring <- logL_monitoring(
      hp_0 = new_hp_0,
      hp_i = new_hp_i,
      db = data,
      m_0 = m_0,
      kern_0 = kern_0,
      kern_i = kern_i,
      post_mean = post$mean,
      post_cov = post$cov,
      pen_diag = pen_diag)

    diff_logL <- new_logL_monitoring - logL_monitoring
    if(diff_logL %>% is.nan()){diff_logL <- -Inf}

    if (diff_logL < - 0.1) {
      warning("The likelihood descreased. Possible numerical issues.")
    }

    ## Update HPs values and the log-likelihood monitoring
    hp_0 <- new_hp_0
    hp_i <- new_hp_i
    logL_monitoring <- new_logL_monitoring

    ## Track the log-likelihood values
    seq_loglikelihood = c(seq_loglikelihood, logL_monitoring)

    ## Compute the convergence ratio
    eps <- diff_logL / abs(logL_monitoring)
    if(eps %>% is.nan()){eps <- 1}

    ## Provide monitoring information
    t_i_2 <- Sys.time()
    paste0(
      "EM algorithm, step ", i, ": ",
      difftime(t_i_2, t_i_1, units = "secs") %>% round(2),
      " seconds \n \n"
    ) %>%
      cat()

    paste0("Value of the likelihood: ",
           logL_monitoring %>% round(5),
           " --- Convergence ratio = ",
           eps %>% round(5),
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
  ## Check for a prematurate ending of the EM algorithm
  if (!cv & (i == n_iter_max)) {
    warning("The EM algorithm has reached the maximum number of iterations ",
            "before convergence, training might be sub-optimal \n \n")
  }

  ## Evaluate the hyper-posterior on the grid of inputs if provided
  if (!is.null(grid_inputs)) {
      cat("Start evaluating hyper-posterior distribution of the mean process",
           "on the provided grid of inputs... \n \n")

      post <- hyperposterior(
        data = data,
        hp_0 = hp_0,
        hp_i = hp_i,
        kern_0 = kern_0,
        kern_i = kern_i,
        prior_mean = prior_mean,
        grid_inputs = grid_inputs,
        pen_diag = pen_diag
      )
      cat("Done!\n \n")
  } else{
    ## Create a variable for directly plotting the mean process' hyper-posterior
    post$pred <- tibble::tibble(
      "Input" = post$mean %>% dplyr::pull(.data$Input),
      "Mean" = post$mean %>% dplyr::pull(.data$Output),
      "Var" = post$cov %>% diag() %>% as.vector()
    )
  }

  ## Create an history list of the initial arguments of the function
  fct_args  = list('data' = data,
                  'prior_mean' = prior_mean,
                  'ini_hp_0' = hp_0_ini,
                  'ini_hp_i' = hp_i_ini,
                  'kern_0'= kern_0,
                  'kern_i' = kern_i,
                  'common_hp' = common_hp,
                  'grid_inputs' = grid_inputs,
                  'pen_diag' = pen_diag,
                  'n_iter_max' = n_iter_max,
                  'cv_threshold' = cv_threshold)

  t_2 <- Sys.time()
  ## Create and return the list of elements from the trained model
  list(
    "hp_0" = hp_0,
    "hp_i" = hp_i,
    "hyperpost" = post,
    "ini_args" = fct_args,
    "seq_loglikelihood" = seq_loglikelihood,
    "converged" = cv,
    "training_time" = difftime(t_2, t_1, units = "secs")
  ) %>%
    return()
}

#' Learning hyper-parameters of a Gaussian Process
#'
#' Learning hyper-parameters of any new individual/task in \code{Magma} is
#' required in the prediction procedure. This function can also be used to learn
#' hyper-parameters of a simple GP (just let the \code{hyperpost} argument set
#' to NULL, and use \code{prior_mean} instead). When using within \code{Magma},
#' by providing data for the new individual/task, the hyper-posterior mean and
#' covariance parameters, and initialisation values for the hyper-parameters,
#' the function computes maximum likelihood estimates of the hyper-parameters.
#'
#' @param data A tibble or data frame. Required columns: \code{Input},
#'    \code{Output}. Additional columns for covariates can be specified.
#'    The \code{Input} column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    \code{Output} column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference \code{Input}.
#' @param prior_mean Mean parameter of the GP. This argument can be
#'    specified under various formats, such as:
#'    - NULL (default). The hyper-posterior mean would be set to 0 everywhere.
#'    - A number. The hyper-posterior mean would be a constant function.
#'    - A vector of the same length as all the distinct Input values in the
#'     \code{data} argument. This vector would be considered as the evaluation
#'     of the hyper-posterior mean function at the training Inputs.
#'    - A function. This function is defined as the hyper-posterior mean.
#'    - A tibble or data frame. Required columns: Input, Output. The Input
#'     values should include at least the same values as in the \code{data}
#'     argument.
#' @param ini_hp A named vector, tibble or data frame of hyper-parameters
#'    associated with the \code{kern} of the new individual/task.
#'    The columns should be named according to the hyper-parameters that are
#'    used in \code{kern}. In cases where the model includes a noise term,
#'    \code{ini_hp} should contain an additional 'noise' column. If NULL
#'    (default), random values are used as initialisation.
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
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the²
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#' @param hyperpost A list, containing the elements 'mean' and 'cov',
#'    the parameters of the hyper-posterior distribution of the mean process.
#'    Typically, this argument should come from a previous learning using
#'    \code{\link{train_magma}}, or from the \code{\link{hyperposterior}}
#'    function. If \code{hyperpost} is provided, the likelihood that is
#'    maximised is the one involved during Magma's prediction step, and the
#'    \code{prior_mean} argument is ignored. For classic GP training, leave
#'    \code{hyperpost} to NULL.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#'
#' @return A tibble, containing the trained hyper-parameters for the kernel of
#'   the new individual/task.
#'
#' @export
#'
#' @examples
#' TRUE
train_gp <- function(data,
                     prior_mean = NULL,
                     ini_hp = NULL,
                     kern = "SE",
                     hyperpost = NULL,
                     pen_diag = 1e-8) {

  ## Extract the reference Input in the data
  input_obs <- unique(data$Input) %>% sort()

  ## Remove the 'ID' column if present
  if ("ID" %in% names(data)) {
    if(dplyr::n_distinct(data$ID) > 1){
      stop(
        "Problem in the 'ID' column: different values are not allowed. ",
        "The prediction can only be performed for one individual/task."
      )
    }
    data <- data %>%
      dplyr::select(-.data$ID)
  }

  ## Check whether 'hyperpost' is provided and thus used for Magma prediction
  if(!is.null(hyperpost)){
    cat(
      "The 'hyperpost' argument is provided. Therefore, this training is",
      "considered to be part of the prediction step in Magma. Hyper-posterior",
      "mean and covariance parameters are used in the likelihood",
      "maximisation. \n \n"
    )
    mean <-  hyperpost$mean %>%
      dplyr::filter(.data$Input %in% input_obs) %>%
      dplyr::arrange(.data$Input) %>%
      dplyr::pull(.data$Output)

    post_cov <- hyperpost$cov[
      as.character(input_obs),
      as.character(input_obs)]

  } else {
    ## Set post_cov to 0 if we are not in Magma but in a classic GP training
    post_cov = 0
    ## Extract the values of the hyper-posterior mean at reference Input
    if (prior_mean %>% is.null()) {
      mean <- rep(0, length(input_obs))
      cat(
        "The 'prior_mean' argument has not been specified. The",
        "mean function is thus set to be 0 everywhere.\n \n"
      )
    }
    else if (prior_mean %>% is.vector()) {
      if (length(prior_mean) == length(input_obs)) {
        mean <- prior_mean
      } else if (length(prior_mean) == 1) {
        mean <- rep(prior_mean, length(input_obs))

        cat(
          "The provided 'prior_mean' argument is of length 1. Thus, the",
          "hyper-posterior mean function has set to be constant everywhere.",
          "\n \n"
        )
      }
      else {
        stop(
          "The 'prior_mean' argument is of length ", length(prior_mean),
          ", whereas the grid of training inputs is of length ",
          length(input_obs)
        )
      }
    }
    else if (prior_mean %>% is.function()) {
      mean <- prior_mean(input_obs)
    }
    else if (prior_mean %>% is.data.frame()) {
      if (all(c("Output", "Input") %in% names(prior_mean))) {
        mean <- prior_mean %>%
          dplyr::filter(.data$Input %in% input_obs) %>%
          dplyr::arrange(.data$Input) %>%
          dplyr::pull(.data$Output)

        if (length(mean) != length(input_obs)) {
          stop(
            "Problem in the length of the prior mean parameter. The ",
            "'prior_mean' argument should provide an Output value for each ",
            "Input value appearing in the training data."
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
        "?train_gp() for details."
      )
    }
  }

  ## Initialise hp according to user's values
  if (kern %>% is.function()) {
    if (ini_hp %>% is.null()) {
      stop(
        "When using a custom kernel function the 'ini_hp' argument is ",
        "mandatory, in order to provide the name of the hyper-parameters. ",
        "You can use the function 'hp()' to easily generate a tibble of random",
        " hyper-parameters with the desired format for initialisation."
      )
    }
    else{hp = ini_hp}
  } else {
    if (ini_hp %>% is.null()) {
      hp <- hp(kern, noise = T)
      cat(
        "The 'ini_hp' argument has not been specified. Random values of",
        "hyper-parameters are used as initialisation.\n \n"
      )
    } else {
      hp <- ini_hp
    }
  }
  ## Extract the names of hyper-parameters
  list_hp <- hp %>% names()

  hp_new <- optimr::opm(
    hp,
    fn = logL_GP,
    gr = gr_GP,
    db = data,
    mean = mean,
    kern = kern,
    post_cov = post_cov,
    pen_diag = pen_diag,
    method = "L-BFGS-B",
    control = list(kkt = FALSE)
  ) %>%
    dplyr::select(list_hp) %>%
    tibble::as_tibble()

  ## If something went wrong during the optimization
  if (hp_new %>% is.na() %>% any()) {
    warning("Training encountered an error and the function returns initial
    values of the hyperparameters")
    hp_new <- hp
  }
  return(hp_new)
}
