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
#' @details The user can specify custom kernel functions for the argument
#'    \code{kern_0} and \code{kern_i}. The hyper-parameters used in the kernel
#'    should have explicit names, and be contained within the \code{hp}
#'    argument. \code{hp} should typically be defined as a named vector or a
#'    data frame. Although it is not mandatory for the \code{train_magma}
#'    function to run, gradients can be provided within kernel function
#'    definition. See for example \code{\link{se_kernel}} to create a custom
#'    kernel function displaying an adequate format to be used in Magma.
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
                        common_hp = TRUE,
                        grid_inputs = NULL,
                        pen_diag = 1e-10,
                        n_iter_max = 25,
                        cv_threshold = 1e-3,
                        fast_approx = FALSE) {

  ## Check for the correct format of the training data
  if (data %>% is.data.frame()) {
    if (!all(c("ID", "Output") %in% names(data))) {
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

  ## Remove possible missing data
  data <- data %>% tidyr::drop_na()
  ## Certify that IDs are of type 'character'
  data$ID <- data$ID %>% as.character()
  ## Extract the list of different IDs
  list_ID <- data$ID %>% unique()

  ## Get input column names
  if (!("Reference" %in% (data %>% names()))) {
    names_col <- data %>%
      dplyr::select(- c(.data$ID, .data$Output)) %>%
      names()
  } else {
    names_col <- data %>%
      dplyr::select(- c(.data$ID, .data$Output, .data$Reference)) %>%
      names()
  }

  ## Keep 6 significant digits for entries to avoid numerical errors and
  ## Add a Reference column for identification
  data <- data %>%
    purrr::modify_at(tidyselect::all_of(names_col), signif) %>%
    tidyr::unite("Reference",
                 tidyselect::all_of(names_col),
                 sep = ":",
                 remove = FALSE) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(.data$ID) %>%
    dplyr::arrange(.data$Reference, .by_group = TRUE) %>%
    dplyr::ungroup()

  ## Check that individuals do not have duplicate inputs
  if(!(setequal(data %>% dplyr::select(-.data$Output),
                data %>% dplyr::select(-.data$Output) %>% unique()
                )
       )
     ){
    stop("At least one individual have several Outputs on the same grid point.",
         " Please read ?train_magma() for further details."
    )
    }

  ## Extract the union of all reference inputs provided in the training data
  all_inputs <- data %>%
    dplyr::select(-c(.data$ID, .data$Output)) %>%
    unique()
  all_input <- all_inputs %>% dplyr::pull(.data$Reference)

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
        "hyper_prior mean function has set to be constant everywhere.\n \n"
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
  } else if (prior_mean %>% is.data.frame()) {
    if (all(c(tidyselect::all_of(names_col), "Output") %in% names(prior_mean))){
      m_0 <- prior_mean %>%
        purrr::modify_at(tidyselect::all_of(names_col), signif) %>%
        tidyr::unite("Reference",
                     tidyselect::all_of(names_col),
                     sep = ":",
                     remove = FALSE) %>%
        dplyr::filter(.data$Reference %in% all_input) %>%
        dplyr::arrange(.data$Reference) %>%
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

  ## Remove ID column if present if hp_0
  if ("ID" %in% names(hp_0)) {
    hp_0 <- hp_0[names(hp_0) != "ID"]
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
      hp_i <- hp(kern_i,
                 list_ID = list_ID,
                 common_hp = common_hp,
                 noise = TRUE
      )
      cat(
        "The 'ini_hp_i' argument has not been specified. Random values of",
        "hyper-parameters for the individal processes are used as",
        "initialisation.\n \n"
      )
    } else if (!("ID" %in% names(ini_hp_i))) {
      ## Create a full tibble of common HPs if the column ID is not specified
      hp_i <- tibble::tibble('ID' = list_ID,
                             dplyr::bind_rows(ini_hp_i)
      )
    } else if (!(all(as.character(ini_hp_i$ID) %in% as.character(list_ID)) &
                 all(as.character(list_ID) %in% as.character(ini_hp_i$ID)))) {
      stop(
        "The 'ID' column in 'ini_hp_i' is different from the 'ID' of the ",
        "'data'."
      )
    } else {
      hp_i <- ini_hp_i
    }
  }

  ## Add a 'noise' hyper-parameter if absent
  if (!("noise" %in% names(hp_i))) {
    if (common_hp) {
      hp_i <- hp_i %>% dplyr::mutate(hp(NULL, noise = T))
    } else {
      hp_i <- hp_i %>%
        dplyr::left_join(hp(NULL,list_ID = hp_i$ID,noise = T),
                         by = "ID"
        )
    }
  }

  ## Keep an history of the (possibly random) initial values of hyper-parameters
  hp_i_ini <- hp_i
  hp_0_ini <- hp_0
  ## Initialise the monitoring information
  cv <- FALSE
  logL_monitoring <- -Inf
  seq_loglikelihood <- c()

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
    if (fast_approx) {
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
        pen_diag = pen_diag
      )

      cv <- FALSE
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
    if (any(is.na(new_hp_0)) | any(is.na(new_hp_i))) {
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
      pen_diag = pen_diag
    )

    diff_logL <- new_logL_monitoring - logL_monitoring
    if (diff_logL %>% is.nan()) {
      diff_logL <- -Inf
    }

    if (diff_logL < 0) {
      warning("The likelihood descreased. Possible numerical issues.")
    }

    ## Update HPs values and the log-likelihood monitoring
    hp_0 <- new_hp_0
    hp_i <- new_hp_i
    logL_monitoring <- new_logL_monitoring

    ## Track the log-likelihood values
    seq_loglikelihood <- c(seq_loglikelihood, logL_monitoring)

    ## Compute the convergence ratio
    eps <- diff_logL / abs(logL_monitoring)
    if (eps %>% is.nan()) {
      eps <- 1
    }

    ## Provide monitoring information
    t_i_2 <- Sys.time()
    paste0(
      "EM algorithm, step ", i, ": ",
      difftime(t_i_2, t_i_1, units = "secs") %>% round(2),
      " seconds \n \n"
    ) %>%
      cat()

    paste0(
      "Value of the likelihood: ",
      logL_monitoring %>% round(5),
      " --- Convergence ratio = ",
      eps %>% round(5),
      "\n \n"
    ) %>%
      cat()

    ## Check the convergence condition
    if (abs(eps) < cv_threshold) {
      cat(
        "The EM algorithm successfully converged, training is completed.",
        "\n \n"
      )
      cv <- TRUE
      break
    }
  }
  ## Check for a prematurate ending of the EM algorithm
  if (!cv & (i == n_iter_max)) {
    warning(
      "The EM algorithm has reached the maximum number of iterations ",
      "before convergence, training might be sub-optimal \n \n"
    )
  }

  ## Evaluate the hyper-posterior on the grid of inputs if provided
  if (!is.null(grid_inputs)) {
    cat(
      "Start evaluating hyper-posterior distribution of the mean process",
      "on the provided grid of inputs... \n \n"
    )

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
  } else {
    ## Create a variable for directly plotting the mean process' hyper-posterior
    post$pred <- tibble::tibble(post$mean %>%
                                  dplyr::rename(Mean = .data$Output),
                                "Var" = post$cov %>% diag() %>% as.vector()
    )
  }

  ## Create an history list of the initial arguments of the function
  fct_args <- list(
    "data" = data %>% dplyr::select(-.data$Reference),
    "prior_mean" = prior_mean,
    "ini_hp_0" = hp_0_ini,
    "ini_hp_i" = hp_i_ini,
    "kern_0" = kern_0,
    "kern_i" = kern_i,
    "common_hp" = common_hp,
    "grid_inputs" = grid_inputs,
    "pen_diag" = pen_diag,
    "n_iter_max" = n_iter_max,
    "cv_threshold" = cv_threshold
  )

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
                     pen_diag = 1e-10) {
  ## Remove the 'ID' column if present
  if ("ID" %in% names(data)) {
    if (dplyr::n_distinct(data$ID) > 1) {
      stop(
        "Problem in the 'ID' column: different values are not allowed. ",
        "The prediction can only be performed for one individual/task."
      )
    }
    data <- data %>% dplyr::select(-.data$ID)
  }
  if (!("Reference" %in% (data %>% names()))) {
    ## Get input column names
    names_col <- data %>%
      dplyr::select(-.data$Output) %>%
      names()
  } else {
    names_col <- data %>%
      dplyr::select(- c(.data$Output, .data$Reference)) %>%
      names()
  }

  ## Keep 6 significant digits for entries to avoid numerical errors and
  ## Add a Reference column for identification
  data <- data %>%
    purrr::modify_at(tidyselect::all_of(names_col), signif) %>%
    tidyr::unite("Reference",
                 tidyselect::all_of(names_col),
                 sep = ":",
                 remove = FALSE) %>%
    tidyr::drop_na() %>%
    dplyr::arrange(.data$Reference)

  ## Extract the union of all reference inputs provided in the training data
  inputs_obs <- data %>%
    dplyr::select(-.data$Output) %>%
    unique()

  input_obs <- inputs_obs %>% dplyr::pull(.data$Reference)

  ## Check whether 'hyperpost' is provided and thus used for Magma prediction
  if (!is.null(hyperpost)) {
    cat(
      "The 'hyperpost' argument is provided. Therefore, this training is",
      "considered to be part of the prediction step in Magma. Hyper-posterior",
      "mean and covariance parameters are used in the likelihood",
      "maximisation. \n \n"
    )
    mean <- hyperpost$mean %>%
      dplyr::filter(.data$Reference %in% input_obs) %>%
      dplyr::arrange(.data$Reference) %>%
      dplyr::pull(.data$Output)

    post_cov <- hyperpost$cov[
      as.character(input_obs),
      as.character(input_obs)
    ]
  } else {
    ## Set post_cov to 0 if we are not in Magma but in a classic GP training
    post_cov <- 0

    ## Extract the values of the hyper-posterior mean at reference Input
    if (prior_mean %>% is.null()) {
      mean <- rep(0, length(input_obs))
      cat(
        "The 'prior_mean' argument has not been specified. The",
        "mean function is thus set to be 0 everywhere.\n \n"
      )
    } else if (prior_mean %>% is.vector()) {
      if (length(prior_mean) == length(input_obs)) {
        mean <- prior_mean
      } else if (length(prior_mean) == 1) {
        mean <- rep(prior_mean, length(input_obs))

        cat(
          "The provided 'prior_mean' argument is of length 1. Thus, the",
          "hyper-posterior mean function has set to be constant everywhere.",
          "\n \n"
        )
      } else {
        stop(
          "The 'prior_mean' argument is of length ", length(prior_mean),
          ", whereas the grid of training inputs is of length ",
          length(input_obs)
        )
      }
    } else if (prior_mean %>% is.function()) {
      mean <- prior_mean(inputs_obs %>%
                           dplyr::select(-.data$Reference)
      )
    } else if (prior_mean %>% is.data.frame()) {
      if (all(c("Output",
                tidyselect::all_of(names_col))
              %in% names(prior_mean))) {
        mean <- prior_mean %>%
          purrr::modify_at(tidyselect::all_of(names_col), signif) %>%
          tidyr::unite("Reference",
                       tidyselect::all_of(names_col),
                       sep = ":",
                       remove = FALSE) %>%
          dplyr::filter(.data$Reference %in% input_obs) %>%
          dplyr::arrange(.data$Reference) %>%
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
    } else {
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
    } else {
      hp <- ini_hp
    }
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

  hp_new <- stats::optim(
    hp,
    fn = logL_GP,
    gr = gr_GP,
    db = data,
    mean = mean,
    kern = kern,
    post_cov = post_cov,
    pen_diag = pen_diag,
    method = "L-BFGS-B",
    control = list(factr = 1e13, maxit = 25)
  )$par %>%
    tibble::as_tibble_row()

  ## If something went wrong during the optimization
  if (hp_new %>% is.na() %>% any()) {
    warning("Training encountered an error and the function returns initial",
            "values of the hyperparameters"
    )
    hp_new <- hp
  }
  return(hp_new)
}


#' Training MagmaClust with a Variational EM algorithm
#'
#' The hyper-parameters and the hyper-posterior distributions involved in
#' MagmaClust can be learned thanks to a VEM algorithm implemented in
#' \code{train_magmaclust}. By providing a dataset, the model hypotheses
#' (hyper-prior mean parameters, covariance kernels and number of clusters) and
#' initialisation values for the hyper-parameters, the function computes
#' maximum likelihood estimates of the HPs as well as the mean and covariance
#' parameters of the Gaussian hyper-posterior distributions of the mean
#' processes.
#'
#' @param data A tibble or data frame. Columns required: \code{ID}, \code{Input}
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
#' @param nb_cluster A number, indicating the number of clusters of
#'    individuals/tasks that are assumed to exist among the dataset.
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
#' @param ini_hp_k A tibble or data frame of hyper-parameters
#'    associated with \code{kern_k}, the mean process' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each cluster. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_k}.
#' @param ini_hp_i A tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}, the individual processes' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each individual/task. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_i}.
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
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#' @param kern_i A kernel function, associated with the individual GPs. (See
#'    details above in \code{kern_k}).
#' @param ini_mixture Initial values of the probability to belong to each
#'    cluster for each individual (\code{\link{ini_mixture}} can be used for
#'    a k-means initialisation. Used by default if NULL).
#' @param common_hp_k A boolean indicating whether hyper-parameters are common
#'    among the mean GPs.
#' @param common_hp_i A boolean indicating whether hyper-parameters are common
#'    among the individual GPs.
#' @param grid_inputs A vector, indicating the grid of additional reference
#'    inputs on which the mean processes' hyper-posteriors should be evaluated.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#' @param n_iter_max A number, indicating the maximum number of iterations of
#'    the VEM algorithm to proceed while not reaching convergence.
#' @param cv_threshold A number, indicating the threshold of the likelihood gain
#'    under which the VEM algorithm will stop. The convergence condition is
#'    defined as the difference of elbo between two consecutive steps,
#'    divided by the absolute value of the last one
#'    ( \eqn{(ELBO_n - ELBO_{n-1}) / |ELBO_n| } ).
#' @param fast_approx A boolean, indicating whether the VEM algorithm should
#'    stop after only one iteration of the VE-step. This advanced feature is
#'    mainly used to provide a faster approximation of the model selection
#'    procedure, by preventing any optimisation over the hyper-parameters.
#'
#' @details The user can specify custom kernel functions for the argument
#'    \code{kern_k} and \code{kern_i}. The hyper-parameters used in the kernel
#'    should have explicit names, and be contained within the \code{hp}
#'    argument. \code{hp} should typically be defined as a named vector or a
#'    data frame. Although it is not mandatory for the \code{train_magmaclust}
#'    function to run, gradients be can provided within kernel function
#'    definition. See for example \code{\link{se_kernel}} to create a custom
#'    kernel function displaying an adequate format to be used in
#'    MagmaClust.
#'
#' @return A list, containing the results of the VEM algorithm used in the
#'    training step of MagmaClust. The elements of the list are:
#'    - hp_k: A tibble containing the trained hyper-parameters for the mean
#'    process' kernel and the mixture proportions for each cluster.
#'    - hp_i: A tibble containing the trained hyper-parameters for the
#'    individual processes' kernels.
#'    - hyperpost: A sub-list containing the parameters of the mean processes'
#'    hyper-posterior distribution, namely:
#'      \itemize{
#'        \item mean: A list of tibbles containing, for each cluster, the
#'              hyper-posterior mean parameters evaluated at each \code{Input}.
#'        \item cov: A list of matrices containing, for each cluster, the
#'              hyper-posterior covariance parameter of the mean process.
#'        \item mixture: A tibble, indicating the mixture probabilities in each
#'              cluster for each individual.
#'      }
#'    - ini_args: A list containing the initial function arguments and values
#'    for the hyper-prior means, the hyper-parameters. In particular, if
#'    those arguments were set to NULL, \code{ini_args} allows us to retrieve
#'    the (randomly chosen) initialisations used during training.
#'    - seq_elbo: A vector, containing the sequence of ELBO values associated
#'    with each iteration.
#'    - converged: A logical value indicated whether the algorithm converged.
#'    - training_time: Total running time of the complete training.
#'
#' @export
#'
#' @examples
#' TRUE
train_magmaclust <- function(data,
                             nb_cluster = NULL,
                             prior_mean_k = NULL,
                             ini_hp_k = NULL,
                             ini_hp_i = NULL,
                             kern_k = "SE",
                             kern_i = "SE",
                             ini_mixture = NULL,
                             common_hp_k = TRUE,
                             common_hp_i = TRUE,
                             grid_inputs = NULL,
                             pen_diag = 1e-10,
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

  ##Convert all non ID columns to double (implicitly throw error if not numeric)
  data = data %>% dplyr::mutate(dplyr::across(- .data$ID, as.double))

  ## Check the number of cluster
  if (nb_cluster %>% is.null()) {
    nb_cluster <- 3
    ID_k <- c("K1", "K2", "K3")
    cat(
      "The number of cluster argument has not been specified. There will",
      "be 3 cluster by default. \n \n"
    )
  }

  ## Retrieve or create the names of the clusters
  if (!is.null(ini_hp_k)) {
    ID_k <- ini_hp_k$ID %>% unique()
    if (length(ID_k) != nb_cluster) {
      stop(
        "The argument 'ini_hp_k' provides hyper-parameters for a number of ",
        "clusters that is different from the 'nb_cluster' argument. "
      )
    }
  } else {
    ID_k <- paste0("K", 1:nb_cluster)
  }

  ## Certify that IDs are of type 'character'
  data$ID <- data$ID %>% as.character()
  ## Extract the list of different IDs
  list_ID <- data$ID %>% unique()

  ## Get input column names
  names_col <- data %>%
    dplyr::select(-c(.data$ID,.data$Output)) %>%
    names()

  ## Keep 6 significant digits for entries to avoid numerical errors and
  ## Add a Reference column for identification and sort according to it
  data <- data %>% purrr::modify_at(tidyselect::all_of(names_col),signif) %>%
    tidyr::unite("Reference",
                 tidyselect::all_of(names_col),
                 sep=":",
                 remove = FALSE) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(.data$ID) %>%
    dplyr::arrange(.data$Reference, .by_group = TRUE) %>%
    dplyr::ungroup()

  ## Check that individuals do not have duplicate inputs
  if(!(setequal(data %>% dplyr::select(-.data$Output),
                data %>% dplyr::select(-.data$Output) %>% unique() ))
     ){
    stop("At least one individual have several Outputs on the same grid point.",
         " Please read ?train_magma() for further details."
    )
  }

  ## Extract the union of all reference inputs provided in the training data
  all_input <- data %>%
    dplyr::pull(.data$Reference) %>%
    unique() %>%
    sort()

  all_inputs <- data %>%
    dplyr::select(-c(.data$ID, .data$Output)) %>%
    unique() %>%
    dplyr::arrange(.data$Reference)

  ## Initialise the individual process' HPs according to user's values
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
      hp_i <- hp(kern_i,
                 list_ID = list_ID,
                 common_hp = common_hp_i,
                 noise = TRUE
      )
      cat(
        "The 'ini_hp_i' argument has not been specified. Random values of",
        "hyper-parameters for the individual processes are used as",
        "initialisation.\n \n"
      )
    } else if (!("ID" %in% names(ini_hp_i))) {
      ## Create a full tibble of common HPs if the column ID is not specified
      hp_i <- tibble::tibble(
        ID = list_ID,
        dplyr::bind_rows(ini_hp_i)
      )
      cat(
        "No 'ID' column in the 'ini_hp_i' argument. The same hyper-parameter",
        "values have been duplicated for every 'ID' present in the 'data'.\n \n"
      )
    } else if (!(all(as.character(ini_hp_i$ID) %in% as.character(list_ID)) &
                 all(as.character(list_ID) %in% as.character(ini_hp_i$ID)))) {
      stop(
        "The 'ID' column in 'ini_hp_i' is different from the 'ID' of the ",
        "'data'."
      )
    } else {
      hp_i <- ini_hp_i
    }
  }

  ## Add a 'noise' hyper-parameter if absent
  if (!("noise" %in% names(hp_i))) {
    if (common_hp_i) {
      hp_i <- hp_i %>% dplyr::mutate(hp(NULL, noise = T))
    } else {
      hp_i <- hp_i %>%
        dplyr::left_join(hp(NULL,
                            list_ID = hp_i$ID,
                            noise = T),
                         by = "ID"
        )
    }
  }

  ## Initialise the cluster process' hp according to user's values
  if (kern_k %>% is.function()) {
    if (ini_hp_k %>% is.null()) {
      stop(
        "When using a custom kernel function the 'ini_hp_k' argument is ",
        "mandatory, in order to provide the name of the hyper-parameters. ",
        "You can use the function 'hp()' to easily generate a tibble of random",
        " hyper-parameters with the desired format for initialisation."
      )
    }
  } else {
    if (ini_hp_k %>% is.null()) {
      hp_k <- hp(kern_k,
                 list_ID = ID_k,
                 common_hp = common_hp_k,
                 noise = F
      )
      cat(
        "The 'ini_hp_k' argument has not been specified. Random values of",
        "hyper-parameters for the mean processes are used as",
        "initialisation.\n \n"
      )
    } else if (!("ID" %in% names(ini_hp_k))) {
      ## Create a full tibble of common HPs if the column ID is not specified
      hp_k <- tibble::tibble(
        'ID' = ID_k,
        dplyr::bind_rows(ini_hp_k)
      )
      cat(
        "No 'ID' column in the 'ini_hp_k' argument. The same hyper-parameter",
        "values have been duplicated for every cluster's 'ID'.\n \n"
      )
    } else {
      hp_k <- ini_hp_k
    }
  }

  ## Initialise m_k according to the value provided by the user
  m_k <- list()
  if (prior_mean_k %>% is.null()) {
    ## Create a list named by cluster with evaluation of the mean at all Input
    for (k in 1:nb_cluster) {
      m_k[[ID_k[k]]] <- rep(0, length(all_input))
    }
    cat(
      "The 'prior_mean' argument has not been specified. The hyper_prior mean",
      "function is thus set to be 0 everywhere.\n \n"
    )
  } else if (prior_mean_k[[1]] %>% is.function()) {
    ## Create a list named by cluster with evaluation of the mean at all Input
    for (k in 1:nb_cluster) {
      all_inputs %>% dplyr::select(-.data$Reference)
      m_k[[ID_k[k]]] <- prior_mean_k[[k]](all_inputs)
    }
  } else if (prior_mean_k %>% is.vector()) {
    if (length(prior_mean_k) == nb_cluster) {
      ## Create a list named by cluster with evaluation of the mean at all Input
      for (k in 1:nb_cluster) {
        m_k[[ID_k[k]]] <- rep(prior_mean_k[[k]], length(all_input))
      }
    } else if (length(prior_mean_k) == 1) {
      ## Create a list named by cluster with evaluation of the mean at all Input
      for (k in 1:nb_cluster) {
        m_k[[ID_k[k]]] <- rep(prior_mean_k, length(all_input))
      }
      cat(
        "The provided 'prior_mean' argument is of length 1. Thus, the same",
        "hyper-prior constant mean function has been set for each",
        "cluster.\n \n "
      )
    } else {
      stop(
        "The 'prior_mean_k' argument is of length ", length(prior_mean_k),
        ", whereas there are ", length(hp_k$ID), " clusters."
      )
    }
  } else {
    stop(
      "Incorrect format for the 'prior_mean_k' argument. Please read ",
      "?train_magmaclust() for details."
    )
  }

  ## Track the total training time
  t1 <- Sys.time()

  ## Initialize the monitoring information
  cv <- FALSE
  elbo_monitoring <- -Inf
  seq_elbo <- c()

  if (is.null(ini_mixture)) {
    mixture <- ini_mixture(data,
                           k = nb_cluster,
                           name_clust = ID_k,
                           50)
  }else if(is.data.frame(ini_mixture)){
    if(!all(c("ID", ID_k) %in% names(ini_mixture))){
      stop("Wrong format for ini_mixture. Make sure that the number of ",
           "clusters are the same both in 'train_magmaclust()' and ",
           "ini_mixture. Please read ?ini_mixture() for further details.")
    }else {
      mixture <- ini_mixture
    }
  }else{
    stop("The 'ini_mixture' argument must be a data frame. Please read ",
         "?ini_mixture() for further details.")
  }

  hp_k[["prop_mixture"]] <- mixture %>%
    dplyr::select(-.data$ID) %>%
    colMeans() %>%
    as.vector()

  ## Keep an history of the (possibly random) initial values of hyper-parameters
  hp_i_ini <- hp_i
  hp_k_ini <- hp_k
  mixture_ini <- mixture
  ## Iterate VE-step and VM-step until convergence
  for (i in 1:n_iter_max)
  {
    ## Track the running time for each iteration of the EM algorithm
    t_i_1 <- Sys.time()

    ## VE-Step of MagmaClust
    post <- ve_step(
      data,
      m_k,
      kern_k,
      kern_i,
      hp_k,
      hp_i,
      mixture,
      iter = i,
      pen_diag
    )

    ## Break after VE-step if we can to compute the fast approximation
    if (fast_approx) {
      ## Track the ELBO values
      seq_elbo <- elbo_monitoring_VEM(
        hp_k,
        hp_i,
        data,
        kern_i,
        kern_k,
        hyperpost = post,
        m_k = m_k,
        pen_diag
      )

      cv <- FALSE
      break
    }

    ## VM-Step of MagmaClsut
    new_hp <- vm_step(data,
                      hp_k,
                      hp_i,
                      list_mu_param = post,
                      kern_k,
                      kern_i,
                      m_k,
                      common_hp_k,
                      common_hp_i,
                      pen_diag
    ) # %>% suppressMessages()

    new_hp_k <- new_hp$hp_k
    new_hp_i <- new_hp$hp_i

    ## In case something went wrong during the optimization
    if (any(is.na(new_hp_k)) | any(is.na(new_hp_i))) {
      warning(paste0("The M-step encountered an error at iteration : ", i))
      warning(
        "Training has stopped and hyper-parameters values from the ",
        "last valid iteration are returned."
      )
      break
    }

    ## Monitoring of the elbo
    new_elbo_monitoring <- elbo_monitoring_VEM(
      new_hp_k,
      new_hp_i,
      data,
      kern_i,
      kern_k,
      hyperpost = post,
      m_k = m_k,
      pen_diag
    )

    diff_moni <- new_elbo_monitoring - elbo_monitoring
    if (diff_moni %>% is.nan()) {
      diff_moni <- -Inf
    }

    if (diff_moni < 0) {
      warning("Likelihood descreased")
    }

    ## Update HPs values and the elbo monitoring
    hp_k <- new_hp_k
    hp_i <- new_hp_i
    mixture <- post$mixture
    elbo_monitoring <- new_elbo_monitoring

    ## Track the ELBO values
    seq_elbo <- c(seq_elbo, elbo_monitoring)

    ## Compute the convergence ratio
    eps <- diff_moni / abs(elbo_monitoring)
    if (eps %>% is.nan()) {
      eps <- 1
    }

    ## Provide monitoring information
    t_i_2 <- Sys.time()
    paste0(
      "VEM algorithm, step ", i, ": ",
      difftime(t_i_2, t_i_1, units = "secs") %>% round(2),
      " seconds \n \n"
    ) %>%
      cat()

    paste0(
      "Value of the elbo: ",
      elbo_monitoring %>% round(5),
      " --- Convergence ratio = ",
      eps %>% round(5),
      "\n \n"
    ) %>%
      cat()

    ## Check the convergence condition
    if (abs(eps) < cv_threshold) {
      cat(
        "The EM algorithm successfully converged, training is completed.",
        "\n \n"
      )
      cv <- TRUE
      break
    }
    ## Check for a prematurate ending of the EM algorithm
    if (!cv & (i == n_iter_max)) {
      warning(
        "The EM algorithm has reached the maximum number of iterations ",
        "before convergence, training might be sub-optimal \n \n"
      )
    }
  }

  ## Evaluate the hyper-posterior on the grid of inputs if provided
  if (!is.null(grid_inputs)) {
    cat(
      "Start evaluating hyper-posterior distributions of the mean processes",
      "on the provided grid of inputs... \n \n"
    )

    post <- hyperposterior_clust(
      data = data,
      post$mixture,
      hp_k = hp_k,
      hp_i = hp_i,
      kern_k = kern_k,
      kern_i = kern_i,
      prior_mean_k = prior_mean_k,
      grid_inputs = grid_inputs,
      pen_diag = pen_diag
    )
    cat("Done!\n \n")
  } else {
    ## Create a variable for directly plotting the mean process' hyper-posterior
    floop_pred <- function(k) {
      tibble::tibble(post$mean[[k]],
                     "Var" =  post$cov[[k]] %>%
                       diag() %>%
                       as.vector()
      ) %>%
        dplyr::rename("Mean" = .data$Output) %>%
        dplyr::select(-.data$Reference) %>%
        return()
    }
    post$pred <- sapply(ID_k, floop_pred, simplify = FALSE, USE.NAMES = TRUE)
  }


  ## Create an history list of the initial arguments of the function
  fct_args <- list(
    "data" = data %>% dplyr::select(-.data$Reference),
    "nb_cluster" = nb_cluster,
    "prior_mean_k" = prior_mean_k,
    "ini_hp_k" = hp_k_ini,
    "ini_hp_i" = hp_i_ini,
    "kern_k" = kern_k,
    "kern_i" = kern_i,
    "ini_mixture" = mixture_ini,
    "common_hp_k" = common_hp_k,
    "common_hp_i" = common_hp_i,
    "n_iter_max" = n_iter_max,
    "pen_diag" = pen_diag,
    "cv_threshold" = cv_threshold
  )
  t2 <- Sys.time()
  ## Create and return the list of elements from the trained model
  list(
    "hp_k" = hp_k,
    "hp_i" = hp_i,
    "hyperpost" = post,
    "ini_args" = fct_args,
    "seq_elbo" = seq_elbo,
    "converged" = cv,
    "training_time" = difftime(t2, t1, units = "secs")
  ) %>%
    return()
}

#' Prediction in MagmaClust: learning new HPs and mixture probabilities
#'
#' Learning hyper-parameters and mixture probabilities of any new
#' individual/task is required in \code{MagmaClust} in the prediction procedure.
#' By providing data for the new individual/task, the hyper-posterior mean and
#' covariance parameters, the mixture proportions, and initialisation values for
#' the hyper-parameters, \code{train_gp_clust} uses an EM algorithm to compute
#' maximum likelihood estimates of the hyper-parameters and hyper-posterior
#' mixture probabilities of the new individual/task.
#'
#' @param data  A tibble or data frame. Required columns: \code{Input},
#'    \code{Output}. Additional columns for covariates can be specified.
#'    The \code{Input} column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    \code{Output} column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference \code{Input}.
#' @param prop_mixture A tibble or a named vector. Each name of column or
#'    element should refer to a cluster. The value associated with each cluster
#'    is a number between 0 and 1, corresponding to the mixture
#'    proportions.
#' @param ini_hp A tibble or data frame of hyper-parameters
#'    associated with \code{kern}, the individual process kernel.
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
#' @param hyperpost A list, containing the elements \code{mean}, \code{cov} and
#'   \code{mixture} the parameters of the hyper-posterior distributions of the
#'    mean processes. Typically, this argument should come from a previous
#'    learning using \code{\link{train_magmaclust}}, or a previous prediction
#'    with \code{\link{pred_magmaclust}}, with the argument \code{get_hyperpost}
#'    set to TRUE.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#' @param n_iter_max A number, indicating the maximum number of iterations of
#'    the EM algorithm to proceed while not reaching convergence.
#' @param cv_threshold A number, indicating the threshold of the likelihood gain
#'    under which the EM algorithm will stop.
#'
#' @return A list, containing the results of the EM algorithm used during the
#'    prediction step of MagmaClust. The elements of the list are:
#'    - hp: A tibble of optimal hyper-parameters for the new individual's GP.
#'    - mixture: A tibble of mixture probabilities for the new individual.
#'
#' @export
#'
#' @examples
#' TRUE
train_gp_clust <- function(data,
                           prop_mixture = NULL,
                           ini_hp = NULL,
                           kern = "SE",
                           hyperpost = NULL,
                           pen_diag = 1e-10,
                           n_iter_max = 25,
                           cv_threshold = 1e-3) {

  ## Get input column names
  if("Reference" %in% names(data)){
    names_col <- data %>%
      dplyr::select(- c(.data$Output, .data$Reference)) %>%
      names()
  }else{
    names_col <- data %>%
      dplyr::select(-.data$Output) %>%
      names()
  }

  if('ID' %in% names(names_col)){
    names_col <- names_col %>% dplyr::select(-.data$ID)
  }

  ## Keep 6 significant digits for entries to avoid numerical errors and
  ## Add Reference column if missing
  data <- data %>%
    purrr::modify_at(tidyselect::all_of(names_col),signif) %>%
    tidyr::unite("Reference",
                 tidyselect::all_of(names_col),
                 sep=":",
                 remove = FALSE) %>%
    tidyr::drop_na() %>%
    dplyr::arrange(.data$Reference)

  ## Extract the observed (reference) Input
  input_obs <- data %>%
    dplyr::arrange(.data$Reference) %>%
    dplyr::pull(.data$Reference)

  ## Initialise the individual process' hp according to user's values
  if (kern %>% is.function()) {
    if (ini_hp %>% is.null()) {
      stop(
        "When using a custom kernel function the 'ini_hp' argument is ",
        "mandatory, in order to provide the name of the hyper-parameters. ",
        "You can use the function 'hp()' to easily generate a tibble of random",
        " hyper-parameters with the desired format for initialisation."
      )
    } else {
      hp <- ini_hp
    }
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

  ## Check whether initial mixture probabilities are provided
  if (prop_mixture %>% is.null()) {
    stop(
      "The 'prop_mixture' argument is mandatory in order to provide clusters' ",
      "number, name and mixture proportions. If a MagmaClust model has been ",
      "previously trained, 'hp_k$prop_mixture' can be used as a default value."
    )
  } else {
    ## Remove the 'ID' column if present
    if ("ID" %in% names(prop_mixture)) {
      prop_mixture <- prop_mixture %>%
        dplyr::select(-.data$ID)
    }
    ## Check clusters' names
    if (!(names(prop_mixture) %>% setequal(names(hyperpost$mean)))) {
      stop(
        "The 'prop_mixture' and 'hyperpost' arguments provide different names ",
        "for clusters."
      )
    }
    ## Check that cluster's probabilities sum to 1
    if (round(sum(prop_mixture), 2) != 1) {
      stop(
        "The initial probabilities in 'prop_mixture' should sum to 1 ",
        "for all clusters."
      )
    }
  }

  ## Check whether 'hyperpost' exists
  if (hyperpost %>% is.null()) {
    stop(
      "The 'hyperpost' argument is necessary. Please read ?train_gp_clust."
    )
  }
  ## Check whether 'hyperpost' is evaluated on the correct Input
  if (!all(input_obs %in% hyperpost$mean[[1]]$Reference)) {
    stop(
      "The 'hyperpost' argument is not evaluated on the same Input location as",
      "the data. Please see ?hyperposterior_clust to recompute the correct",
      "evaluation."
    )
  }

  ## Check whether column 'ID' exists in 'data' and 'hp' or create if necessary
  if (("ID" %in% names(data)) & ("ID" %in% names(hp))) {
    if (dplyr::n_distinct(data$ID) > 1) {
      stop("The 'data' argument cannot have multiple 'ID' values.")
    }
    if (dplyr::n_distinct(hp$ID) > 1) {
      stop("The 'hp' argument cannot have multiple 'ID' values.")
    }
    if (unique(data$ID) != unique(hp$ID)) {
      stop("The 'data' and 'hp' arguments have different 'ID' values.")
    }
  } else if ("ID" %in% names(data)) {
    if (dplyr::n_distinct(data$ID) > 1) {
      stop("The 'data' argument cannot have multiple 'ID' values.")
    }
    hp <- hp %>% dplyr::mutate("ID" = unique(data$ID), .before = 1)
  } else if ("ID" %in% names(hp)) {
    if (dplyr::n_distinct(hp$ID) > 1) {
      stop("The 'data' argument cannot have multiple 'ID' values.")
    }
    data <- data %>% dplyr::mutate("ID" = unique(hp$ID), .before = 1)
  } else {
    data <- data %>% dplyr::mutate("ID" = "ID_pred", .before = 1)
    hp <- hp %>% dplyr::mutate("ID" = "ID_pred", .before = 1)
  }

  ## Certify that IDs are of type 'character'
  data$ID <- data$ID %>% as.character()
  hp$ID <- hp$ID %>% as.character()
  ## Collect hyper-parameters' names
  list_hp <- hp %>%
    dplyr::select(-.data$ID) %>%
    names()
  ID_hp <- hp$ID %>% unique()

  ## Initialise the monitoring information
  cv <- FALSE
  logL_monitoring <- -Inf

  ## Initialisation
  mixture <- prop_mixture

  for (i in 1:n_iter_max)
  {
    ## Track the running time for each iteration of the EM algorithm
    t_i_1 <- Sys.time()

    ## Format the hyper-parameters for optimisation
    par <- hp %>% dplyr::select(-.data$ID)

    ## We start with a M-step to take advantage of the initial 'prop_mixture'
    ## M step
    new_hp <- stats::optim(
      par = par,
      fn = sum_logL_GP_clust,
      gr = gr_sum_logL_GP_clust,
      db = data,
      mixture = mixture,
      mean = hyperpost$mean,
      kern = kern,
      post_cov = hyperpost$cov,
      pen_diag = pen_diag,
      method = "L-BFGS-B",
      control = list(factr = 1e13, maxit = 25)
    )$par %>%
      tibble::as_tibble_row() %>%
      dplyr::mutate("ID" = ID_hp, .before = 1)

    ## In case something went wrong during the optimization
    if (any(is.na(new_hp))) {
      warning(paste0("The M-step encountered an error at iteration : ", i))
      warning(
        "Training has stopped and hyper-parameters values from the ",
        "last valid iteration are returned."
      )
      break
    }

    ## E step
    new_mixture <- update_mixture(
      data,
      hyperpost$mean,
      hyperpost$cov,
      new_hp,
      kern,
      prop_mixture,
      pen_diag
    )

    ## Monitoring the complete log-likelihood
    new_logL_monitoring <- sum_logL_GP_clust(
      hp = new_hp,
      db = data,
      mixture = new_mixture,
      mean = hyperpost$mean,
      kern = kern,
      post_cov = hyperpost$cov,
      prop_mixture = prop_mixture,
      pen_diag = pen_diag
    )

    diff_logL <- new_logL_monitoring - logL_monitoring
    if (diff_logL %>% is.nan()) {
      diff_logL <- -Inf
    }

    if (diff_logL < 0) {
      warning("The likelihood descreased. Possible numerical issues.")
    }

    ## Update HPs values and the log-likelihood monitoring
    hp <- new_hp
    mixture <- new_mixture
    logL_monitoring <- new_logL_monitoring

    ## Compute the convergence ratio
    eps <- diff_logL / abs(logL_monitoring)
    if (eps %>% is.nan()) {
      eps <- 1
    }

    ## Provide monitoring information
    t_i_2 <- Sys.time()
    paste0(
      "EM algorithm, step ", i, ": ",
      difftime(t_i_2, t_i_1, units = "secs") %>% round(2),
      " seconds \n \n"
    ) %>%
      cat()

    paste0(
      "Value of the likelihood: ",
      logL_monitoring %>% round(5),
      " --- Convergence ratio = ",
      eps %>% round(5),
      "\n \n"
    ) %>%
      cat()

    ## Check the convergence condition
    if (abs(eps) < cv_threshold) {
      cat(
        "The EM algorithm successfully converged, training is completed.",
        "\n \n"
      )
      cv <- TRUE
      break
    }
    ## Check for a prematurate ending of the EM algorithm
    if (!cv & (i == n_iter_max)) {
      warning(
        "The EM algorithm has reached the maximum number of ",
        "iterations before convergence, the training might be ",
        "sub-optimal. \n \n"
      )
    }
  }

  list("hp" = hp, "mixture" = mixture) %>%
    return()
}

#' Select the optimal number of clusters
#'
#' In MagmaClust, as for any clustering method, the number K of clusters has to
#' be provided as an hypothesis of the model. This function implements a model
#' selection procedure, by maximising a variational BIC criterion, computed
#' for different values of K. A heuristic for a fast approximation of the
#' procedure is proposed as well, although the corresponding models would not
#' be properly trained.
#'
#' @param data A tibble or data frame. Columns required: \code{ID}, \code{Input}
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
#' @param fast_approx A boolean, indicating whether a fast approximation should
#'    be used for selecting the number of clusters. If TRUE, each Magma or
#'    MagmaClust model will perform only one E-step of the training, using
#'    the same fixed values for the hyper-parameters (\code{ini_hp_k} and
#'    \code{ini_hp_i}, or random values if not provided) in all models. The
#'    resulting models should not be considered as trained, but this approach
#'    provides an convenient heuristic to avoid a cumbersome model selection
#'    procedure.
#' @param grid_nb_cluster A vector of integer, corresponding to grid of values
#'    that will be tested for the number of clusters.
#' @param ini_hp_k A tibble or data frame of hyper-parameters associated with
#'    \code{kern_k}.
#' @param ini_hp_i A tibble or data frame of hyper-parameters associated with
#'    \code{kern_i}.
#' @param kern_k A kernel function associated to the mean processes.
#' @param kern_i A kernel function associated to the individuals/tasks.
#' @param plot A boolean indicating whether the plot of V-BIC values for all
#'    numbers of clusters should displayed.
#' @param ... Any additional argument that could be passed to
#'    \code{\link{train_magmaclust}}.
#'
#' @return A list, containing the results of model selection procedure for
#'    selecting the optimal number of clusters thanks to a V-BIC criterion
#'    maximisation. The elements of the list are:
#'    - best_k: An integer, indicating the resulting optimal number of clusters
#'    - seq_vbic: A vector, corresponding to the sequence of the V-BIC values
#'    associated with the models trained for each provided cluster's number in
#'    \code{grid_nb_cluster}.
#'    - trained_models: A list, named by associated number of clusters, of
#'    Magma or MagmaClust models that have been trained (or approximated if
#'    \code{fast_approx} = T) during the model selection procedure.
#'
#' @export
#'
#' @examples
#' TRUE
select_nb_cluster <- function(data,
                              fast_approx = TRUE,
                              grid_nb_cluster = 1:10,
                              ini_hp_k = NULL,
                              ini_hp_i = NULL,
                              kern_k = "SE",
                              kern_i = "SE",
                              plot = TRUE,
                              ...) {

  ## Remove possible missing data
  data <- data %>% tidyr::drop_na()
  ## Compute the number of different individuals/tasks
  nb_i <- data$ID %>% dplyr::n_distinct()

  ## Draw common initialisation for hyper-parameters for all values of K
  if (ini_hp_k %>% is.null()) {
    ini_hp_k <- hp(kern = kern_k)
  }

  if (ini_hp_i %>% is.null()) {
    hp_i <- hp(
      kern = kern_i, list_ID = unique(data$ID),
      common_hp = T, noise = T
    )
  } else {
    hp_i <- ini_hp_i
  }

  ## Initialise tracking of V-BIC values
  seq_vbic <- c()

  floop <- function(k) {
    t_1 = Sys.time()

    ## Attribute the initial common hyper-parameters to all clusters
    hp_k <- tibble::tibble("ID" = paste0("K", 1:k)) %>%
      dplyr::mutate(ini_hp_k)

    cat("Model selection: K = ", k,"\n \n")

    ## Train Magma if k = 1 and MagmaClust otherwise
    if (k == 1) {
      mod <-
        train_magma(
          data = data,
          ini_hp_0 = hp_k,
          ini_hp_i = hp_i,
          fast_approx = fast_approx,
          ...)

      ## Extract the value of the log-likelihood at convergence
      elbo <- mod$seq_loglikelihood %>% dplyr::last()
    } else {
      mod <-
        train_magmaclust(
          data = data,
          nb_cluster = k,
          ini_hp_k = hp_k,
          ini_hp_i = hp_i,
          fast_approx = fast_approx,
          ...)
      ## Extract the value of the ELBO at convergence
      elbo <- mod$seq_elbo %>% dplyr::last()
    }

    ## Define the adequate BIC penalty according to the hypotheses on HPs
    nb_hp_k <- hp_k %>%
      dplyr::select(-.data$ID) %>%
      unlist() %>%
      dplyr::n_distinct()

    nb_hp_i <- hp_i %>%
      dplyr::select(-.data$ID) %>%
      unlist() %>%
      dplyr::n_distinct()

    pen_bic <- 0.5 * (nb_hp_k + nb_hp_i + k - 1) * log(nb_i)

    mod[["V-BIC"]] <- elbo - pen_bic
    seq_vbic <<- c(seq_vbic, elbo - pen_bic)

    t_2 = Sys.time()
    ## Track training time
    paste0("Training time: ",
           difftime(t_2, t_1, units = "secs") %>% round(2),
           " seconds \n \n") %>%
      cat()

    return(mod)
  }
  mod_k <- sapply(grid_nb_cluster, floop, simplify = FALSE, USE.NAMES = TRUE)

  tib_vbic <- tibble::tibble("Nb_clust" = grid_nb_cluster, "VBIC" = seq_vbic)

  best_k <- tib_vbic %>%
    dplyr::filter(.data$VBIC == max(.data$VBIC)) %>%
    dplyr::pull(.data$Nb_clust)

  if (plot) {
    (ggplot2::ggplot(
      tib_vbic,
      ggplot2::aes(x = .data$Nb_clust, y = .data$VBIC)
    ) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::theme_classic() +
      ggplot2::xlab("Number of clusters") +
      ggplot2::ylab("Value of the variational BIC")) %>%
      print()
  }

  list("best_k" = best_k, "seq_vbic" = tib_vbic, "trained_models" = mod_k) %>%
    return()
}
