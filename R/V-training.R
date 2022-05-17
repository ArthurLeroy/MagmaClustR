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
#'    - Converged: A logical value indicated whether the algorithm converged.
#'    - Training_time: Total running time of the complete training.
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
                             common_hp_k = T,
                             common_hp_i = T,
                             grid_inputs = NULL,
                             pen_diag = 1e-8,
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

  ## Check the number of cluster
  if (nb_cluster %>% is.null()) {
    nb_cluster <- 3
    ID_k <- c("K1", "K2", "K3")
    cat(
      "The number of cluster argument has not been specified. There will",
      "be 3 cluster by default. \n \n"
    )
  }

  if (!is.null(ini_hp_k)) {
    ID_k <- ini_hp_k$ID
    if (length(ID_k) != nb_cluster) {
      stop(
        "The argument 'ini_hp_k' provides hyper-parameters for a number of ",
        "clusters that is different from the 'nb_cluster' argument. "
      )
    }
  }

  ## Certify that IDs are of type 'character'
  data$ID <- data$ID %>% as.character()
  ## Extract the list of different IDs
  list_ID <- unique(data$ID)
  ## Extract the union of all reference inputs provided in the training data
  all_input <- data$Input %>%
    unique() %>%
    sort()

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
      hp_i <- hp(kern_i, list_ID = list_ID, common_hp = common_hp_i, noise = TRUE)
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
        dplyr::left_join(hp(NULL, list_ID = hp_i$ID, noise = T), by = "ID")
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
      hp_k <- hp(kern_k, list_ID = ID_k, common_hp = common_hp_k, noise = F)
      cat(
        "The 'ini_hp_k' argument has not been specified. Random values of",
        "hyper-parameters for the mean processes are used as",
        "initialisation.\n \n"
      )
    } else if (!("ID" %in% names(ini_hp_k))) {
      ## Create a full tibble of common HPs if the column ID is not specified
      hp_k <- tibble::tibble(
        ID = ID_k,
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
      m_k[[ID_k[k]]] <- prior_mean_k[[k]](all_input)
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

  if (is.null(ini_mixture)) {
    mixture <- ini_mixture(data, k = nb_cluster, name_clust = ID_k, 50)
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
      pen_diag
    )

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
    ) #%>% suppressMessages()

    new_hp_k <- new_hp$hp_k
    new_hp_i <- new_hp$hp_i

    ## In case something went wrong during the optimization
    if (any(is.na(new_hp_k)) | any(is.na(new_hp_i)) ) {
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
      post_k = post,
      m_k = m_k,
      pen_diag
    )

    diff_moni <- new_elbo_monitoring - elbo_monitoring
    if(diff_moni %>% is.nan()){diff_moni <- -Inf}

    if (diff_moni < -0.1) {
      warning("Likelihood descreased")
    }

    ## Update HPs values and the elbo monitoring
    hp_k <- new_hp_k
    hp_i <- new_hp_i
    mixture <- post$mixture
    elbo_monitoring <- new_elbo_monitoring

    ## Compute the convergence ratio
    eps <- diff_moni / abs(elbo_monitoring)
    if(eps %>% is.nan()){eps <- 1}

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
    cat("Start evaluating hyper-posterior distributions of the mean processes",
        "on the provided grid of inputs... \n \n")

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
  }
  else{
    ## Create a variable for directly plotting the mean process' hyper-posterior
    floop_pred = function(k){
      tibble::tibble(
        "Input" = post$mean[[k]] %>% dplyr::pull(.data$Input),
        "Mean" = post$mean[[k]] %>% dplyr::pull(.data$Output),
        "Var" = post$cov[[k]] %>% diag() %>% as.vector()
      ) %>%
        return()
    }
    post$pred <- sapply(ID_k, floop_pred, simplify = FALSE, USE.NAMES = TRUE)
  }


  ## Create an history list of the initial arguments of the function
  fct_args <- list(
    "data" = data,
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
                           pen_diag = 1e-8,
                           n_iter_max = 25,
                           cv_threshold = 1e-3) {
  ## Extract the observed (reference) Input
  input_obs <- data %>%
    dplyr::arrange(.data$Input) %>%
    dplyr::pull(.data$Input)

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
  if(prop_mixture %>% is.null()){
    stop(
      "The 'prop_mixture' argument is mandatory in order to provide clusters' ",
      "number, name and mixture proportions. If a MagmaClust model has been ",
      "previously trained, 'hp_k$prop_mixture' can be used as a default value."
    )
  } else{
    ## Remove the 'ID' column if present
    if ("ID" %in% names(prop_mixture)) {
      prop_mixture <- prop_mixture %>%
        dplyr::select(-.data$ID)
    }
    ## Check clusters' names
    if( !(names(prop_mixture) %>% setequal(names(hyperpost$mean))) ){
      stop(
        "The 'prop_mixture' and 'hyperpost' arguments provide different names ",
        "for clusters."
      )
    }
    ## Check that cluster's probabilities sum to 1
    if(round(sum(prop_mixture), 2) != 1){
      stop(
        "The initial probabilities in 'prop_mixture' should sum to 1 ",
        "for all clusters."
      )
    }
  }

  ## Check whether 'hyperpost' exists
  if (hyperpost %>% is.null()){
    stop(
      "The 'hyperpost' argument is necessary. Please read ?train_gp_clust."
    )
  }
  ## Check whether 'hyperpost' is evaluated on the correct Input
  if( !all(input_obs %in% hyperpost$mean[[1]]$Input) ){
    stop(
      "The 'hyperpost' argument is not evaluated on the same Input location as",
      "the data. Please see ?hyperposterior_clust to recompute the correct",
      "evaluation."
    )
  }

  ## Check whether column 'ID' exists in 'data' and 'hp' or create if necessary
  if ( ("ID" %in% names(data)) & ("ID" %in% names(hp)) ) {
    if( dplyr::n_distinct(data$ID) > 1){
      stop("The 'data' argument cannot have multiple 'ID' values.")
    }
    if( dplyr::n_distinct(hp$ID) > 1){
      stop("The 'hp' argument cannot have multiple 'ID' values.")
    }
    if( unique(data$ID) != unique(hp$ID) ){
      stop("The 'data' and 'hp' arguments have different 'ID' values.")
    }
  } else if ("ID" %in% names(data)){
    if( dplyr::n_distinct(data$ID) > 1){
      stop("The 'data' argument cannot have multiple 'ID' values.")
    }
    hp = hp %>% dplyr::mutate('ID' = unique(data$ID), .before = 1)
  } else if ("ID" %in% names(hp)){
    if( dplyr::n_distinct(hp$ID) > 1){
      stop("The 'data' argument cannot have multiple 'ID' values.")
    }
    data = data %>% dplyr::mutate('ID' = unique(hp$ID), .before = 1)
  } else{
    data = data %>% dplyr::mutate('ID' = 'ID_pred', .before = 1)
    hp = hp %>% dplyr::mutate('ID' = 'ID_pred', .before = 1)
  }

  ## Certify that IDs are of type 'character'
  data$ID <- data$ID %>% as.character()
  hp$ID <- hp$ID %>% as.character()
  ## Collect hyper-parameters' names
  list_hp <- hp %>% dplyr::select(-.data$ID) %>% names()
  ID_hp <- hp$ID %>% unique()

  ## Initialise the monitoring information
  cv <- FALSE
  logL_monitoring <- -Inf
  ## Initialisation
  mixture = prop_mixture

  for (i in 1:n_iter_max)
  {
    ## Track the running time for each iteration of the EM algorithm
    t_i_1 <- Sys.time()

    ## Format the hyper-parameters for optimisation with opm()
    par = hp %>% dplyr::select(- .data$ID)
    ## We start with a M-step to take advantage of the initial 'prop_mixture'
    ## M step
    new_hp <- optimr::opm(
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
      control = list(kkt = FALSE)
    ) %>%
      dplyr::select(list_hp) %>%
      tibble::as_tibble() %>%
      dplyr::mutate('ID' = ID_hp, .before = 1)

    ## In case something went wrong during the optimization
    if (any(is.na(new_hp)) ) {
      warning(paste0("The M-step encountered an error at iteration : ", i))
      warning(
        "Training has stopped and hyper-parameters values from the ",
        "last valid iteration are returned."
      )
      break
    }

    ## E step
    new_mixture = update_mixture(
      data,
      hyperpost$mean,
      hyperpost$cov,
      new_hp,
      kern,
      prop_mixture,
      pen_diag)

    ## Monitoring the complete log-likelihood
    new_logL_monitoring <- sum_logL_GP_clust(
      hp = new_hp,
      db = data,
      mixture = new_mixture,
      mean = hyperpost$mean,
      kern = kern,
      post_cov = hyperpost$cov,
      prop_mixture = prop_mixture,
      pen_diag = pen_diag)

    diff_logL <- new_logL_monitoring - logL_monitoring
    if(diff_logL %>% is.nan()){diff_logL <- -Inf}

    if (diff_logL < - 0.1) {
      warning("The likelihood descreased. Possible numerical issues.")
    }

    ## Update HPs values and the log-likelihood monitoring
    hp <- new_hp
    mixture <- new_mixture
    logL_monitoring <- new_logL_monitoring

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

  list("hp" = hp, "mixture" = mixture) %>% return()
}
