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
#' @param prior_mean_k The set of hyper-prior mean parameters for the K
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
#'    ( \eqn{(ELBO_n - ELBO_{n-1}) / \abs{ELBO_n}} ). Default is \eqn{10^{-3}}.
#'
#' @return A list, containing the results of the VEM algorithm used in the
#'    training step of MagmaClust. The elements of the list are:
#'    - hp_k: A tibble containing the trained hyper-parameters for the mean
#'    process' kernel and the mixture proportions for each cluster.
#'    - hp_i: A tibble containing the trained hyper-parameters for the
#'    individual processes' kernels.
#'    - param: A sub-list containing the parameters of the mean processes'
#'    hyper-posterior distributions, namely:
#'        -> mean: A tibble containing the values of hyper-posterior's mean
#'           parameters (\code{Output}) evaluated at each \code{Input} reference
#'        -> cov: A matrix, covariance parameter of the hyper-posterior
#'           distribution of the mean process.
#'        -> mixture: A tibble, indicating the mixture probabilities in each
#'           cluster for each individual.
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
      "Incorrect format for the 'prior_mean' argument. Please read ",
      "?hyperposterior() for details."
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
      hp_k = hp_k,
      hp_i = hp_i,
      kern_k = kern_k,
      kern_i = kern_i,
      prior_mean_k = m_k,
      grid_inputs = grid_inputs,
      pen_diag = pen_diag
    )
    cat("Done!\n \n")
  }

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

  ## Create an history list of the initial arguments of the function
  fct_args <- list(
    "data" = data,
    "nb_cluster" = nb_cluster,
    "prior_mean_k" = m_k,
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
    "post_param" = post,
    "ini_args" = fct_args,
    "converged" = cv,
    "training_time" = difftime(t2, t1, units = "secs")
  ) %>%
    return()
}

#' Update by cluster the mixture star
#'
#' @param db data A tibble or data frame. Columns required: \code{ID},
#'    \code{Input}, \code{Output}.
#'    Additional columns for covariates can be specified.
#' @param hp A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_k}, the mean process' kernel. The
#'    columns/elements should be named according to the hyper-parameters
#'    that are used in \code{kern_k}.
#' @param prop_mixture_k A tibble containing the hyper-parameters associated
#'    with each individual, indicating in which cluster it belongs.
#' @param mean_k mean
#' @param cov_k cov
#' @param kern kernel A kernel function, defining the covariance structure of
#'    the GP.
#'
#' @return update mixture star
#'
#' @examples
#' TRUE
update_hp_k_mixture_star_EM <- function(db,
                                        mean_k,
                                        cov_k,
                                        kern,
                                        hp,
                                        prop_mixture_k) {
  # browser()
  all_input <- unique(db$Input) %>% sort()
  names_k <- names(mean_k)
  c_k <- 0
  mat_elbo <- rep(NA, length(names_k))

  ## Extract the specific inputs
  input_i <- db %>%
    dplyr::pull(.data$Input) %>%
    round(digits = 5)


  if ("ID" %in% names(db)) {
    db_1 <- db %>% dplyr::select(-.data$ID)
  } else {
    db_1 <- db
  }

  prop_mixture <- unlist(prop_mixture_k)

  for (k in names_k)
  {
    c_k <- c_k + 1

    round_mean <- mean_k[[k]] %>% round(digits = 5)
    unique_mean <- round_mean %>% dplyr::distinct(.data$Input, .keep_all = TRUE)

    mean <- unique_mean %>%
      dplyr::filter(.data$Input %in% input_i) %>%
      dplyr::pull(.data$Output)

    cov <- (kern_to_cov(db_1, kern, hp) +
      cov_k[[k]][as.character(input_i), as.character(input_i)])
    inv <- tryCatch(cov %>% chol() %>% chol2inv(),
      error = function(e) {
        MASS::ginv(cov)
      }
    ) %>%
      `rownames<-`(input_i) %>%
      `colnames<-`(input_i)
    # `rownames<-`(all_input) %>%
    # `colnames<-`(all_input)

    ## Classic gaussian loglikelihood
    mat_elbo[c_k] <- dmnorm(db %>% dplyr::pull(.data$Output), mean, inv,log = T)
  }
  ## We need the 'log-sum-exp' trick: exp(x - max(x)) / sum exp(x - max(x))
  ## to remain numerically stable
  mat_L <- exp(mat_elbo - max(mat_elbo))

  ((prop_mixture * mat_L) / sum(prop_mixture * mat_L)) %>%
    as.vector() %>%
    split(names_k) %>%
    tibble::as_tibble() %>%
    return()
}


#' Learning hyper-parameters of a Gaussian Process
#'
#' Learning hyper-parameters of any new individual/task in \code{magmaclust} is
#' required in the prediction procedure. When using within \code{magma},
#' by providing data for the new individual/task, the trained model
#' (hyper-posterior mean and  covariance parameters) and
#' initialization values for the hyper-parameters, the function computes
#' maximum likelihood estimates of the hyper-parameters.
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
#' @param db_train  A tibble or data frame on wich we want the training.
#'    Required columns: \code{Input}, \code{Output}.
#'    Additional columns for covariates can be specified.
#'    The \code{Input} column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    \code{Output} column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference \code{Input}.
#' @param grid_inputs The grid of inputs (reference Input and covariates) values
#'    on which the GP should be evaluated. Ideally, this argument should be a
#'    tibble or a data frame, providing the same columns as \code{data}, except
#'    'Output'. Nonetheless, in cases where \code{data} provides only one
#'    'Input' column, the \code{grid_inputs} argument can be NULL (default) or a
#'    vector. This vector would be used as reference input for prediction and if
#'    NULL, a vector of length 500 is defined, ranging between the min and max
#'    Input values of \code{data}.
#' @param nb_cluster The number of cluster wanted.
#' @param param_mu_k list of parameters for the K mean Gaussian processes
#' @param ini_hp_k A tibble or data frame of hyper-parameters
#'    associated with \code{kern_k}, the cluster processes' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each individual/task. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_i}.
#' @param ini_hp_i A tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}, the individual processes' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each individual/task. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_i}.
#' @param kern_i A kernel function, defining the covariance structure of the GP.
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
#' @param trained_magmaclust A tibble of list containing 'hp_i',
#'    a named vector, tibble or data frame of hyper-parameters associated with
#'    \code{kern_i}.
#' @param n_iter_max A number, indicating the maximum number of iterations of
#'    the EM algorithm to proceed while not reaching convergence.
#' @param cv_threshold A number, indicating the threshold of the likelihood gain
#'    under which the EM algorithm will stop.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#'
#' @return A list, containing the results of the EM algorithm used for training
#'    in MagmaClust. The elements of the list are:
#'    - theta_new :
#'    - hp_k_mixture :
#' @export
#' @examples
#' \dontrun{
#' k <- seq_len(2)
#' m_k <- c("K1" = 0, "K2" = 0, "K3" = 0)
#'
#' db <- simu_db()
#' hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
#' hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))
#'
#' ini_hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))
#' old_mixture <- MagmaClustR:::ini_mixture(
#'   db = db, k = length(k),
#'   nstart = 50
#' )
#'
#' training_test <- train_magmaclust(db)
#'
#' timestamps <- seq(0.01, 10, 0.01)
#' mu_k <- hyperposterior_clust(db, timestamps, m_k, "SE", "SE", training_test,
#'   pen_diag = 0.01
#' )
#'
#' train_new_gp_EM(simu_db(M = 1),
#'   param_mu_k = mu_k, ini_hp_i = ini_hp_i,
#'   kern_i = "SE", trained_magmaclust = training_test
#' )
#'
#' ###########################
#' db <- simu_db()
#' training_test <- train_magmaclust(db)
#' train_new_gp_EM(simu_db(M = 1), trained_magmaclust = training_test)
#'
#' ##########################
#' train_new_gp_EM(simu_db(M = 1))
#' }
train_new_gp_EM <- function(data,
                            db_train = NULL,
                            grid_inputs = NULL,
                            nb_cluster = NULL,
                            param_mu_k = NULL,
                            ini_hp_k = NULL,
                            ini_hp_i = NULL,
                            kern_i = "SE",
                            trained_magmaclust = NULL,
                            n_iter_max = 25,
                            cv_threshold = 1e-3,
                            pen_diag = 1e-8) {
  # browser()
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
      input_pred <- c(
        seq(min(data$Input),
          max(data$Input),
          length.out = 500
        ),
        data$Input
      ) %>%
        unique()
      inputs_pred <- tibble::tibble("Input" = input_pred)
    } else if (inputs_obs %>% names() %>% length() == 2) {
      ## Define a default grid for 'Input'
      input_pred <- rep(
        c(
          seq(min(data$Input),
            max(data$Input),
            length.out = 20
          ),
          data$Input
        ) %>%
          unique(),
        each = 20
      )
      inputs_pred <- tibble::tibble("Input" = input_pred)
      ## Add a grid for the covariate
      name_cova <- inputs_obs %>%
        dplyr::select(-.data$Input) %>%
        names()
      cova <- inputs_obs[name_cova]
      inputs_pred[name_cova] <- rep(
        c(seq(min(cova), max(cova), length.out = 20), data$Covariate) %>% unique(),
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
  ## Define the union of all distinct reference Input
  all_input <- union(input_obs, input_pred) %>% sort()
  ## Initialise the individual process' hp according to user's values
  if (kern_i %>% is.function()) {
    if (ini_hp_i %>% is.null()) {
      stop(
        "When using a custom kernel function the 'ini_hp_i' argument is ",
        "mandatory, in order to provide the name of the hyper-parameters. ",
        "You can use the function 'hp()' to easily generate a tibble of random",
        " hyper-parameters with the desired format for initialisation."
      )
    } else {
      hp <- ini_hp_i
    }
  } else {
    if (ini_hp_i %>% is.null()) {
      hp <- hp(kern_i, noise = T, list_ID = unique(data$ID))
      cat(
        "The 'ini_hp_i' argument has not been specified. Random values of",
        "hyper-parameters are used as initialisation.\n \n"
      )
    } else {
      hp <- ini_hp_i
    }
  }
  hp <- hp %>% dplyr::slice(1)

  ## Extract the names of hyper-parameters
  list_hp_new <- hp %>% names()
  ## Extract the values of the hyper-posterior mean and covariance parameters
  ## at reference Input
  trained_hp_i <- trained_magmaclust$hp_i
  ## Define a training database
  if (db_train %>% is.null()) {
    db_train <- simu_db()
  }
  ## Define a default prediction grid
  t_pred <- input_pred %>%
    dplyr::union(unique(db_train$Input)) %>%
    dplyr::union(unique(data$Input)) %>%
    sort()
  if (param_mu_k %>% is.null()) {
    if (trained_magmaclust %>% is.null() |
      trained_magmaclust$prop_mixture_k %>% is.null() |
      trained_magmaclust$ini_args %>% is.null()) {
      cat(
        "Neither the 'param_mu_k' nor 'trained_magmaclust' argument
        has been specified. The 'train_magmaclust()' function",
        "(with random initialisation) has been used to learn ML estimators",
        "for the hyper-parameters associated with the 'kern' argument.\n \n"
      )
      trained_magmaclust <- train_magmaclust(db_train,
        kern_i = kern_i,
        kern_k = kern_i,
        nb_cluster = nb_cluster,
        ini_hp_k = ini_hp_k,
        ini_hp_i = ini_hp_i
      )
    }
    param_mu_k <- hyperposterior_clust(db_train,
      grid_inputs = t_pred,
      trained_magmaclust$prop_mixture_k,
      trained_magmaclust$ini_args$kern_k,
      trained_magmaclust$ini_args$kern_i,
      trained_magmaclust,
      pen_diag = pen_diag
    )
  }
  mean_mu_k <- param_mu_k$mean
  cov_mu_k <- param_mu_k$cov
  ## Remove the ID column
  hp_1 <- param_mu_k$mixture %>% dplyr::select(-.data$ID)

  ## Initialize the monitoring information
  cv <- FALSE
  prop_mixture_k <- lapply(hp_1, function(x) Reduce("+", x) / length(x))
  ## Extract the reference Input in the data
  all_input <- unique(data$Input) %>% sort()

  if (is.null(trained_hp_i)) {
    for (i in 1:n_iter_max)
    {
      ## Track the running time for each iteration of the EM algorithm
      t_i_1 <- Sys.time()

      ## E step
      hp_k_mixture <- update_hp_k_mixture_star_EM(
        data,
        mean_mu_k,
        cov_mu_k,
        kern_i,
        hp,
        prop_mixture_k
      )

      ## M step
      LL_GP <- function(hp, data, kern_i) {
        inputs <- data %>%
          dplyr::pull(.data$Input) %>%
          unique() %>%
          round(digits = 10)

        floop <- function(k) {
          round_mean <- mean_mu_k[[k]] %>% round(digits = 10)

          mean <- round_mean %>%
            dplyr::filter(.data$Input %in% inputs) %>%
            dplyr::pull(.data$Output)
          cov <- (kern_to_cov(data$Input, kern_i, hp) +
            cov_mu_k[[k]][as.character(inputs), as.character(inputs)])
          inv <- tryCatch(cov %>% chol() %>% chol2inv(),
            error = function(e) {
              MASS::ginv(cov)
            }
          ) %>%
            `rownames<-`(all_input) %>%
            `colnames<-`(all_input)

          (data$Input - hp_k_mixture[[k]] *
            dmnorm(data$Output, mean, inv, log = T)) %>%
            return()
        }
        sapply(names(mean_mu_k), floop) %>%
          sum() %>%
          return()
      }
      ## Extract the hyper-parameters associated
      par_i <- hp %>%
        dplyr::select(-.data$ID)

      new_hp <- optimr::opm(
        par_i,
        LL_GP,
        data = data,
        kern = kern_i,
        method = "L-BFGS-B",
        control = list(kkt = FALSE)
      ) %>%
        dplyr::select(hp %>%
          dplyr::select(-.data$ID) %>%
          names()) %>%
        tibble::as_tibble() %>%
        dplyr::mutate("ID" = unique(hp$ID), .before = 1)

      if (new_hp %>% anyNA(recursive = T)) {
        print(paste0("The M-step encountered an error at iteration : ", i))
        print("Training has stopped and the function returns
              values from the last valid iteration")
        break
      }


      ## Provide monitoring information
      t_i_2 <- Sys.time()
      paste0(
        "EM algorithm on the trainig of new, step ", i, ": ",
        difftime(t_i_2, t_i_1, units = "secs") %>% round(2),
        " seconds \n \n"
      ) %>%
        cat()


      names_hp <- new_hp %>%
        dplyr::select(-.data$ID) %>%
        names()
      eps <- (new_hp[names_hp] - hp[names_hp]) %>%
        abs() %>%
        sum()

      print("Value of the hp mixture: ")
      print(hp_k_mixture)
      paste0(
        "Convergence ratio = ",
        eps,
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
      hp <- new_hp
    }
  } else {
    new_hp <- trained_hp_i %>% dplyr::slice(1)

    hp_k_mixture <- update_hp_k_mixture_star_EM(
      data,
      mean_mu_k,
      cov_mu_k,
      kern_i,
      new_hp,
      prop_mixture_k
    )
  }

  ## If something went wrong during the optimization
  if (new_hp %>% is.na() %>% any()) {
    warning("Training encountered an error and the function returns initial
    values of the hyperparameters")
    new_hp <- hp
  }

  list("theta_new" = new_hp, "hp_k_mixture" = hp_k_mixture) %>% return()
}
