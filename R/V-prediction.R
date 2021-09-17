#' Posterior mu of the cluster
#'
#' @param db data A tibble or data frame. Required columns: \code{ID}, \code{Input}
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
#' @param timestamps timestamps on which we want a prediction
#' @param kern_0 kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0).
#' @param kern_i kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i).
#' @param m_k prior value of the mean parameter of the mean GPs (mu_k). Length = 1 or nrow(db)
#' @param list_hp list of your hyperparameters
#'
#' @return Pamameters of the mean GP at timestamps chosen.
#' @export
#'
#' @examples
posterior_mu_k = function(db, grid_inputs, m_k, kern_0, kern_i, list_hp, pen_diag)
{
  #browser()
  hp_i = list_hp$hp_i
  hp_k = list_hp$hp_k
  hp_mixture = list_hp$param$hp_mixture
  t_clust = tibble::tibble('ID' = rep(hp_k$ID, each = length(grid_inputs)) , 'Input' = rep(grid_inputs, length(hp_k$ID)))
  inv_k = list_kern_to_inv(t_clust, kern_0, hp_k, pen_diag = pen_diag)
  list_inv_i = list_kern_to_inv(db, kern_i, hp_i, pen_diag = pen_diag)
  value_i = base::split(db$Output, list(db$ID))

  ## Update each mu_k parameters
  floop = function(k)
  {
    new_inv = inv_k[[k]]
    for(x in list_inv_i %>% names())
    {
      inv_i = list_inv_i[[x]]
      common_times = intersect(row.names(inv_i), row.names(new_inv))
      new_inv[common_times, common_times] = new_inv[common_times, common_times] +
        as.double(hp_mixture[k][x,]) * inv_i[common_times, common_times]
    }
    s_inv <- tryCatch(new_inv %>% chol() %>% chol2inv(),
                      error = function(e){MASS::ginv(new_inv)})

    colnames(s_inv) <- colnames(new_inv)
    rownames(s_inv) <- rownames(new_inv)

    s_inv %>% return()
  }
  cov_k = sapply(hp_k$ID, floop, simplify = FALSE, USE.NAMES = TRUE)

  floop2 = function(k)
  {
    prior_mean <- m_k[[k]]
    if(length(prior_mean) == 1){prior_mean = rep(prior_mean, ncol(inv_k[[k]]))}
    weighted_mean = inv_k[[k]] %*% prior_mean
    #row.names(weithed_mean) = row.names(inv_k[[k]])

    for(i in list_inv_i %>% names())
    {
      weighted_i = as.double(hp_mixture[k][i,]) * list_inv_i[[i]] %*% value_i[[i]]
      #row.names(weithed_i) = row.names(list_inv_i[[j]])

      common_times = intersect(row.names(weighted_i), row.names(weighted_mean))
      weighted_mean[common_times,] = weighted_mean[common_times,] + weighted_i[common_times,]
    }

    new_mean = cov_k[[k]] %*% weighted_mean %>% as.vector()
    tibble::tibble('Input' = grid_inputs, 'Output' = new_mean) %>% return()
  }
  #browser()
  mean_k = sapply(hp_k$ID, floop2, simplify = FALSE, USE.NAMES = TRUE)

  list('mean' = mean_k, 'cov' = cov_k, 'hp_mixture' = hp_mixture) %>% return()
}

#' Prediction Gaussian Process on the clustering
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
#' @param timestamps timestamps on which we want a prediction
#' @param hp_new_indiv A tibble of the Trainied results of a variation of the EM algorithm used for training
#'    in MagmaClust. Containing columns 'theta_new' and 'hp_k_mixture'.
#'    'theta_new' being a new set of hyperparameters with 1 individual.
#'    'hp_k_mixture' the probability of being in the clusters.
#'    Can be compute with the function \code{train_new_gp_EM}.
#'    train_new_gp_EM(simu_db(M=1))
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
#' @param list_mu List containing mean and cov of the K mean GPs.
#' - mean   : value of mean GP at timestamps (obs + pred) (matrix dim: timestamps x 1, with Input rownames)
#' - cov_mu : covariance of mean GP at timestamps (obs + pred) (square matrix, with Input row/colnames)
#' @param nb_cluster The number of clusters wanted.
#' @param ini_hp_k named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern}, the mean process' kernel. The
#'    columns/elements should be named according to the hyper-parameters
#'    that are used in \code{kern_k}.
#' @param ini_hp_i A tibble or data frame of hyper-parameters
#'    associated with \code{kern}, the individual processes' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each individual/task. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_i}.
#' @param trained_magmaclust A tibble of the Trainied results of a variation of the EM algorithm used for training
#'    in MagmaClust. Can be compute with the fonction \code{train_magma_VEM}.
#'    db <- simu_db()
#'    trained_magmaclust <- train_magma_VEM(db)
#'
#' @return pamameters of the gaussian density predicted at timestamps
#' @export
#'
#' @examples
#' \donttest{
#' db <- simu_db()
#' training_test = train_magma_VEM(db)
#' pred_gp_clust(simu_db(M=1), trained_magmaclust = training_test)
#'}
pred_gp_clust = function(data,
                         db_train = NULL,
                         grid_inputs = NULL,
                         list_mu = NULL,
                         kern = "SE",
                         hp_new_indiv = NULL,
                         nb_cluster = NULL,
                         ini_hp_k = NULL,
                         ini_hp_i = NULL,
                         trained_magmaclust = NULL,
                         pen_diag = 0.01)
{
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
    } else if (inputs_obs %>% names() %>% length() == 2) {
      ## Define a default grid for 'Input'
      input_pred <- rep(
        seq(min(data$Input), max(data$Input), length.out = 20),
        each = 20)
      inputs_pred <- tibble::tibble("Input" = input_pred)
      ## Add a grid for the covariate
      name_cova = inputs_obs %>% dplyr::select(-.data$Input) %>% names()
      cova = inputs_obs[name_cova]
      inputs_pred[name_cova] <- rep(
        seq(min(cova), max(cova), length.out = 20),
        times = 20)
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
  all_input <- union(input_obs, input_pred) %>% sort()
  ## Define a training database
  if(db_train %>% is.null()){db_train <- simu_db()}
  ## Initialize the individual process' hp according to user's values
  if(hp_new_indiv %>% is.null()){
    cat(
      "The 'hp_new_indiv' argument has not been specified. Random values of",
      "hyper-parameters are used as initialisation.\n \n"
    )
    if(trained_magmaclust %>% is.null()){
    trained_magmaclust <- train_magma_VEM(db_train,
                               kern_i = kern,
                               kern_k = kern,
                               nb_cluster = nb_cluster,
                               ini_hp_k = ini_hp_k,
                               ini_hp_i = ini_hp_i)
  }

    hp <- train_new_gp_EM(data, kern_i = kern, trained_magmaclust = trained_magmaclust, grid_inputs = grid_inputs)
  }
  else{
    names_a <- hp(kern) %>% names
    names_b <- hp_new_indiv$theta_new %>% names()
    if(names_a %in% names_b %>% sum == length(names_a) ||
       hp_new_indiv$hp_k_mixture %>% sum() != 1){
      cat(
        "The 'hp_new_indiv' argument has not been specified correctly.",
        " 'hp_new_indiv' must have 'theta_new' arguments related to your kernel",
        " 'hp_new_indiv' must have 'hp_k_mixture' arguments with the sum equal to 1 (probability by cluster)",
        "Random values of hyper-parameters are used as initialisation.\n \n"
      )
      if(trained_magmaclust %>% is.null()){
      trained_magmaclust <- train_magma_VEM(db_train,
                                            kern_i = kern,
                                            kern_k = kern,
                                            nb_cluster = nb_cluster,
                                            ini_hp_k = ini_hp_k,
                                            ini_hp_i = ini_hp_i)
      }

      hp <- train_new_gp_EM(data, kern_i = kern, trained_magmaclust = trained_magmaclust)
    }
  }
  hp_new = hp$theta_new
  hp_k_mixture = hp$hp_k_mixture


  ## Extract the observed Output
  yn = data %>% dplyr::pull(.data$Output)
  # ## Define the target inputs to predict
  # if(is.null(timestamps)){timestamps = seq(min(data$Input), max(data$Input), length.out = 500)}
  t_pred <- input_pred %>%
    dplyr::union(unique(db_train$Input)) %>%
    dplyr::union(unique(data$Input)) %>%
    sort()
  ## Extract the values of the hyper-posterior mean and covariance parameter at reference Input
  if(list_mu %>% is.null()){
    if(trained_magmaclust %>% is.null()
       | trained_magmaclust$prop_mixture_k %>% is.null()
       | trained_magmaclust$ini_args %>% is.null()){
      cat(
        "Neither the 'list_mu' nor 'trained_magmaclust' argument has been specified. The 'train_magma_VEM()' function",
        "(with random initialisation) has been used to learn ML estimators",
        "for the hyper-parameters associated with the 'kern' argument.\n \n"
      )
      trained_magmaclust = train_magma_VEM(db_train,
                                           kern_i = kern,
                                           kern_k = kern,
                                           nb_cluster = nb_cluster,
                                           ini_hp_k = ini_hp_k,
                                           ini_hp_i = ini_hp_i)
    }
    list_mu <- posterior_mu_k(db_train,
                              grid_inputs = t_pred,
                              trained_magmaclust$prop_mixture_k,
                              trained_magmaclust$ini_args$kern_k,
                              trained_magmaclust$ini_args$kern_i,
                              trained_magmaclust,
                              pen_diag = pen_diag)
  }


  ## Remove the noise of the hp for evaluating some of the sub-matrix
  if("noise" %in% names(hp_new)){
    hp_rm_noi = hp_new %>% dplyr::select(- .data$noise)
    noise = exp(hp_new[['noise']])
  } else {
    hp_rm_noi = hp_new
    noise = 0
  }

  ## Round the observed (reference) Input
  input_round <- input_obs %>% round(digits = 5)

  floop = function(k)
  {
    mean <- list_mu$mean[[k]] %>% round(digits = 5) %>% dplyr::distinct(.data$Input, .keep_all = TRUE)
    mean_mu_obs = mean %>% dplyr::filter(.data$Input %in% input_round) %>% dplyr::pull(.data$Output)
    timestamps_round <- input_pred %>% round(digits = 5)
    mean_mu_pred = mean %>% dplyr::filter(.data$Input %in% timestamps_round) %>% dplyr::pull(.data$Output)
    cov_mu = list_mu$cov[[k]]

    cov_tn_tn = (kern_to_cov(input_obs, kern, hp_rm_noi) + cov_mu[as.character(input_obs), as.character(input_obs)])
    inv_mat = tryCatch(cov_tn_tn %>% chol() %>% chol2inv(), error = function(e){MASS::ginv(cov_tn_tn)})
    cov_tn_t = kern_to_cov(inputs_obs, kern, hp_rm_noi, input_2 = inputs_pred) + cov_mu[as.character(input_obs), as.character(input_pred)]
    cov_t_t = kern_to_cov(inputs_pred, kern, hp_rm_noi) + cov_mu[as.character(input_pred), as.character(input_pred)]

    tibble::tibble(inputs_pred,
           'Mean' = (mean_mu_pred + t(cov_tn_t) %*% inv_mat %*% (yn - mean_mu_obs)) %>% as.vector(),
           'Var' =  (cov_t_t - t(cov_tn_t) %*% inv_mat %*% cov_tn_t) %>% diag %>% as.vector + noise,
           'hp_k_mixture' = hp_k_mixture[[k]]) %>% return()
  }
  pred = sapply(names(list_mu$mean), floop, simplify = FALSE, USE.NAMES = TRUE) %>% tibble::as_tibble()

  return(pred)
}


#' Prediction animate of the Gaussian Process on the clustering
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
#' @param timestamps timestamps on which we want a prediction
#' @param hp_new_indiv A tibble of the Trainied results of a variation of the EM algorithm used for training
#'    in MagmaClust. Containing columns 'theta_new' and 'hp_k_mixture'.
#'    'theta_new' being a new set of hyperparameters with 1 individual.
#'    'hp_k_mixture' the probability of being in the clusters.
#'    Can be compute with the function \code{train_new_gp_EM}.
#'    train_new_gp_EM(simu_db(M=1))
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
#' @param list_mu List containing mean and cov of the K mean GPs.
#' - mean   : value of mean GP at timestamps (obs + pred) (matrix dim: timestamps x 1, with Input rownames)
#' - cov_mu : covariance of mean GP at timestamps (obs + pred) (square matrix, with Input row/colnames)
#' @param nb_cluster The number of clusters wanted.
#' @param ini_hp_k named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern}, the mean process' kernel. The
#'    columns/elements should be named according to the hyper-parameters
#'    that are used in \code{kern_k}.
#' @param ini_hp_i A tibble or data frame of hyper-parameters
#'    associated with \code{kern}, the individual processes' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each individual/task. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_i}.
#' @param trained_magmaclust A tibble of the Trainied results of a variation of the EM algorithm used for training
#'    in MagmaClust. Can be compute with the fonction \code{train_magma_VEM}.
#'    db <- simu_db()
#'    trained_magmaclust <- train_magma_VEM(db)
#'
#' @return tibble of classic GP predictions but with an inscreasing number of data points considered as 'observed'
#' @export
#'
#' @examples
#' \donttest{
#' db <- simu_db()
#' training_test = train_magma_VEM(db)
#' pred_gp_clust_animate(simu_db(M=1), trained_magmaclust = training_test)
#' }
pred_gp_clust_animate = function(data,
                                 timestamps = NULL,
                                 list_mu = NULL,
                                 kern = "SE",
                                 hp_new_indiv = NULL,
                                 nb_cluster = NULL,
                                 ini_hp_k = NULL,
                                 ini_hp_i = NULL,
                                 trained_magmaclust = NULL)
{
  #browser()
  data <- data %>% dplyr::arrange(.data$Input)
  all_pred = tibble::tibble()

  for(j in 1:nrow(data))
  {
    pred_j = pred_gp_clust(data[1:j,],
                            timestamps,
                            list_mu,
                            kern,
                            hp_new_indiv,
                            nb_cluster,
                            ini_hp_k,
                            ini_hp_i,
                            trained_magmaclust)

    for(k in pred_j %>% names)
    {
      pred_j[[k]] = pred_j[[k]] %>% dplyr::mutate(Nb_data = j, Cluster = k)
    }
    #browser()
    #pred_j_all = pred_j %>% dplyr::bind_rows
    pred_j_all = dplyr::bind_rows(pred_j,pred_j)

    all_pred = all_pred %>% rbind(pred_j_all)
  }
  return(all_pred)
}

#' Prediction of the maximum of cluster
#'
#' @param hp_mixture hp_mixture
#'
#' @return Prediction of the maximum of cluster
#' @export
#'
#' @examples
pred_max_cluster = function(hp_mixture)
{
  hp_mixture %>% tibble::as_tibble %>% tidyr::unnest %>% apply(1, which.max) %>% return()
}
