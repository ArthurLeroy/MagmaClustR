#' Full algorithm MAGMACLUST
#'
#' @param data_train A tibble or data frame. Required columns: \code{ID}, \code{Input}
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
#' @param data_obs Database containing data for a new individual we want a prediction on.
#' @param grid_inputs The grid of inputs (reference Input and covariates) values
#'    on which the GP should be evaluated. Ideally, this argument should be a
#'    tibble or a data frame, providing the same columns as \code{data}, except
#'    'Output'. Nonetheless, in cases where \code{data} provides only one
#'    'Input' column, the \code{grid_inputs} argument can be NULL (default) or a
#'    vector. This vector would be used as reference input for prediction and if
#'    NULL, a vector of length 500 is defined, ranging between the min and max
#'    Input values of \code{data}.
#' @param prior_mean Prior arbitrary value for the mean process. Optional, not needed if 'mu' is given.
#' @param kern_k A kernel function, associated with the mean GP.
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
#' @param list_hp Hyper-parameters for all individuals in training set. Optional, computed if NULL.
#' @param ini_hp_mixture Initial probability to belong to each cluster for individuals. Optional, if list_hp is given
#' @param common_hp_k A boolean indicating whether hp are common among mean GPs (for each mu_k)
#' @param common_hp_i A boolean indicating whether hp are common among individual GPs (for each y_i)
#' @param mu_k list containing parameters of the mean GPs at all prediction timestamps. Optional, computed if NULL.
#' @param ini_hp_k named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_k}, the mean process' kernel. The
#'    columns/elements should be named according to the hyper-parameters
#'    that are used in \code{kern_k}.
#' @param ini_hp_i A tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}, the individual processes' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each individual/task. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_i}.#' @param hp_new_i Hyper-pameters for the new individual to predict. Optional, computed if NULL.
#' @param hp_new_i Hyper-pameters for the new individual to predict. Optional, computed if NULL.
#' @param nb_cluster The number of clusters wanted.
#' @param n_iter_max A number, indicating the maximum number of iterations of
#'    the EM algorithm to proceed while not reaching convergence.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#' @param cv_threshold A number, indicating the threshold of the likelihood gain
#'    under which the EM algorithm will stop. The convergence condition is
#'    defined as the difference of likelihoods between two consecutive steps,
#'    divided by the absolute value of the last one
#'    ( (LL_n - LL_n-1) / |LL_n| ).
#' @param plot A logical value, indicating whether a plot of the results is
#'    automatically displayed.
#' @param prior_mean Prior mean parameter of the K mean GPs (mu_k).
#'
#' @return predicted GP parameters | posterior mean process | all trained hyperparameters
#' @return
#' @export
#'
#' @examples
#' data_train <- simu_db()
#' data_obs <- simu_db(M=1)
#' training_test <- train_magma_VEM(data_train)
#' full_algo_clust(data_train, data_obs, list_hp = training_test)
full_algo_clust = function(data_train,
                           data_obs,
                           grid_inputs = NULL,
                           kern_i = "SE",
                           kern_k = "SE",
                           ini_hp_mixture = NULL,
                           common_hp_k = T,
                           common_hp_i = T,
                           prior_mean = NULL,
                           list_hp = NULL,
                           mu_k = NULL,
                           ini_hp_k = NULL,
                           ini_hp_i = NULL,
                           hp_new_i = NULL,
                           nb_cluster = NULL,
                           n_iter_max = 25,
                           pen_diag = 0.01,
                           cv_threshold = 1e-3,
                           plot = TRUE)
{
  if(is.null(list_hp))
  {
    list_hp = train_magma_VEM(data_train,
                              nb_cluster = nb_cluster,
                              prior_mean_k = prior_mean,
                              ini_hp_k = ini_hp_k,
                              ini_hp_i = ini_hp_i,
                              kern_k = kern_k,
                              kern_i = kern_i,
                              ini_hp_mixture = ini_hp_mixture,
                              common_hp_k = common_hp_k,
                              common_hp_i = common_hp_i,
                              n_iter_max = n_iter_max,
                              pen_diag = pen_diag,
                              cv_threshold = cv_threshold)
  }

  if(prior_mean %>% is.null()){prior_mean <- list_hp$prop_mixture_k}
  if(is.null(hp_new_i))
  {
    if(common_hp_i)
    {
      hp_new_i = train_new_gp_EM(data_obs,
                                 db_train = data_train,
                                 grid_inputs = grid_inputs,
                                 nb_cluster = nb_cluster,
                                 param_mu_k = mu_k,
                                 ini_hp_k = ini_hp_k,
                                 ini_hp_i = ini_hp_i,
                                 kern_i = kern_i,
                                 trained_magmaclust = list_hp,
                                 n_iter_max = n_iter_max,
                                 cv_threshold = cv_threshold)

    }
    else{
      hp_new_i = train_new_gp_EM(data_obs,
                                 db_train = data_train,
                                 grid_inputs = grid_inputs,
                                 nb_cluster = nb_cluster,
                                 param_mu_k = mu_k,
                                 ini_hp_k = ini_hp_k,
                                 ini_hp_i = ini_hp_i,
                                 kern_i = kern_i,
                                 trained_magmaclust = NULL,
                                 n_iter_max = n_iter_max,
                                 cv_threshold = cv_threshold)}
  }

  #browser()
  pred = pred_magma_clust(data_obs,
                       data_train = data_train,
                       grid_inputs = grid_inputs,
                       list_mu = mu_k,
                       kern = kern_i,
                       hp_new_indiv = hp_new_i,
                       nb_cluster = nb_cluster,
                       ini_hp_k = ini_hp_k,
                       ini_hp_i = ini_hp_i,
                       trained_magmaclust = list_hp)

  ## Display the graph of the prediction if expected
  if(plot){plot_magma_clust(pred,
                         data = data_obs,
                         data_train = data_train) %>%
      print()

    }



  list('Prediction' = pred, 'Mean_processes' = mu_k,
       'Hyperparameters' = list('other_i' = list_hp, 'new_i' = hp_new_i)) %>%
    return()
}
