#' Full algorithm MAGMACLUST
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
#' @param new_db Database containing data for a new individual we want a prediction on.
#' @param timestamps Timestamps we want to predict at.
#' @param prior_mean Prior arbitrary value for the mean process. Optional, not needed if 'mu' is given.
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
#' @param list_hp Hyper-parameters for all individuals in training set. Optional, computed if NULL.
#' @param ini_tau_i_k Initial probability to belong to each cluster for individuals. Optional, if list_hp is given
#' @param common_hp_k A boolean indicating whether hp are common among mean GPs (for each mu_k)
#' @param common_hp_i A boolean indicating whether hp are common among individual GPs (for each y_i)
#' @param mu_k list containing parameters of the mean GPs at all prediction timestamps. Optional, computed if NULL.
#' @param ini_hp_k named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_0}, the mean process' kernel. The
#'    columns/elements should be named according to the hyper-parameters
#'    that are used in \code{kern_0}.
#' @param ini_hp_i A tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}, the individual processes' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each individual/task. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_i}.#' @param hp_new_i Hyper-pameters for the new individual to predict. Optional, computed if NULL.
#' @param hp_new_i Hyper-pameters for the new individual to predict. Optional, computed if NULL.
#'
#' @return predicted GP parameters | posterior mean process | all trained hyperparameters
#' @return
#' @export
#'
#' @examples
full_algo_clust = function(db, new_db, timestamps, kern_i, ini_tau_i_k = NULL, common_hp_k = T, common_hp_i = T,
                           prior_mean = NULL, kern_0 = NULL, list_hp = NULL, mu_k = NULL, ini_hp_k = NULL, ini_hp_i = NULL, hp_new_i = NULL)
{
  if(is.null(list_hp))
  {
    list_hp = train_magma_VEM(db, prior_mean, ini_hp_k, ini_hp_i, kern_0, kern_i, ini_tau_i_k, common_hp_k, common_hp_i)
  }

  t_pred = timestamps %>% union(unique(db$Timestamp)) %>% union(unique(new_db$Timestamp)) %>% sort()
  if(is.null(mu_k)){mu_k = posterior_mu_k(db, t_pred, prior_mean, kern_0, kern_i, list_hp)}

  if(is.null(hp_new_i))
  {
    if(common_hp_i)
    {
      hp_new_i = train_new_gp_EM(new_db, mu_k, ini_hp_i, kern_i, hp_i = list_hp)
    }
    else{hp_new_i = train_new_gp_EM(new_db, mu_k, ini_hp_i, kern_i)}
  }

  pred = pred_gp_clust(new_db, timestamps, mu_k, kern_i, hp_new_i)

  # if(plot)
  # {
  #   (pred %>% plot_gp(data = rbind(new_db, db)) + geom_point(mu, aes(Timestamp, Output, col = ID))) %>% print()
  # }

  list('Prediction' = pred, 'Mean_processes' = mu_k,
       'Hyperparameters' = list('other_i' = list_hp, 'new_i' = hp_new_i)) %>%
    return()
}
