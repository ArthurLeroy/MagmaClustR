#' Full algorithm MAGMA
#'
#' @param db Database containing all training data from all individuals. Column: ID - input - Output.
#' @param new_db Database containing data for a new individual we want a prediction on.
#' @param timestamps Timestamps we want to predict at.
#' @param kern_i Kernel associated to individual GPs.
#' @param common_hp your common Hyperparameters
#' @param plot Boolean indicating whether we want to display a graph at the end of computations.
#' @param prior_mean Prior arbitrary value for the mean process. Optional, not needed if 'mu' is given.
#' @param kern_0 Kernel associated to the mean GPs. Optional, not needed if 'mu' is given.
#' @param list_hp Hyper-parameters for all individuals in training set. Optional, computed if NULL.
#' @param mu Database containing parameters of the mean GP at all prediction timestamps. Optional, computed if NULL.
#' @param hp_new_i Hyper-pameters for the new individual to predict. Optional, computed if NULL.
#' @param ini_hp_0 Initial values of the HP to start the training for corresponding input. Optional, not needed if 'list_hp' is given.
#' @param ini_hp_i Initial values of the HP to start the training for individual. Optional, not needed if 'list_hp' is given.
#'
#' @return predicted GP parameters | posterior mean process | all trained hyperparameters
#' @export
#'
#' @examples
full_algo = function(db, new_db, timestamps, kern_i, common_hp = T, plot = T, prior_mean,
                     kern_0 = NULL, list_hp = NULL, mu = NULL, ini_hp_0 = NULL, ini_hp_i = NULL, hp_new_i = NULL)
{
  ## If hp are not provided, train the model
  if(is.null(list_hp)){list_hp = training(db, prior_mean, ini_hp_0, ini_hp_i, kern_0, kern_i, common_hp)$hp}

  ## If mean GP (mu_0) paramaters at prediction timestamps are not provided , compute them
  if(is.null(mu)){mu = posterior_mu(db, new_db, timestamps, prior_mean, kern_0, kern_i, list_hp)}

  ## If hyperparameters of the GP for the new individuals are not provided, learn them
  ## If hyperparameters are common across individuals by hypothesis, simply pick them up from the trained model
  if(is.null(hp_new_i) & common_hp){hp_new_i = list_hp}
  else if(is.null(hp_new_i)){hp_new_i = train_new_gp(new_db, mu$mean, mu$cov, ini_hp_i, kern_i)}

  pred = pred_gp(new_db, timestamps, mean_mu = mu$mean, cov_mu = mu$cov, kern_i, hp_new_i)

  ## If True, display a plot of the resulting GP prediction
  if(plot)
  {
    (plot_gp(pred, data = new_db, data_train = db) +
       ggplot2::geom_line(data = mu$pred_GP, ggplot2::aes(x = mu$pred_GP$input, y = mu$pred_GPMean), color = 'black', linetype = 'dashed')) %>% print()
  }

  list('Prediction' = pred, 'Mean_process' = mu, 'Hyperparameters' = list('other_i' = list_hp, 'new_i' = hp_new_i)) %>%
    return()
}
