#' Posterior mu
#'
#' @param db matrix of data columns required ('input', 'Output')
#' @param new_db Database containing data for a new individual we want a prediction on.
#' @param timestamps timestamps on which we want a prediction
#' @param m_0 prior mean value of the mean GP (scalar value or vector of same length as 'timestamps')
#' @param kern_0 kernel associated to the covariance function of the mean GP
#' @param hp_0 Set of hyperparameters of the mean GP
#' @param hp_i Set of hyperparameters for each indicidual GPs
#' @param kern_i Kernel associated to individual GPs.
#'
#' @return pamameters of the mean GP at timestamps
#' @export
#'
#' @examples
posterior_mu = function(db, new_db, timestamps, m_0, kern_0, kern_i, hp_0,hp_i)
{
  t_pred = timestamps %>% union(unique(db$input)) %>% union(unique(new_db$input)) %>% sort()
  ## Mean GP (mu_0) is noiseless and thus has only 2 hp. We add a penalty on diag for numerical stability
  pen_diag = 0.1 #sapply(hp_i, function(x) x[[3]]) %>% mean

  inv_0 = kern_to_inv(t_pred, kern_0, hp_0, pen_diag)
  inv_i = kern_to_inv(db, kern_i, hp_i)
  value_i = base::split(db$Output, list(db$ID))

  new_inv = update_inv(prior_inv = inv_0, list_inv_i = inv_i)
  new_cov = tryCatch(solve(new_inv), error = function(e){MASS::ginv(new_inv)}) ## fast or slow matrix inversion if singular
  rownames(new_cov) = rownames(new_inv)
  colnames(new_cov) = colnames(new_inv)

  weighted_mean = update_mean(prior_mean = m_0, prior_inv = inv_0, list_inv_i = inv_i, list_value_i = value_i)
  new_mean = (new_cov %*% weighted_mean) %>% as.vector

  #names(mean_mu) = paste0('X', t_mu)
  list('mean' = tibble::tibble('input' = t_pred, 'Output' = new_mean) , 'cov' = new_cov,
       'pred_GP' = tibble::tibble('input' = t_pred, 'Mean' = new_mean, 'Var' = diag(new_cov))) %>% return()
}

#' Prediction Gaussian Process
#'
#' @param db tibble of data columns required ('input', 'Output')
#' @param timestamps timestamps on which we want a prediction
#' @param mean_mu mean value of mean GP at timestamps (obs + pred) (matrix dim: timestamps x 1, with Input rownames)
#' @param cov_mu covariance of mean GP at timestamps (obs + pred) (square matrix, with Input row/colnames)
#' @param hp list of hyperparameters for the kernel of the GP
#' @param kern kernel associated to the covariance function of the GP
#'
#' @return pamameters of the gaussian density predicted at timestamps
#' @export
#'
#' @examples
pred_gp = function(db, timestamps = NULL, mean_mu = 0, cov_mu = NULL, kern, hp)
{
  tn = db %>% dplyr::pull(db$input)
  #input = db %>% dplyr::pull(Input)
  input = paste0('X', db$input)
  yn = db %>% dplyr::pull(db$Output)

  ## Define a default prediction grid
  if(is.null(timestamps)){timestamps = seq(min(tn), max(tn), length.out = 500)}
  input_t = paste0('X', timestamps)
  all_times = union(tn,timestamps)

  if(is.null(cov_mu))
  { ## Case of standard GP regression. Without trained posterior mean process.
    cov_mu = matrix(0, length(all_times), length(all_times),
                    dimnames = list(paste0('X', all_times), paste0('X', all_times)))
  }
  if(length(mean_mu) == 1)
  { ## If the provided mean is a constant function
    mean_mu_obs = rep(mean_mu, length(tn))
    mean_mu_pred = rep(mean_mu, length(timestamps))
  }
  else
  { ## If the provided mean has defined values at timestamps, typically from training. Format : input, Output
    mean_mu_obs = mean_mu %>%  dplyr::filter(db$input %in% tn) %>% dplyr::pull(db$Output)
    mean_mu_pred = mean_mu %>% dplyr::filter(db$input %in% timestamps) %>% dplyr::pull(db$Output)
  }

  cov_tn_tn = kern_to_cov(tn, kern, hp) + cov_mu[input, input]
  inv_mat = tryCatch(solve(cov_tn_tn), error = function(e){MASS::ginv(cov_tn_tn)})
  #cov_tn_t = kern(mat_dist(tn, timestamps), theta) + cov_mu[input,input_t]
  cov_tn_t = kern(tn,timestamps,hp) + cov_mu[input,input_t]


  cov_t_t = kern_to_cov(timestamps, kern, hp) + cov_mu[input_t ,input_t]

  tibble::tibble('input' = timestamps,
                 'Mean' = (mean_mu_pred + t(cov_tn_t) %*% inv_mat %*% (yn - mean_mu_obs)) %>% as.vector(),
                 'Var' = (cov_t_t - t(cov_tn_t) %*% inv_mat %*% cov_tn_t) %>% diag()) %>% return()
}

#' Prediction Gaussian Process Animate
#'
#' Function used to generate data compatible with a GIF ploting of the results
#'
#' @param db Your database, the same as for a classic GP prediction
#' @param timestamps Timestamps we want to predict at.
#' @param mean_mu mean value of mean GP at timestamps (obs + pred)
#' @param cov_mu covariance of mean GP at timestamps (obs + pred)
#' @param hp list of hyperparameters of your model
#' @param kern kernel of your choice use with your hyperparameters
#'
#' @return tibble of classic GP predictions but with an inscreasing number of data points considered as 'observed'
#' @export
#'
#' @examples
pred_gp_animate = function(db, timestamps = NULL, mean_mu = 0, cov_mu = NULL,
                           kern, hp)
{
  db %>% dplyr::arrange(db$input)
  all_pred = tibble::tibble()

  if(is.null(timestamps)){timestamps = seq(min(db$input), max(db$input), length.out = 500)}

  for(j in 1:nrow(db))
  {
    pred_j = pred_gp(db[1:j,], timestamps, mean_mu, cov_mu, kern, hp) %>% dplyr::mutate(Nb_data = j)
    all_pred = all_pred %>% rbind(pred_j)
  }
  return(all_pred)
}
