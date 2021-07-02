#' Log Likelihood function of the Gaussian Process
#'
#' @param hp The hyperparameters of your kernel
#' @param db tibble containing values we want to compute. Required columns : input, Output
#' @param mean mean of the GP at corresponding inputs
#' @param kern kernel used to compute the covariance matrix at corresponding inputs
#' @param new_cov posterior covariance matrix of the mean GP (mu_0). Used to compute correction term (cor_term)
#'
#' @return value of the Gaussian log-likelihood for one GP as it appears in the model
#' @export
#'
#' @examples
#'logL_GP(tibble::tibble(sigma = 1, lengthscale = 0.5),
#'tibble::tibble(input= 1:5,Output= 2:6),rep(3,5),kernel_sqrd_exp,0)
#'
logL_GP<- function(hp, db, mean, kern, new_cov)
{
  cov = kern_to_cov(db$input, kern, hp) + new_cov
  inv = tryCatch(solve(cov), error = function(e){MASS::ginv(cov)})
  (-mvtnorm::dmvnorm(db$Output, mean, inv , log = T)) %>%  return()
}


#' Modified Log Likelihood function of the Gaussian Process
#'
#' @param hp vector or list of parameters of the kernel
#' @param db tibble containing values we want to compute logL on. Required columns : input, Output
#' @param mean mean of the GP at corresponding inputs
#' @param kern kernel used to compute the covariance matrix at corresponding inputs
#' @param new_cov posterior covariance matrix of the mean GP (mu_0). Used to compute correction term (cor_term)
#' @param pen_diag value of the penalization of the diagonal
#'
#' @return value of the modified Gaussian log-likelihood for one GP as it appears in the model
#' @export
#'
#' @examples
#'logL_GP(tibble::tibble(sigma = 1, lengthscale = 0.5),
#'tibble::tibble(input= 1:5,Output= 2:6),rep(3,5),kernel_sqrd_exp,0)
#'
logL_GP_mod = function(hp, db, mean, kern, new_cov, pen_diag = NULL)
{
  t1 = Sys.time()
  if(length(mean) == 1){mean = rep(mean, nrow(db))} ## mean is equal for all inputs
  ## Mean GP (mu_0) is noiseless and thus has only 2 hp. We add a penalty on diag for numerical stability
  if(length(hp) != 3){hp$sigma<- pen_diag}


  inv =  kern_to_inv(db$input, kern, hp)

  LL_norm = - mvtnorm::dmvnorm(db$Output, mean, inv, log = T) ## classic gaussian loglikelihood
  cor_term =  0.5 * (inv * new_cov) %>% sum() ## correction term (0.5 * Trace(inv %*% new_cov))
  t2 = Sys.time()
  #print(paste0('LogL_0 iteration ', t2 - t1))
  return(LL_norm + cor_term)
}

#' Modified Gaussian log-likelihood for the sum of all indiv
#'
#' @param hp vector of common hyperparameters for all individuals.
#' @param db tibble of data. Required columns : ID, input, Output
#' @param mean mean of the GP at union of observed inputs
#' @param kern kernel used to compute the covariance matrix at corresponding inputs
#' @param new_cov posterior covariance matrix of the mean GP (mu_0). Used to compute correction term (cor_term)
#'
#' @return value of the modified Gaussian log-likelihood for the sum of all indiv with same HPs
#' @export
#'
#' @examples
#'logL_GP(tibble::tibble(sigma = 1, lengthscale = 0.5),
#'tibble::tibble(ID = 1:5, input= 1:5,Output= 2:6),rep(3,5),kernel_sqrd_exp,0)
logL_GP_mod_common_hp = function(hp, db, mean, kern, new_cov)
{
  LL_norm = 0
  cor_term = 0
  t_i_old = NULL
  t1 = Sys.time()
  for(i in unique(db$ID))
  {
    t_i = db %>% dplyr::filter(db$ID == i) %>% dplyr::pull(db$input)
    input_i = paste0('X', t_i)
    y_i = db %>% dplyr::filter(db$ID == i) %>% dplyr::pull(db$Output)

    if( !identical(t_i, t_i_old) )
    { ## We update the inverse cov matrix only if necessary (if different inputs)
      inv =  kern_to_inv(t_i, kern, hp)
    }

    LL_norm = LL_norm - mvtnorm::dmvnorm(y_i, mean %>% dplyr::filter(db$input %in% t_i) %>% dplyr::pull(db$Output), inv, log = T)
    cor_term = cor_term + 0.5 * (inv * new_cov[input_i, input_i]) %>% sum()  ##(0.5 * Trace(inv %*% new_cov))

    t_i_old = t_i
  }
  t2 = Sys.time()
  #print(paste0('LogL_i iteration ', t2 - t1))
  return(LL_norm + cor_term)
}

#' Title
#'
#' @param db tibble containing values we want to compute logL on. Required columns : ID, input, Output
#' @param kern_i kernel used to compute the covariance matrix of individuals GP at corresponding inputs (Psi_i)
#' @param kern_0 kernel used to compute the covariance matrix of the mean GP at corresponding inputs (K_0)
#' @param mean_mu posterior mean of the mean GP (mu_0). Needed to compute the log-likelihood
#' @param cov_mu posterior covariance matrix of the mean GP (mu_0). Needed to compute correction term
#' @param m_0 prior value of the mean parameter of the mean GP (mu_0). Length = 1 or nrow(db)
#' @param hp_0 hyperparameters for the kernel 0
#' @param hp_i hyperparameters for each kernel i
#' @param pen_diag value of the penalization of the diagonal
#'
#' @return value of expectation of joint log-likelihood of the model. The function to be maximised in step M
#' The full likelihood is composed of M+1 independent parts, depending on only theta_0, or theta_i respectively
#' for each i. The following code computes and sums these M+1 (modified) gaussian likelihoods.
#'
#' @examples
#' kern_i = kern_0 = kernel_sqrd_exp
#' hp_0 = tibble::tibble(sigma = 1, lengthscale = 0.5)
#' hp_i = list(hp_0,hp_0,hp_0,hp_0,hp_0)
#' db = mean_mu = tibble::tibble(ID=1:5,input= 1:5,Output= 2:6)
#'
#'
logL_monitoring = function(hp_0,hp_i, db, kern_i, kern_0, mean_mu, cov_mu, m_0, pen_diag)
{
  ## Mean GP (mu_0) is noiseless and thus has only 2 hp. We add a penalty on diag for numerical stability
  #pen_diag = sapply(hp_i, function(x) x$sigma) %>% mean

  ll_0 = logL_GP_mod(hp_0, db = mean_mu, mean = m_0, kern_0, cov_mu, pen_diag = pen_diag)

  funloop = function(i)
  {
    t_i = db %>% dplyr::filter(db$ID == i) %>% dplyr::pull(db$input)
    logL_GP_mod(hp_i[[i]], db %>% dplyr::filter(db$ID == i),
                mean = mean_mu %>% dplyr::filter(db$input %in% t_i) %>% dplyr::pull(db$Output),
                kern_i, cov_mu[paste0('X', t_i), paste0('X', t_i)]) %>% return()
  }
  sum_ll_i = sapply(unique(db$ID), funloop) %>% sum()

  return(-ll_0 - sum_ll_i)
}
