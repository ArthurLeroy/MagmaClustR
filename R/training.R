#' Training
#'
#' @param db Database with all individuals in training set. Column required : 'ID', input', 'Output'
#' @param prior_mean prior mean parameter of the mean GP (mu_0)
#' @param kern_0 kernel associated to covariance functions of the mean GP
#' @param kern_i kernel associated to common covariance functions of all individuals GPs
#' @param ini_hp_0 your initial hyperparameters for your kernel 0
#' @param ini_hp_i your initial hyperparameters for your individuals
#' @param common_hp your common Hyperparameters
#'
#' @return list of trained HP, boolean to indicate convergence
#' @export
#'
#' @examples
#' TRUE
training = function(db, prior_mean, ini_hp_0, ini_hp_i, kern_0, kern_i, common_hp = T)
{
  n_loop_max = 25
  db$ID = db$ID %>% as.character
  list_ID = unique(db$ID)

  ## TODO: Set m_0 to a coherent size := Ut_i

  cv = FALSE
  logLL_monitoring = - Inf
  list_plot = list()
  t1 = Sys.time()
  for(i in 1:n_loop_max)
  {
    #print(i)
    ## E-Step
    param = e_step(db, prior_mean, kern_0, kern_i, ini_hp_0,ini_hp_i)

    ## For visualising successive values of \mu_0
    #list_plot[[i]] = param$pred_GP %>% plot_gp(data_train = db)
    ## Return list_plot if you want monitoring graphs of the mean process' learning

    ## M-Step
    new_hp = m_step(db, ini_hp_0,ini_hp_i, mean = param$mean, cov = param$cov, kern_0, kern_i, prior_mean, common_hp)

    ## If something went wrong during the optimization
    if(new_hp %>% anyNA(recursive = T))
    {
      print(paste0('The M-step encountered an error at iteration : ', i))
      print('Training has stopped and the function returns values from the last valid iteration')
      break
    }

    ## Monitoring of the LL
    new_logLL_monitoring = logL_monitoring(ini_hp_0, ini_hp_i, db, kern_i, kern_0, param$mean, param$cov, prior_mean)
    # + 0.5 * log(det(param$cov)) for an exact likelihood but constant respectively to HPs
    paste0('logLL = ', new_logLL_monitoring) %>% print()
    diff_moni = new_logLL_monitoring - logLL_monitoring
    if(diff_moni < - 0.1){warning('Likelihood descreased')}

    logL_new_hp = logL_monitoring(new_hp, db, kern_i, kern_0, param$mean, param$cov, prior_mean)
    # + 0.5 * log(det(param$cov)) for an exact likelihood but constant respectively to HPs

    ## Testing the stoping condition
    eps = (logL_new_hp - new_logLL_monitoring) / abs(logL_new_hp)
    paste0('eps = ', eps) %>% print()
    if(eps < 1e-3)
    {
      if(eps > 0){hp = new_hp}
      cv = TRUE
      break
    }

    ## Update HP values and loglikelihood monitoring
    hp = new_hp
    logLL_monitoring = new_logLL_monitoring
  }
  t2 = Sys.time()
  list('hp' = hp, 'convergence' = cv, 'param' = param,
       'Time_train' =  difftime(t2, t1, units = "secs"),
       'plot' = list_plot) %>%
    return()
}

#' Training new gaussian process
#'
#' @param db Database with all individuals in training set
#' @param mean_mu mean value of mean GP at timestamps (obs + pred)
#' @param cov_mu covariance value of mean GP at timestamps (obs + pred)
#' @param ini_hp_i Initial values of each parameters of the HP to start the training.
#' @param kern_i Kernel associated to individual GPs.
#'
#' @return list of trained HP
#' @export
#'
#' @examples
#' TRUE
train_new_gp = function(db, mean_mu, cov_mu, ini_hp_i, kern_i)
{
  if(is.vector(mean_mu)){mean = mean_mu}
  else {mean = mean_mu %>% dplyr::filter(db$Timestamp %in% db$Timestamp) %>% dplyr::pull(db$Output) %>% as.vector}
  if(length(mean) == 1){mean = rep(mean, length(db$Timestamp))}

  if(is.matrix(cov_mu)){new_cov = cov_mu[paste0('X', db$Timestamp), paste0('X', db$Timestamp)]}
  else {new_cov = 0}

  new_hp = optimr::opm(ini_hp_i, fn = logL_GP, gr = gr_GP, db = db, mean = mean, kern = kern_i, new_cov = new_cov,
                       method = "L-BFGS-B", control = list(kkt = FALSE))

  ## If something went wrong during the optimization
  if(new_hp[1,] %>% anyNA())
  {
    print('Training has stopped and the function returns initial values of hyperparameters')
    new_hp = ini_hp_i
  }

  tibble::as_tibble(new_hp) %>% return()
}
