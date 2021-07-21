#' Title
#'
#' @param db
#' @param prior_mean_k
#' @param ini_hp
#' @param kern_0
#' @param kern_i
#' @param ini_tau_i_k
#' @param common_hp_k
#' @param common_hp_i
#'
#' @return
#' @export
#'
#' @examples
training_VEM = function(db, prior_mean_k, ini_hp = list('theta_k' = c(1, 1, 0.2), 'theta_i' = c(1, 1, 0.2)),
                        kern_0 = kernel_mu, kern_i = kernel, ini_tau_i_k = NULL, common_hp_k = T, common_hp_i = T)
{ ## db : database with all individuals in training set. Column required : 'ID', Timestamp',  'Output'
  ## prior_mean : prior mean parameter of the K mean GPs (mu_k)
  ## ini_hp : initial values of HP for the kernels
  ## kern_0 : kernel associated to covariance functions of the mean GP
  ## kern_i : kernel associated to common covariance functions of all individuals GPs
  ## ini_tau_i_k : initial values of probabiliy to belong to each cluster for each individuals.
  ####
  ## return : list of trained HP, boolean to indicate convergence
  #browser()
  n_loop_max = 15
  list_ID = unique(db$ID)
  ID_k = names(prior_mean_k)
  hp = list('theta_k' = ini_hp$theta_k %>% list() %>% rep(length(ID_k))  %>% setNames(nm = ID_k),
            'theta_i' = ini_hp$theta_i %>% list() %>% rep(length(list_ID))  %>% setNames(nm = list_ID))
  cv = 'FALSE'
  if(is.null(ini_tau_i_k)){ini_tau_i_k = ini_tau_i_k(db, k = length(ID_k), nstart = 50)}
  tau_i_k = ini_tau_i_k
  hp[['pi_k']] = sapply( tau_i_k, function(x) x %>% unlist() %>% mean() )
  logLL_monitoring = - Inf
  t1 = Sys.time()

  for(i in 1:n_loop_max)
  {
    print(i)
    ## E-Step
    param = e_step_VEM(db, prior_mean_k, kern_0, kern_i, hp, tau_i_k)

    ## Monitoring of the LL
    new_logLL_monitoring = logL_monitoring_VEM(hp, db, kern_i, kern_0, mu_k_param = param , m_k = prior_mean_k)
    #0.5 * (length(param$cov) * nrow(param$cov[[1]]) +
    #Reduce('+', lapply(param$cov, function(x) log(det(x)))) )
    print(new_logLL_monitoring)
    diff_moni = new_logLL_monitoring - logLL_monitoring

    if(diff_moni < - 0.1){warning('Likelihood descreased')}

    ## M-Step
    new_hp = m_step_VEM(db, hp, list_mu_param = param, kern_0, kern_i, prior_mean_k, common_hp_k, common_hp_i)

    if(new_hp %>% anyNA(recursive = T))
    {
      print(paste0('The M-step encountered an error at iteration : ', i))
      print('Training has stopped and the function returns values from the last valid iteration')
      break
    }

    ## Testing the stoping condition
    logL_new = logL_monitoring_VEM(new_hp, db, kern_i, kern_0, mu_k_param = param, m_k = prior_mean_k)
    eps = (logL_new - logL_monitoring_VEM(hp, db, kern_i, kern_0, mu_k_param = param, m_k = prior_mean_k)) /
      abs(logL_new)

    print(c('eps', eps))
    if(eps < 1e-1)
    {
      if(eps > 0){hp = new_hp}
      cv = TRUE
      break
    }
    hp = new_hp
    tau_i_k = param$tau_i_k
    logLL_monitoring = new_logLL_monitoring
  }
  t2 = Sys.time()
  list('hp' = new_hp, 'convergence' = cv,  'param' = param,
       'Time_train' =  difftime(t2, t1, units = "secs")) %>%
    return()
}
