#' Training Magma with a Variation EM algorithm
#'
#' @param data A tibble or data frame. Columns required: \code{ID}, \code{Input}
#'    , \code{Output}.
#'    Additional columns for covariates can be specified.
#'    The \code{ID} column contains the unique names/codes used to identify each
#'    individual/task (or batch of data).
#'    The \code{Input} column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    \code{Output} column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at each
#'    reference \code{Input}.
#' @param ini_hp_k named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_0}, the mean process' kernel. The
#'    columns/elements should be named according to the hyper-parameters
#'    that are used in \code{kern_0}.
#' @param ini_hp_i A tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}, the individual processes' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each individual/task. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_i}.
#' @param kern_0 A kernel function, associated with the mean GP.
#' @param kern_i A kernel function, associated with the individual GPs.
#' @param ini_tau_i_k initial values of probabiliy to belong to each cluster for each individuals.
#' @param common_hp_k boolean indicating whether hp are common among mean GPs (for each mu_k)
#' @param common_hp_i boolean indicating whether hp are common among individual GPs (for each y_i)
#' @param prior_mean_k prior mean parameter of the K mean GPs (mu_k)
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#'
#' @return A list, containing the results of the EM algorithm used for training
#'    in Magma. The elements of the list are:
#'    - hp_k: A tibble containing the trained hyper-parameters for the mean
#'    process' kernel.
#'    - hp_i: A tibble containing all the trained hyper-parameters for the
#'    individual processes' kernels.
#'    - pi_k :
#'
#'    - param :
#'
#'    - Converged: A logical value indicated whether the EM algorithm converged
#'    or not.
#'    - Training_time: Total running time of the complete training.
#' @export
#'
#' @examples
#' k = seq_len(3)
#' m_k <- c("K1" = 0, "K2" = 0, "K3" = 0)
#'
#' db <- simu_db(N = 10, common_input = TRUE)
#' hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
#' hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))
#' old_tau_i_k = MagmaClustR:::ini_tau_i_k(db = db, k = length(k), nstart = 50)
#'
#' train_magma_VEM(db, m_k, hp_k, hp_i, "SE", "SE", old_tau_i_k, FALSE, FALSE, 0.1)

train_magma_VEM = function(data, prior_mean_k, ini_hp_k, ini_hp_i,
                        kern_0, kern_i, ini_tau_i_k = NULL,
                        common_hp_k = F, common_hp_i = F, pen_diag)
{
  browser()
  n_loop_max = 15
  list_ID = unique(data$ID)
  ID_k = names(prior_mean_k)
  hp_k = ini_hp_k %>% list() %>% rep(length(ID_k))  %>% stats::setNames(nm = ID_k)
  hp_i = ini_hp_i %>% list() %>% rep(length(list_ID))  %>% stats::setNames(nm = list_ID)

  cv = 'FALSE'
  if(is.null(ini_tau_i_k)){ini_tau_i_k = ini_tau_i_k(data, k = length(ID_k), nstart = 50)}
  tau_i_k = ini_tau_i_k
  hp_k[['pi']] = sapply( tau_i_k, function(x) x %>% unlist() %>% mean() )
  logLL_monitoring = - Inf
  t1 = Sys.time()

  for(i in 1:n_loop_max)
  {
    print(i)
    ## E-Step
    param = e_step_VEM(data, prior_mean_k, kern_0, kern_i, hp_k, hp_i, tau_i_k, pen_diag)

    ## Monitoring of the LL
    new_logLL_monitoring = logL_monitoring_VEM(hp_k, hp_i, data, kern_i, kern_0, mu_k_param = param , m_k = prior_mean_k)
    #0.5 * (length(param$cov) * nrow(param$cov[[1]]) +
    #Reduce('+', lapply(param$cov, function(x) log(det(x)))) )
    print(new_logLL_monitoring)
    diff_moni = new_logLL_monitoring - logLL_monitoring

    if(diff_moni < - 0.1){warning('Likelihood descreased')}

    ## M-Step
    new_hp = m_step_VEM(data, hp_k, hp_i, list_mu_param = param, kern_0, kern_i, prior_mean_k, common_hp_k, common_hp_i, pen_diag)

    if(new_hp %>% anyNA(recursive = T))
    {
      print(paste0('The M-step encountered an error at iteration : ', i))
      print('Training has stopped and the function returns values from the last valid iteration')
      break
    }

    ## Testing the stoping condition
    logL_new = logL_monitoring_VEM(new_hp$hp_k, new_hp$hp_i, data, kern_i, kern_0, mu_k_param = param, m_k = prior_mean_k, pen_diag)
    eps = (logL_new - logL_monitoring_VEM(hp_k, hp_i, data, kern_i, kern_0, mu_k_param = param, m_k = prior_mean_k, pen_diag)) /
      abs(logL_new)

    print(c('eps', eps))
    if(eps < 1e-1)
    {
      if(eps > 0){hp = new_hp}
      cv = TRUE
      break
    }

    hp_i = new_hp$hp_i
    hp_k = new_hp$hp_k
    pi_k = new_hp$pi_k

    tau_i_k = param$tau_i_k
    logLL_monitoring = new_logLL_monitoring
  }
  t2 = Sys.time()
  list('hp_0' = hp_k, 'hp_i' = hp_i, 'pi_k' = pi_k,
       'convergence' = cv,  'param' = param,
       'Training_time' =  difftime(t2, t1, units = "secs")) %>%
    return()
}

#' Training new GP with the EM algorithm
#'
#' @param data A tibble or data frame. Columns required: \code{ID}, \code{Input}
#'    , \code{Output}.
#'    Additional columns for covariates can be specified.
#'    The \code{ID} column contains the unique names/codes used to identify each
#'    individual/task (or batch of data).
#'    The \code{Input} column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    \code{Output} column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at each
#'    reference \code{Input}.
#' @param param_mu_k list of parameters for the K mean Gaussian processes
#' @param ini_hp_i your initial hyperparameters for your individuals
#' @param kern_i A kernel function, associated with the individual GPs.
#' @param hp_i A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}.
#'
#' @return list of trained HP, tau_k
#' @export
#'
#' @examples
train_new_gp_EM = function(data, param_mu_k, ini_hp_i, kern_i, hp_i = NULL)
{
  mean_mu_k = param_mu_k$mean
  cov_mu_k = param_mu_k$cov
  pi_k = lapply(param_mu_k$tau_i_k, function(x) Reduce("+", x)/ length(x))
  if(is.null(hp_i))
  {
    n_loop_max = 25
    hp = ini_hp_i

    for(i in 1:n_loop_max)
    {
      ## E step
      names_k = names(mean_mu_k)
      c_k = 0
      mat_logL = rep(NA, length(names_k))
      t_i = data %>% dplyr::pull(.data$Input)
      pi = unlist(pi_k)

      for(k in names_k)
      {
        c_k = c_k + 1

        mean = mean_mu_k[[k]] %>% dplyr::filter(.data$Input %in% t_i) %>% dplyr::pull(.data$Output)
        cov =  (kern_to_cov(data$Input, kern_i, hp) + cov_mu_k[[k]][paste0('X',t_i), paste0('X',t_i)] )
        inv = tryCatch(solve(cov), error = function(e){MASS::ginv(cov)})

        mat_logL[c_k] =  dmnorm(data %>% dplyr::pull(.data$Output) , mean, inv, log = T) ## classic gaussian loglikelihood
      }
      ## We need to use the 'log-sum-exp' trick: exp(x - max(x)) / sum exp(x - max(x)) to remain numerically stable
      mat_L = exp(mat_logL - max(mat_logL))

      tau_k = ( ((pi * mat_L)/ sum(pi * mat_L)) %>% as.vector %>% split(names_k) )



      ## M step
      LL_GP<- function(hp, data, kern_i)
      {
        floop = function(k)
        {
          mean = mean_mu_k[[k]] %>% dplyr::filter(.data$Input %in% t) %>% dplyr::pull(.data$Output)
          cov = (kern_to_cov(data$Input, kern_i, hp) +
                   cov_mu_k[[k]][paste0('X',t), paste0('X',t)])
          inv = tryCatch(solve(cov), error = function(e){MASS::ginv(cov)})

          (data$Input - tau_k[[k]] * dmnorm(data$Output, mean, inv, log = T)) %>%
            return()
        }
        sapply(names(mean_mu_k), floop) %>% sum() %>% return()
      }
      new_hp = optimr::opm(hp, LL_GP, db = data, kern = kern_i, method = "L-BFGS-B", control = list(kkt = FALSE))[1,1:3]
      if(new_hp %>% anyNA(recursive = T))
      {
        print(paste0('The M-step encountered an error at iteration : ', i))
        print('Training has stopped and the function returns values from the last valid iteration')
        break
      }

      ## Testing the stoping condition
      eps = (new_hp - hp) %>% abs %>% sum
      print(c('tau_k', tau_k %>% unlist))
      print(c('eps', eps))
      if(eps>0 & eps < 1e-3)
      {
        cv = 'TRUE'
        break
      }
      hp = new_hp
    }
  }
  else
  {
    new_hp = hp_i$hp_i[[1]]

    ## E step
    names_k = names(mean_mu_k)
    c_k = 0
    mat_logL = rep(NA, length(names_k))
    t_i = data %>% dplyr::pull(.data$Input)
    pi = unlist(pi_k)

    for(k in names_k)
    {
      c_k = c_k + 1

      mean = mean_mu_k[[k]] %>% dplyr::filter(.data$Input %in% t_i) %>% dplyr::pull(.data$Output)
      cov =  (kern_to_cov(data$Input, kern_i, new_hp) + cov_mu_k[[k]][paste0('X',t_i), paste0('X',t_i)] )
      inv = tryCatch(solve(cov), error = function(e){MASS::ginv(cov)})

      mat_logL[c_k] =  dmnorm(data %>% dplyr::pull(.data$Output) , mean, inv, log = T) ## classic gaussian loglikelihood
    }
    ## We need to use the 'log-sum-exp' trick: exp(x - max(x)) / sum exp(x - max(x)) to remain numerically stable
    mat_L = exp(mat_logL - max(mat_logL))

    tau_k = ( ((pi * mat_L)/ sum(pi * mat_L)) %>% as.vector %>% split(names_k) )


  }
  list('theta_new' = new_hp , 'tau_k' = tau_k) %>% return()
}


