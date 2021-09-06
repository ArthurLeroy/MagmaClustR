#' Training Magma with a Variation of the EM algorithm
#'
#' The hyper-parameters and the hyper-posterior distribution involved in Magma
#' can be learned thanks to an EM algorithm implemented in \code{train_magma_VEM}.
#' By providing a dataset, the model hypotheses (hyper-prior mean parameter and
#' covariance kernels) and initialisation values for the hyper-parameters, the
#' function computes maximum likelihood estimates of the HPs as well as the
#' mean and covariance parameters of the Gaussian hyper-posterior distribution
#' of the mean process.
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
#'    "PERIO" and "RQ" are also available here)
#' @param ini_tau_i_k initial values of probabiliy to belong to each cluster for each individuals.
#' @param common_hp_k A boolean indicating whether hp are common among mean GPs (for each mu_k).
#' @param common_hp_i A boolean indicating whether hp are common among individual GPs (for each y_i).
#' @param prior_mean_k prior mean parameter of the K mean GPs (mu_k)
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#'
#' @return A list, containing the results of the EM algorithm used for training
#'    in MagmaClust. The elements of the list are:
#'    - hp_k: A tibble containing the trained hyper-parameters for the mean
#'    process' kernel.
#'    - hp_i: A tibble containing all the trained hyper-parameters for the
#'    individual processes' kernels.
#'    - pi_k :
#'
#'    - param :
#'
#'    - Converance: A logical value indicated whether the EM algorithm converged
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
  #browser()
  n_loop_max = 15
  list_ID = unique(data$ID)
  ID_k = names(prior_mean_k)
  hp_k = ini_hp_k # %>% list() %>% rep(length(ID_k))  %>% stats::setNames(nm = ID_k)
  hp_i = ini_hp_i # %>% list() %>% rep(length(list_ID))  %>% stats::setNames(nm = list_ID)

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
    new_logLL_monitoring = logL_monitoring_VEM(hp_k, hp_i, data, kern_i, kern_0, mu_k_param = param , m_k = prior_mean_k, pen_diag)
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
  list('hp_k' = hp_k, 'hp_i' = hp_i, 'pi_k' = pi_k,
       'convergence' = cv,  'param' = param,
       'Training_time' =  difftime(t2, t1, units = "secs")) %>%
    return()
}

#' Update tau star k EM
#'
#' @param db data
#' @param hp hp
#' @param pi_k pi
#' @param mean_k mean
#' @param cov_k cov
#' @param kern kernel
#'
#' @return update tau star
#' @export
#'
#' @examples
#' TRUE
update_tau_star_k_EM <- function(db, mean_k, cov_k, kern, hp, pi_k)
{
  #browser()
  names_k = names(mean_k)
  c_k = 0
  mat_logL = rep(NA, length(names_k))

  ## Extract the specific inputs
    input_i <- db %>%
      dplyr::pull(.data$Input) %>%
      round(digits = 10)

    #unique_input <- input_i %>% unique

  pi = unlist(pi_k)

  for(k in names_k)
  {
    c_k = c_k + 1

    round_mean <- mean_k[[k]] %>% round(digits = 10)

    mean = round_mean %>% dplyr::filter(.data$Input %in% input_i) %>% dplyr::pull(.data$Output)
    cov =  (kern_to_cov(db$Input, kern, hp) + cov_k[[k]][as.character(input_i), as.character(input_i)] )
    inv = tryCatch(solve(cov), error = function(e){MASS::ginv(cov)})

    mat_logL[c_k] =  dmnorm(db %>% dplyr::pull(.data$Output) , mean, inv, log = T) ## classic gaussian loglikelihood
  }
  ## We need to use the 'log-sum-exp' trick: exp(x - max(x)) / sum exp(x - max(x)) to remain numerically stable
  mat_L = exp(mat_logL - max(mat_logL))

  ((pi * mat_L)/ sum(pi * mat_L)) %>% as.vector %>% split(names_k) %>% return()

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
#' @param ini_hp_i A tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}, the individual processes' kernel.
#'    Required column : \code{ID}. The \code{ID} column contains the unique
#'    names/codes used to identify each individual/task. The other columns
#'    should be named according to the hyper-parameters that are used in
#'    \code{kern_i}.
#' @param kern_i A kernel function, associated with the individual GPs. ("SE",
#'    "PERIO" and "RQ" are also available here).
#' @param hp_i A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}.
#'
#' @return A list, containing the results of the EM algorithm used for training
#'    in MagmaClust. The elements of the list are:
#'    - theta_new :
#'    - tau_k :
#' @export
#'
#' @examples
#' k = seq_len(2)
#' m_k <- c("K1" = 0, "K2" = 0)
#'
#' db <- simu_db(N = 2, common_input = FALSE)
#' hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
#' hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))

#' ini_hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))
#' old_tau_i_k = MagmaClustR:::ini_tau_i_k(db = db, k = length(k), nstart = 50)
#'
#' training_test = train_magma_VEM(db, m_k, hp_k, ini_hp_i, "SE", "SE", old_tau_i_k, FALSE, FALSE, 0.1)
#'
#' timestamps = seq(0.01, 10, 0.01)
#' mu_k <- posterior_mu_k(db, timestamps, m_k, "SE", "SE", training_test)
#'
#'
#' list_hp <- train_magma_VEM(db, m_k, hp_k, hp_i, "SE", "SE", old_tau_i_k, FALSE, FALSE, 0.1)
#'
#' train_new_gp_EM(simu_db(M=1, covariate = FALSE), mu_k, ini_hp_i, "SE", hp_i = list_hp$hp_i)
train_new_gp_EM = function(data, param_mu_k, ini_hp_i, kern_i, hp_i = NULL)
{
  #browser()

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

      tau_k <- update_tau_star_k_EM(data, mean_mu_k, cov_mu_k, kern_i, hp, pi_k)



      ## M step
      LL_GP<- function(hp, data, kern_i)
      {
        inputs <- data %>% dplyr::pull(.data$Input) %>% unique %>% round(digits = 10)

        floop = function(k)
        {
          round_mean <- mean_mu_k[[k]] %>% round(digits = 10)

          mean = round_mean %>% dplyr::filter(.data$Input %in% inputs) %>% dplyr::pull(.data$Output)
          cov = (kern_to_cov(data$Input, kern_i, hp) +
                   cov_mu_k[[k]][as.character(inputs), as.character(inputs)])
          inv = tryCatch(solve(cov), error = function(e){MASS::ginv(cov)})

          (data$Input - tau_k[[k]] * dmnorm(data$Output, mean, inv, log = T)) %>%
            return()
        }
        sapply(names(mean_mu_k), floop) %>% sum() %>% return()
      }
      ## Extract the hyper-parameters associated
      par_i <- hp %>%
        dplyr::select(-.data$ID)

      new_hp = optimr::opm(
        par_i,
        LL_GP,
        db = data,
        kern = kern_i,
        method = "L-BFGS-B",
        control = list(kkt = FALSE)
        )[1,1:3]

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
    new_hp = hp_i %>% dplyr::slice(1)

    tau_k <- update_tau_star_k_EM(data, mean_mu_k, cov_mu_k, kern_i, new_hp, pi_k)

  }
  list('theta_new' = new_hp , 'tau_k' = tau_k) %>% return()
}


