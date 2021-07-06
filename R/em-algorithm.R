#' E Step
#'
#' @param db full database with all individuals. Columns required : ID, input, Output
#' @param m_0 mean of your Gaussian Process initial (mu_0)
#' @param kern_0 kernel used to compute the covariance matrix of the mean GP at corresponding inputs (K_0)
#' @param hp_0 hyperparameters for the kernel 0
#' @param hp_i hyperparameters for each kernel i
#' @param kern_i kernel used to compute the covariance matrix of individuals GP at corresponding inputs (Psi_i)
#' @param pen_diag value of the penalization of the diagonal
#'
#' @return mean and covariance parameters of the mean GP (mu_0)
#' @export
#'
#' @examples
e_step = function(db, m_0, kern_0, kern_i, hp_0,hp_i,pen_diag)
{
  all_t = unique(db$input) %>% sort()
  ## Mean GP (mu_0) is noiseless and thus has only 2 hp. We add a penalty on diag for numerical stability
  #pen_diag = sapply(hp$theta_i, function(x) x[[3]]) %>% mean

  #hp_0$sigma = pen_diag ## ou sigma + pen_diag
  inv_0 = kern_to_inv(all_t, kern_0, hp_0)
  inv_i = kern_to_inv(db, kern_i, hp_i)
  value_i = base::split(db$Output, list(db$ID))

  new_inv = update_inv(prior_inv = inv_0, list_inv_i = inv_i)
  new_cov = tryCatch(solve(new_inv), error = function(e){MASS::ginv(new_inv)}) ## fast or slow matrix inversion if singular
  #new_cov = solve(new_inv)

  weighted_mean = update_mean(prior_mean = m_0, prior_inv = inv_0, list_inv_i = inv_i, list_value_i = value_i)
  new_mean = new_cov %*% weighted_mean %>% as.vector()

  list('mean' = tibble::tibble('input' = all_t, 'Output' = new_mean), 'cov' = new_cov,
       'pred_GP' = tibble::tibble('input' = all_t, 'Mean' = new_mean, 'Var' = diag(new_cov)) ) %>% return()
}

#' Update Inv
#'
#' @param prior_inv inverse of the covariance matrix of the prior mean GP (mu_0). dim = all inputs
#' @param list_inv_i list of inverse of the covariance matrices of each individuals. dim = inputs of i
#'
#' @return inverse of the covariance of the posterior mean GP (mu_0 | (y_i)_i). dim = (all inputs)^2
#' @export
#'
#' @examples
update_inv = function(prior_inv, list_inv_i)
{
  new_inv = prior_inv

  for(x in list_inv_i)
  {
    inv_i = x
    common_times = intersect(row.names(inv_i), row.names(new_inv))
    new_inv[common_times, common_times] = new_inv[common_times, common_times] + inv_i[common_times, common_times]
  }
  return(new_inv)
}


#' update Mean
#'
#' @param prior_mean mean parameter of the prior mean GP (mu_0)
#' @param prior_inv inverse of the covariance matrix of the prior mean GP (mu_0). dim = (all inputs)^2
#' @param list_inv_i list of inverse of the covariance matrices of each individuals. dim = (inputs of i)^2
#' @param list_value_i list of outputs (y_i) for each individuals. dim = (inputs of i) x 1list of outputs (y_i) for each individuals. dim = (inputs of i) x 1
#'
#' @return mean parameter of the posterior mean GP (mu_0 | (y_i)_i). dim = (all inputs) x 1
#' @export
#'
#' @examples
update_mean = function(prior_mean, prior_inv, list_inv_i, list_value_i)
{
  if(length(prior_mean) == 1){prior_mean = rep(prior_mean, ncol(prior_inv))}
  weighted_mean = prior_inv %*% prior_mean
  #row.names(weithed_mean) = row.names(prior_inv)

  for(i in list_inv_i %>% names())
  {
    weighted_i = list_inv_i[[i]] %*% list_value_i[[i]]
    #row.names(weithed_i) = row.names(list_inv_i[[i]])

    common_times = intersect(row.names(weighted_i), row.names(weighted_mean))
    weighted_mean[common_times,] = weighted_mean[common_times,] + weighted_i[common_times,]
  }
  return(weighted_mean)
}


#' M-step of the training procedure for Magma
#'
#' @param db full database with all individuals. Columns required : ID, input, Output
#' @param old_hp the set of hyper-parameters from the previous step of the EM
#' @param mean mean parameter of the mean GP (mu_0), computed during the E step
#' @param cov covariance parameter of the mean GP (mu_0), computed during the E step
#' @param kern_0 kernel used to compute the covariance matrix of individuals GP at corresponding inputs (Psi_i)
#' @param kern_i kernel used to compute the covariance matrix of the mean GP at corresponding inputs (K_0)
#' @param m_0 prior value of the mean parameter of the mean GP (mu_0). Length = 1 or nrow(db)
#' @param common_hp your common Hyperparameters
#'
#' @return set of optimised hyper parameters for the different kernels of the model
#' @export
#'
#' @examples
m_step = function(db, old_hp, mean, cov, kern_0, kern_i, m_0, common_hp)
{
  list_ID = unique(db$ID)
  ## Mean GP (mu_0) is noiseless and thus has only 2 hp. We add a penalty on diag for numerical stability
  pen_diag = sapply(old_hp$theta_i, function(x) 2*x[[3]]) %>% mean

  t1 = Sys.time()
  new_theta_0 = optimr::opm(old_hp$theta_0, logL_GP_mod, gr = gr_GP_mod, db = mean, mean = m_0, kern = kern_0,
                            new_cov = cov, pen_diag = pen_diag, method = "L-BFGS-B", control = list(kkt = FALSE))[1,1:2]

  if(common_hp)
  {
    param = optimr::opm(old_hp$theta_i[[1]], logL_GP_mod_common_hp, gr = gr_GP_mod_common_hp , db = db, mean = mean,
                        kern = kern_i, new_cov = cov, method = "L-BFGS-B", control = list(kkt = F))[1,1:3]
    new_theta_i = param %>% list() %>% rep(length(list_ID))  %>% stats::setNames(nm = list_ID)
  }
  else
  {
    floop = function(i)
    {
      t_i = db %>% dplyr::filter(db$ID == i) %>% dplyr::pull(db$input)
      return(optimr::opm(old_hp$theta_i[[i]] %>% unlist(), logL_GP_mod, gr = gr_GP_mod , db = db %>% dplyr::filter(db$ID == i),
                         mean = mean %>% dplyr::filter(db$input %in% t_i) %>% dplyr::pull(db$Output), kern = kern_i,
                         new_cov = cov[paste0('X', t_i), paste0('X', t_i)], method = "L-BFGS-B",
                         control = list(kkt = F))[1,1:3])
    }
    new_theta_i = sapply(list_ID, floop, simplify = FALSE, USE.NAMES = TRUE)
  }

  t2 = Sys.time()
  print(t2-t1)
  list('theta_0' = new_theta_0, 'theta_i' = new_theta_i) %>% return()
}
