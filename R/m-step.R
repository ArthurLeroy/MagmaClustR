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
