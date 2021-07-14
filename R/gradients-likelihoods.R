#' Gradient Gaussian Process modif
#'
#' @param hp set of hyper-parameter
#' @param db full database with all individuals. Columns required : ID, input, Output
#' @param mean mean of your Gaussian Process
#' @param kern kernel used to compute the covariance matrix
#' @param new_cov posterior covariance matrix of the mean GP (mu_0).
#' @param pen_diag value of the penalization of the diagonal
#'
#' @return Gradient of the Gaussian Process modified
#' @export
#'
#' @examples
#' TRUE
gr_GP_mod = function(hp, db, mean, kern, new_cov, pen_diag = NULL)
{
  y = db$Output
  t = db$Input

  att = attributes(kern)
  deriv_hp1 = att$derivative_sigma()
  deriv_hp2 = att$derivative_lengthscale()

  ## Mean GP (mu_0) is noiseless and thus has only 2 hp. We add a penalty on diag for numerical stability
  sigma = ifelse((length(hp) == 3), hp$sigma, pen_diag)
  inv = kern_to_inv(t, kern, hp)
  prod_inv = inv %*% (y - mean)
  cste_term = prod_inv %*% t(prod_inv) + inv %*% ( new_cov %*% inv - diag(1, length(t)) )


  g_1 = 1/2 * (cste_term %*% kern_to_cov(t, deriv_hp1, hp)) %>%  diag() %>% sum()
  if(length(t) == 1)
  { ## Second hp has a 0 diagonal, and dist() return an error for only one observation
    g_2 = 0
  }
  else
  {
    g_2 = 1/2 * (cste_term %*% as.matrix(deriv_hp2(t,t, hp) ))  %>%  diag() %>% sum()
  }

  if(length(hp) == 3)
  {
    g_3 = hp$sigma * (cste_term %>% diag() %>% sum() )
    (- c(g_1, g_2, g_3)) %>% return()
  }
  else   (- c(g_1, g_2)) %>% return()
}

#' Gradient common Gaussian Process
#'
#' @param hp set of hyper-parameter
#' @param db full database with all individuals. Columns required : ID, input, Output
#' @param mean mean of the GP at union of observed input
#' @param kern kernel used to compute the covariance matrix at corresponding inputs
#' @param new_cov posterior covariance matrix of the mean GP (mu_0).
#'
#' @return Gradient of the common Gaussian Process
#' @export
#'
#' @examples
#' TRUE
gr_GP_mod_common_hp = function(hp, db, mean, kern, new_cov)
{
  g_1 = 0
  g_2 = 0
  g_3 = 0
  t_i_old = NULL
  hp_2 <- hp
  hp_2$sigma = 0

  att = attributes(kern)
  deriv_hp1 = att$derivative_sigma()
  deriv_hp2 = att$derivative_lengthscale()

  for(i in unique(db$ID))
  {
    t_i = db %>% dplyr::filter(db$ID == i) %>% dplyr::pull(db$input)
    input_i = paste0('X', t_i)
    y_i = db %>% dplyr::filter(db$ID == i) %>% dplyr::pull(db$output)

    if( !identical(t_i, t_i_old) )
    { ## We update the inverse cov matrix only if necessary (if different inputs)
      inv = kern_to_inv(t_i, kern, hp)
    }
    prod_inv = inv %*% (y_i - mean %>% dplyr::filter(db$input %in% t_i) %>% dplyr::pull(db$output))
    cste_term = prod_inv %*% t(prod_inv) + inv %*%
      ( new_cov[input_i,input_i] %*% inv - diag(1, length(t_i)) )

    g_1 = g_1 + 1/2 * (cste_term %*% kern_to_cov(t_i, deriv_hp1, hp_2)) %>%  diag() %>% sum()
    if(length(t_i) == 1)
    { ## Second hp has a 0 diagonal, and dist() return an error for only one observation
      g_2 = g_2 + 0
    }
    else
    {
      g_2 = g_2 + 1/2 * (cste_term %*% as.matrix(deriv_hp2(t_i,t_i, hp) ))  %>%  diag() %>% sum()
    }
    g_3 = g_3 + hp$sigma * (cste_term %>% diag() %>% sum() )
    t_i_old = t_i
  }
  return(- c(g_1, g_2, g_3))
}
