#' Compute the Multivariate Gaussian likelihood
#'
#' Modification of the function \code{dmvnorm()} from the package \code{mvtnorm},
#' providing an implementation of the Multivariate Gaussian likelihood. This
#' version uses inverse of the covariance function as argument instead of the
#' traditional covariance.
#'
#' @param x A vector, containing values the likelihood is evaluated on.
#' @param mu A vector or matrix, specifying the mean parameter.
#' @param inv_Sigma A matrix, specifying the inverse of covariance parameter.
#' @param log A logical value, indicating whether we return the log-likelihood.
#'
#' @return A number, corresponding to the Multivariate Gaussian (log)-likelihood.
#'
#' @examples
#' MagmaClustR:::dmnorm(c(1, 2), c(0, 0), cbind(c(1, 0), c(0, 1)), TRUE)
dmnorm <- function(x, mu, inv_Sigma, log = FALSE) {
  if (is.vector(x)) {
    x <- t(as.matrix(x))
  }
  n <- length(mu)
  if (is.vector(mu)) {
    p <- length(mu)
    if (is.matrix(x)) {
      mu <- matrix(rep(mu, nrow(x)), ncol = p, byrow = TRUE)
    }
  }
  else {
    p <- ncol(mu)
  }
  if (!all(dim(inv_Sigma) == c(p, p)) || nrow(x) != nrow(mu)) {
    stop("incompatible arguments")
  }

  z <- t(x - mu)
  logdetS <- try(-determinant(inv_Sigma, logarithm = TRUE)$modulus,
    silent = TRUE
  )
  attributes(logdetS) <- NULL

  ssq <- t(z) %*% inv_Sigma %*% z
  loglik <- (-(n * (log(2 * pi)) + logdetS + ssq) / 2) %>% as.vector()
  if (log) {
    return(loglik)
  } else {
    return(exp(loglik))
  }
}

#' Log Likelihood function of the Gaussian Process
#'
#' @param hp The hyperparameters of your kernel
#' @param db tibble containing the values we want to compute. Required columns : input, Output
#' @param mean mean of the GP at corresponding inputs
#' @param kern kernel used to compute the covariance matrix at corresponding inputs
#' @param new_cov posterior covariance matrix of the mean GP (mu_0). Used to compute correction term (cor_term)
#'
#' @return value of the Gaussian log-likelihood for one GP as it appears in the model
#'
#' @examples
#' logL_GP(
#'   tibble::tibble(variance = 1, lengthscale = 0.5),
#'   tibble::tibble(Input = 1:5, Output = 2:6), rep(3, 5), se_kernel, 0
#' )
logL_GP <- function(hp, db, mean, kern, new_cov) {
  ## Extract the input variables (reference Input + Covariates)
  input = db %>% dplyr::select(- .data$Output)

  ## Sum the two covariance matrices and inverse the result
  cov <- kern_to_cov(input, kern, hp) + new_cov
  inv <- tryCatch(cov %>% chol() %>% chol2inv(),
    error = function(e) {
      MASS::ginv(cov)
    }
  )
  (-dmnorm(db$Output, mean, inv, log = T)) %>%
    return()
}


#' Modified log-Likelihood function for GPs
#'
#' Log-Likelihood function involved in Magma during the maximisation step of
#' the training. The log-Likelihood is defined as a simple Gaussian likelihood
#' added with correction trace term.
#'
#' @param hp A tibble, data frame or name vector of hyper-parameters.
#' @param db A tibble containing values we want to compute logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mean A vector, specifying the mean of the GP at the reference inputs.
#' @param kern A kernel function.
#' @param new_cov A matrix, covariance parameter of the hyper-posterior.
#'    Used to compute the correction term.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A number, corresponding to the value of the modified Gaussian
#' log-Likelihood defined in Magma.
#'
#' @examples
#' db <- tibble(Input = 1:5, Output = 2:6)
#' mean <- rep(0, 5)
#' hp <- tibble(variance = 1, lengthscale = 0.5)
#' logL_GP_mod(hp, db, mean, "SE", 0, 0.001)
logL_GP_mod <- function(hp, db, mean, kern, new_cov, pen_diag) {
  ## Extract the input variables (reference Input + Covariates)
  inputs <- db %>% dplyr::select(-.data$Output)
  ## Compute the inverse of the covariance matrix
  inv <- kern_to_inv(inputs, kern, hp, pen_diag)

  ## Classical Gaussian log-likelihood
  LL_norm <- - dmnorm(db$Output, mean, inv, log = T)
  ## Correction trace term (- 1/2 * Trace(inv %*% new_cov))
  cor_term <- 0.5 * (inv * new_cov) %>% sum()

  return(LL_norm + cor_term)
}

#' Modified log-Likelihood function with common HPs for GPs
#'
#' Log-Likelihood function involved in Magma during the maximisation step of
#' the training, in the particular case where the hyper-parameters are shared by
#' all individuals. The log-Likelihood is defined as a sum over all individuals
#' of Gaussian likelihoods added with correction trace terms.
#'
#' @param hp A tibble, data frame or name vector of hyper-parameters.
#' @param db A tibble containing the values we want to compute the logL on.
#'    Required columns: ID, Input, Output. Additional covariate columns are
#'    allowed.
#' @param mean A vector, specifying the mean of the GP at the reference inputs.
#' @param kern A kernel function.
#' @param new_cov A matrix, covariance parameter of the hyper-posterior.
#'    Used to compute the correction term.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A number, corresponding to the value of the modified Gaussian
#' log-Likelihood with common hyper-parameters defined in Magma.
#'
#' @examples
#' db <- simu_db(N = 10, common_input = T)
#' mean <- tibble::tibble(Input = unique(db$Input), Output = 0)
#' hp <- tibble::tibble(variance = 1, lengthscale = 0.5)
#' cov <- kern_to_cov(unique(db$Input), 'SE', hp)
#' logL_GP_mod_common_hp(hp, db, mean, "SE", cov, 0.001)
logL_GP_mod_common_hp <- function(hp, db, mean, kern, new_cov, pen_diag) {

  ## Initialise the value of sums and inputs
  LL_norm = 0
  cor_term = 0
  input_i_old <- NULL
  ## Loop over individuals to compute the sum of log-Likelihoods
  for (i in unique(db$ID)){
    ## Extract the i-th specific reference inputs
    input_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Input)
    ## Extract the Output values associated with the i-th specific inputs
    output_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Output)
    ## Extract the mean values associated with the i-th specific inputs
    mean_i = mean %>%
      dplyr::filter(.data$Input %in% input_i) %>%
      dplyr::pull(.data$Output)

    ## Update the inverse cov matrix only if necessary (if different inputs)
    if (!identical(input_i, input_i_old)) {
      ## Extract the i-th specific inputs (reference + covariates)
      inputs_i = db %>%
        dplyr::filter(.data$ID == i) %>%
        dplyr::select(- c(.data$ID, .data$Output))
      inv <- kern_to_inv(inputs_i, kern, hp, pen_diag)
    }

    ## Add a new term to the sum of individual log-Likelihoods
    LL_norm <- LL_norm - dmnorm(output_i, mean_i, inv, log = T)
    ## Correction trace term
    cor_term <- cor_term +
      0.5 * sum(inv * new_cov[as.character(input_i), as.character(input_i)])
    ## Keep track of the reference inputs
    input_i_old <- input_i
  }

  return(LL_norm + cor_term)
}

#' Log-Likelihood for monitoring the EM algorithm
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
#' kern_i <- kern_0 <- kernel_sqrd_exp
#' hp_0 <- tibble::tibble(sigma = 1, lengthscale = 0.5)
#' hp_i <- list(hp_0, hp_0, hp_0, hp_0, hp_0)
#' db <- mean_mu <- tibble::tibble(ID = 1:5, input = 1:5, Output = 2:6)
logL_monitoring <- function(hp_0, hp_i, db, kern_i, kern_0, mean_mu, cov_mu, m_0, pen_diag) {
  ## Mean GP (mu_0) is noiseless and thus has only 2 hp. We add a penalty on diag for numerical stability
  # pen_diag = sapply(hp_i, function(x) x$sigma) %>% mean

  ll_0 <- logL_GP_mod(hp_0, db = mean_mu, mean = m_0, kern_0, cov_mu, pen_diag = pen_diag)

  funloop <- function(i) {
    t_i <- db %>%
      dplyr::filter(db$ID == i) %>%
      dplyr::pull(db$input)
    logL_GP_mod(hp_i[[i]], db %>% dplyr::filter(db$ID == i),
      mean = mean_mu %>% dplyr::filter(db$input %in% t_i) %>% dplyr::pull(db$Output),
      kern_i, cov_mu[paste0("X", t_i), paste0("X", t_i)]
    ) %>% return()
  }
  sum_ll_i <- sapply(unique(db$ID), funloop) %>% sum()

  return(-ll_0 - sum_ll_i)
}
