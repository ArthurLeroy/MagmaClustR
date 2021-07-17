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

#' Log-Likelihood function of a Gaussian Process
#'
#' @param hp A tibble, data frame or named vector containing hyper-parameters.
#' @param db A tibble containing the values we want to compute the logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mean A vector, specifying the mean of the GP at the reference inputs.
#' @param kern A kernel function.
#' @param new_cov (optional) A matrix, corresponding to covariance parameter of
#'    the hyper-posterior. Used to compute the hyper-prior distribution of a new
#'    individual in Magma.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A number, corresponding to the value of Gaussian
#'    log-Likelihood (where the covariance can be the sum of the individual and
#'    the hyper-posterior's mean process covariances).
#' @examples
#' db <- tibble::tibble(Input = 1:5, Output = 2:6)
#' mean <- rep(0, 5)
#' hp <- tibble::tibble(variance = 1, lengthscale = 0.5)
#' new_cov = kern_to_cov(1:5, 'SE', hp)
#' MagmaClustR:::logL_GP(hp, db, mean, "SE", new_cov, 0.001)
logL_GP <- function(hp, db, mean, kern, new_cov, pen_diag) {
  ## Extract the input variables (reference Input + Covariates)
  input = db %>% dplyr::select(- .data$Output)

  ## Sum the two covariance matrices and inverse the result
  cov <- kern_to_cov(input, kern, hp) + new_cov
  diag <- diag(x = pen_diag, ncol = ncol(cov), nrow = nrow(cov))
  inv <- tryCatch((cov + diag) %>% chol() %>% chol2inv(),
                  error = function(e) {
                    MASS::ginv(cov + diag)
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
#' db <- tibble::tibble(Input = 1:5, Output = 2:6)
#' mean <- rep(0, 5)
#' hp <- tibble::tibble(variance = 1, lengthscale = 0.5)
#' new_cov = kern_to_cov(1:5, 'SE', hp)
#' MagmaClustR:::logL_GP_mod(hp, db, mean, "SE", new_cov, 0.001)
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
#' db <- simu_db(N = 10, common_input = TRUE)
#' mean <- tibble::tibble(Input = unique(db$Input), Output = 0)
#' hp <- tibble::tibble(variance = 1, lengthscale = 0.5)
#' new_cov <- kern_to_cov(unique(db$Input), 'SE', hp)
#' MagmaClustR:::logL_GP_mod_common_hp(hp, db, mean, "SE", new_cov, 0.001)
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

#' Log-Likelihood for monitoring the EM algorithm in Magma
#'
#' @param hp_0 A tibble or data frame, containing the hyper-parameters
#'    associated with the mean GP.
#' @param hp_i A tibble or data frame, containing the hyper-parameters with the
#'    individual GPs.
#' @param db A tibble or data frame. Columns required: ID, Input, Output.
#'    Additional columns for covariates can be specified.
#' @param m_0 A vector, corresponding to the prior mean of the mean GP.
#' @param kern_0 A kernel function, associated with the mean GP.
#' @param kern_i A kernel function, associated with the individual GPs.
#' @param mean A tibble, coming out of the E step, containing the Input and
#'    associated Output of the hyper-posterior mean parameter.
#' @param cov A matrix, coming out of the E step, being the hyper-posterior
#'    covariance parameter.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A number, expectation of joint log-likelihood of the model. This
#'    quantity is supposed to increase at each step of the EM algorithm, and
#'    thus used for monitoring the procedure.
#'
#' @examples
#' db <- simu_db(N = 10, common_input = TRUE)
#' m_0 <- rep(0, 10)
#' hp_0 <- tibble::tibble(variance = 1, lengthscale = 0.5)
#' hp_i <- MagmaClustR:::hp('SE', unique(db$ID))
#' mean <- tibble::tibble(Input = unique(db$Input), Output = 5)
#' cov <- kern_to_cov(unique(db$Input), 'SE', hp_0)
#' MagmaClustR:::logL_monitoring(hp_0, hp_i, db, m_0,
#'  "SE", "SE", mean, cov, 0.001)
logL_monitoring <- function(hp_0, hp_i, db, m_0, kern_0, kern_i,
                            mean, cov, pen_diag) {
  ## Compute the modified logL for the mean process
  ll_0 <- logL_GP_mod(hp_0, db = mean, mean = m_0, kern_0, cov, pen_diag)

  ## Sum over the individuals
  funloop <- function(i) {
    ## Extract the i-th specific reference inputs
    input_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Input)
    ## Extract the i-th specific hyper-parameters
    hp_i_i <- hp_i %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(- .data$ID)
    ## Extract the i-th specific Inputs and Output
    db_i = db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(- .data$ID)
    ## Extract the mean values associated with the i-th specific inputs
    mean_i = mean %>%
      dplyr::filter(.data$Input %in% input_i) %>%
      dplyr::pull(.data$Output)
    ## Extract the covariance values associated with the i-th specific inputs
    cov_i = cov[as.character(input_i), as.character(input_i)]

    ## Compute the modified logL for the individual processes
    logL_GP_mod(hp_i_i, db_i, mean_i, kern_i, cov_i, pen_diag) %>%
      return()
  }
  sum_ll_i <- sapply(unique(db$ID), funloop) %>% sum()
  ## Since the logL_GP_* functions return negative likelihoods for minimisation
  ## in the M-step, we need to x(-1) once more to retrieve the correct logL
  return(-ll_0 - sum_ll_i)
}
