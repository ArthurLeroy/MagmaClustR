#' Compute the Multivariate Gaussian likelihood
#'
#' Modification of the function \code{dmvnorm()} from the package
#' \code{mvtnorm}, providing an implementation of the Multivariate Gaussian
#' likelihood. This version uses inverse of the covariance function as argument
#' instead of the traditional covariance.
#'
#' @param x A vector, containing values the likelihood is evaluated on.
#' @param mu A vector or matrix, specifying the mean parameter.
#' @param inv_Sigma A matrix, specifying the inverse of covariance parameter.
#' @param log A logical value, indicating whether we return the log-likelihood.
#'
#' @return A number, corresponding to the Multivariate Gaussian log-likelihood.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
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
  } else {
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
#' @param post_cov (optional) A matrix, corresponding to covariance parameter of
#'    the hyper-posterior. Used to compute the hyper-prior distribution of a new
#'    individual in Magma.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A number, corresponding to the value of Gaussian
#'    log-Likelihood (where the covariance can be the sum of the individual and
#'    the hyper-posterior's mean process covariances).
#'
#' @keywords internal
#'
#' @examples
#' TRUE
logL_GP <- function(hp, db, mean, kern, post_cov, pen_diag) {
  ## Extract the input variables (reference Input + Covariates)
  input <- db %>% dplyr::select(-.data$Output)

  ## Sum the two covariance matrices and inverse the result
  cov <- kern_to_cov(input, kern, hp) + post_cov

  inv <- cov %>% chol_inv_jitter(pen_diag = pen_diag)

  (-dmnorm(db$Output, mean, inv, log = T)) %>%
    return()
}


#' Modified log-Likelihood function for GPs
#'
#' Log-Likelihood function involved in Magma during the maximisation step of
#' the training. The log-Likelihood is defined as a simple Gaussian likelihood
#' added with correction trace term.
#'
#' @param hp A tibble, data frame or named vector of hyper-parameters.
#' @param db A tibble containing values we want to compute logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mean A vector, specifying the mean of the GP at the reference inputs.
#' @param kern A kernel function.
#' @param post_cov A matrix, covariance parameter of the hyper-posterior.
#'    Used to compute the correction term.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A number, corresponding to the value of the modified Gaussian
#' log-Likelihood defined in Magma.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
logL_GP_mod <- function(hp, db, mean, kern, post_cov, pen_diag) {
  if (length(mean) == 1) {
    mean <- rep(mean, nrow(db))
  }
  ## mean is equal for all timestamps

  ## Extract the input variables (reference Input + Covariates)
  inputs <- db %>% dplyr::select(-.data$Output)
  ## Compute the inverse of the covariance matrix
  inv <- kern_to_inv(inputs, kern, hp, pen_diag)

  ## Classical Gaussian log-likelihood
  LL_norm <- -dmnorm(db$Output, mean, inv, log = T)
  ## Correction trace term (- 1/2 * Trace(inv %*% post_cov))
  cor_term <- 0.5 * sum(inv * post_cov)

  return(LL_norm + cor_term)
}

#' Modified log-Likelihood function with common HPs for GPs
#'
#' Log-Likelihood function involved in Magma during the maximisation step of
#' the training, in the particular case where the hyper-parameters are shared by
#' all individuals. The log-Likelihood is defined as a sum over all individuals
#' of Gaussian likelihoods added with correction trace terms.
#'
#' @param hp A tibble, data frame of hyper-parameters.
#' @param db A tibble containing the values we want to compute the logL on.
#'    Required columns: ID, Input, Output. Additional covariate columns are
#'    allowed.
#' @param mean A vector, specifying the mean of the GP at the reference inputs.
#' @param kern A kernel function.
#' @param post_cov A matrix, covariance parameter of the hyper-posterior.
#'    Used to compute the correction term.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A number, corresponding to the value of the modified Gaussian
#' log-Likelihood with common hyper-parameters defined in Magma.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
logL_GP_mod_common_hp <- function(hp, db, mean, kern, post_cov, pen_diag) {
  funloop <- function(i) {
    ## Extract the i-th specific reference inputs
    input_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Input)
    ## Extract the i-th specific Inputs and Output
    db_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(-.data$ID)
    ## Extract the mean values associated with the i-th specific inputs
    mean_i <- mean %>%
      dplyr::filter(.data$Input %in% input_i) %>%
      dplyr::pull(.data$Output)
    ## Extract the covariance values associated with the i-th specific inputs
    post_cov_i <- post_cov[as.character(input_i), as.character(input_i)]

    ## Compute the modified logL for the individual processes
    logL_GP_mod(hp, db_i, mean_i, kern, post_cov_i, pen_diag) %>%
      return()
  }
  sapply(unique(db$ID), funloop) %>%
    sum() %>%
    return()
}

#' Log-Likelihood for monitoring the EM algorithm in Magma
#'
#' @param hp_0 A named vector, tibble or data frame, containing the
#' hyper-parameters associated with the mean GP.
#' @param hp_i A tibble or data frame, containing the hyper-parameters with the
#'    individual GPs.
#' @param db A tibble or data frame. Columns required: ID, Input, Output.
#'    Additional columns for covariates can be specified.
#' @param m_0 A vector, corresponding to the prior mean of the mean GP.
#' @param kern_0 A kernel function, associated with the mean GP.
#' @param kern_i A kernel function, associated with the individual GPs.
#' @param post_mean A tibble, coming out of the E step, containing the Input and
#'    associated Output of the hyper-posterior mean parameter.
#' @param post_cov A matrix, coming out of the E step, being the hyper-posterior
#'    covariance parameter.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A number, expectation of joint log-likelihood of the model. This
#'    quantity is supposed to increase at each step of the EM algorithm, and
#'    thus used for monitoring the procedure.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
logL_monitoring <- function(hp_0,
                            hp_i,
                            db,
                            m_0,
                            kern_0,
                            kern_i,
                            post_mean,
                            post_cov,
                            pen_diag) {
  ## Compute the modified logL for the mean process
  ll_0 <- logL_GP_mod(hp = hp_0,
                      db = post_mean,
                      mean = m_0,
                      kern = kern_0,
                      post_cov = post_cov,
                      pen_diag = pen_diag)

  ## Sum over the individuals
  funloop <- function(i) {
    ## Extract the i-th specific reference inputs
    input_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Input)
    ## Extract the i-th specific hyper-parameters
    hp_i_i <- hp_i %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(-.data$ID)
    ## Extract the i-th specific Inputs and Output
    db_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(-.data$ID)
    ## Extract the mean values associated with the i-th specific inputs
    post_mean_i <- post_mean %>%
      dplyr::filter(.data$Input %in% input_i) %>%
      dplyr::pull(.data$Output)
    ## Extract the covariance values associated with the i-th specific inputs
    post_cov_i <- post_cov[as.character(input_i), as.character(input_i)]

    ## Compute the modified logL for the individual processes
    logL_GP_mod(hp_i_i, db_i, post_mean_i, kern_i, post_cov_i, pen_diag) %>%
      return()
  }
  sum_ll_i <- sapply(unique(db$ID), funloop) %>% sum()

  ## Compute the log-determinant term using Cholesky decomposition
  ## log(det(A)) = 2*sum(log(diag(chol(A))))
  det <- post_cov %>%
    chol() %>%
    diag() %>%
    log() %>%
    sum()

  ## Since the logL_GP_* functions return negative likelihoods for minimisation
  ## in the M-step, we need to x(-1) once more to retrieve the correct logL
  return(-ll_0 - sum_ll_i + det)
}

#' Compute a mixture of Gaussian log-likelihoods
#'
#' During the prediction step of MagmaClust, an EM algorithm is used to compute
#' the maximum likelihood estimator of the hyper-parameters along with
#' mixture probabilities for the new individual/task. This function implements
#' the quantity that is maximised (i.e. a sum of Gaussian log-likelihoods,
#' weighted by their mixture probabilities). It can also be used to monitor the
#' EM algorithm when providing the 'prop_mixture' argument, for proper
#' penalisation of the full log-likelihood.
#'
#' @param hp A tibble, data frame or named vector of hyper-parameters.
#' @param db A tibble containing data we want to evaluate the logL on.
#'    Required columns: Input, Output. Additional covariate columns are allowed.
#' @param mixture A tibble or data frame, indicating the mixture probabilities
#'    of each cluster for the new individual/task.
#' @param mean A list of hyper-posterior mean parameters for all clusters.
#' @param kern A kernel function.
#' @param post_cov A list of hyper-posterior covariance parameters for all
#'    clusters.
#' @param prop_mixture A tibble or a named vector. Each name of column or
#'    element should refer to a cluster. The value associated with each cluster
#'    is a number between 0 and 1, corresponding to the mixture
#'    proportions.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#'
#' @return A number, expectation of mixture of Gaussian log-likelihoods in
#'    the prediction step of MagmaClust. This quantity is supposed to increase
#'    at each step of the EM algorithm, and can be used for monitoring the
#'    procedure.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
sum_logL_GP_clust <- function(hp,
                              db,
                              mixture,
                              mean,
                              kern,
                              post_cov,
                              prop_mixture = NULL,
                              pen_diag) {
  ## Extract the observed (reference) Input
  input_obs <- db %>%
    dplyr::arrange(.data$Input) %>%
    dplyr::pull(.data$Input)

  ## Remove 'ID' if present in 'db'
  if ("ID" %in% names(db)) {
    db <- db %>% dplyr::select(-.data$ID)
  }

  ## Loop over the K clusters
  floop <- function(k) {
    tau_k <- mixture[[k]]
    mean_k <- mean[[k]] %>%
      dplyr::filter(.data$Input %in% input_obs) %>%
      dplyr::pull(.data$Output)

    cov_k <- post_cov[[k]][
      as.character(input_obs),
      as.character(input_obs)
    ]

    sum_LL <- (tau_k * logL_GP(hp, db, mean_k, kern, cov_k, pen_diag))

    ## If prop_mixture is provided, compute full likelihood for monitoring EM
    if (!is.null(prop_mixture)) {
      pi_k <- prop_mixture[[k]]
      ## To avoid numerical issues if evaluating log(0/0)
      log_frac <- ifelse((tau_k == 0 | (pi_k == 0)), 0, log(pi_k / tau_k))

      ## Return -sum_LL because its a quantity that is minimised otherwise
      sum_LL <- -sum_LL + tau_k * log_frac
    }

    return(sum_LL)
  }
  sapply(names(mean), floop) %>%
    sum() %>%
    return()
}
