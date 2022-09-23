#' E-Step of the EM algorithm
#'
#' Expectation step of the EM algorithm to compute the parameters of the
#' hyper-posterior Gaussian distribution of the mean process in Magma.
#'
#' @param db A tibble or data frame. Columns required: ID, Input, Output.
#'    Additional columns for covariates can be specified.
#' @param m_0 A vector, corresponding to the prior mean of the mean GP.
#' @param kern_0 A kernel function, associated with the mean GP.
#' @param kern_i A kernel function, associated with the individual GPs.
#' @param hp_0 A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_0}.
#' @param hp_i A tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#' @return A named list, containing the elements \code{mean}, a tibble
#' containing the Input and associated Output of the hyper-posterior's mean
#' parameter, and \code{cov}, the hyper-posterior's covariance matrix.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
e_step <- function(db,
                   m_0,
                   kern_0,
                   kern_i,
                   hp_0,
                   hp_i,
                   pen_diag) {
  ## Extract the union of all reference inputs provided in the training data
  all_inputs <- db %>%
    dplyr::select(-.data$ID, -.data$Output) %>%
    unique() %>%
    dplyr::arrange(.data$Reference)

  ## Compute all the inverse covariance matrices
  inv_0 <- kern_to_inv(all_inputs, kern_0, hp_0, pen_diag)
  list_inv_i <- list_kern_to_inv(db, kern_i, hp_i, pen_diag)
  ## Create a named list of Output values for all individuals
  list_output_i <- base::split(db$Output, list(db$ID))

  ## Update the posterior inverse covariance ##
  post_inv <- inv_0
  for (inv_i in list_inv_i)
  {
    ## Collect the input's common indices between mean and individual processes
    co_input <- intersect(row.names(inv_i), row.names(post_inv))
    ## Sum the common inverse covariance's terms
    post_inv[co_input, co_input] <- post_inv[co_input, co_input] +
      inv_i[co_input, co_input]
  }

  post_cov <- post_inv %>%
    chol_inv_jitter(pen_diag = pen_diag) %>%
    `rownames<-`(all_inputs %>%
                   dplyr::pull(.data$Reference)
    ) %>%
    `colnames<-`(all_inputs %>%
                   dplyr::pull(.data$Reference)
    )
  ##############################################

  ## Update the posterior mean ##
  weighted_0 <- inv_0 %*% m_0
  for (i in names(list_inv_i))
  {
    ## Compute the weighted mean for the i-th individual
    weighted_i <- list_inv_i[[i]] %*% list_output_i[[i]]
    ## Collect the input's common indices between mean and individual processes
    co_input <- intersect(row.names(weighted_i), row.names(weighted_0))
    ## Sum the common weighted mean's terms
    weighted_0[co_input, ] <- weighted_0[co_input, ] +
      weighted_i[co_input, ]
  }
  ## Compute the updated mean parameter
  post_mean <- post_cov %*% weighted_0 %>% as.vector()
  ##############################################

  ## Format the mean parameter of the hyper-posterior distribution
  tib_mean <- tibble::tibble(all_inputs,
                             "Output" = post_mean
  )
  list(
    "mean" = tib_mean,
    "cov" = post_cov
  ) %>%
    return()
}


#' M-Step of the EM algorithm
#'
#' Maximisation step of the EM algorithm to compute hyper-parameters of all the
#' kernels involved in Magma.
#'
#' @param db A tibble or data frame. Columns required: ID, Input, Output.
#'    Additional columns for covariates can be specified.
#' @param m_0 A vector, corresponding to the prior mean of the mean GP.
#' @param kern_0 A kernel function, associated with the mean GP.
#' @param kern_i A kernel function, associated with the individual GPs.
#' @param old_hp_0 A named vector, tibble or data frame, containing the
#'    hyper-parameters from the previous M-step (or initialisation) associated
#'     with the mean GP.
#' @param old_hp_i A tibble or data frame, containing the hyper-parameters
#'    from the previous M-step (or initialisation) associated with the
#'    individual GPs.
#' @param post_mean A tibble, coming out of the E step, containing the Input and
#'    associated Output of the hyper-posterior mean parameter.
#' @param post_cov A matrix, coming out of the E step, being the hyper-posterior
#'    covariance parameter.
#' @param common_hp A logical value, indicating whether the set of
#'    hyper-parameters is assumed to be common to all indiviuals.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#' numerical issues when inverting nearly singular matrices.
#'
#' @return A named list, containing the elements \code{hp_0}, a tibble
#'    containing the hyper-parameters associated with the mean GP,
#'    \code{hp_i}, a tibble containing the hyper-parameters
#'    associated with the individual GPs.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
m_step <- function(db,
                   m_0,
                   kern_0,
                   kern_i,
                   old_hp_0,
                   old_hp_i,
                   post_mean,
                   post_cov,
                   common_hp,
                   pen_diag) {

  list_ID <- unique(db$ID)
  list_hp_0 <- old_hp_0 %>% names()
  list_hp_i <- old_hp_i %>%
    dplyr::select(-.data$ID) %>%
    names()

  ## Detect whether the kernel_0 provides derivatives for its hyper-parameters
  if (kern_0 %>% is.function()) {
    if (!("deriv" %in% methods::formalArgs(kern_0))) {
      gr_GP_mod <- NULL
    }
  }

  ## Detect whether the kernel_i provides derivatives for its hyper-parameters
  if (kern_i %>% is.function()) {
    if (!("deriv" %in% methods::formalArgs(kern_i))) {
      gr_GP_mod <- NULL
      gr_GP_mod_common_hp <- NULL
    }
  }

  ## Optimise hyper-parameters of the mean process
  new_hp_0 <- stats::optim(
    par = old_hp_0,
    fn = logL_GP_mod,
    gr = gr_GP_mod,
    db = post_mean,
    mean = m_0,
    kern = kern_0,
    post_cov = post_cov,
    pen_diag = pen_diag,
    method = "L-BFGS-B",
    control = list(factr = 1e13, maxit = 25)
  )$par %>%
    tibble::as_tibble_row()

  ## Check whether hyper-parameters are common to all individuals
  if (common_hp) {
    ## Retrieve the adequate initial hyper-parameters
    par_i <- old_hp_i %>%
      dplyr::select(-.data$ID) %>%
      dplyr::slice(1)
    ## Optimise hyper-parameters of the individual processes
    new_hp_i <- stats::optim(
      par = par_i,
      fn = logL_GP_mod_common_hp,
      gr = gr_GP_mod_common_hp,
      db = db,
      mean = post_mean,
      kern = kern_i,
      post_cov = post_cov,
      pen_diag = pen_diag,
      method = "L-BFGS-B",
      control = list(factr = 1e13, maxit = 25)
    )$par %>%
      tibble::as_tibble_row() %>%
      tidyr::uncount(weights = length(list_ID)) %>%
      dplyr::mutate("ID" = list_ID, .before = 1)
  } else {
    floop <- function(i) {
      ## Extract the i-th specific inputs
      input_i <- db %>%
        dplyr::filter(.data$ID == i) %>%
        dplyr::pull(.data$Reference)
      ## Extract the mean values associated with the i-th specific inputs
      post_mean_i <- post_mean %>%
        dplyr::filter(.data$Reference %in% input_i) %>%
        dplyr::pull(.data$Output)
      ## Extract the covariance values associated with the i-th specific inputs
      post_cov_i <- post_cov[as.character(input_i), as.character(input_i)]
      ## Extract the data associated with the i-th individual
      db_i <- db %>%
        dplyr::filter(.data$ID == i) %>%
        dplyr::select(-.data$ID)
      ## Extract the hyper-parameters associated with the i-th individual
      par_i <- old_hp_i %>%
        dplyr::filter(.data$ID == i) %>%
        dplyr::select(-.data$ID)

      ## Optimise hyper-parameters of the individual processes
      stats::optim(
        par = par_i,
        fn = logL_GP_mod,
        gr = gr_GP_mod,
        db = db_i,
        mean = post_mean_i,
        kern = kern_i,
        post_cov = post_cov_i,
        pen_diag = pen_diag,
        method = "L-BFGS-B",
        control = list(factr = 1e13, maxit = 25)
      )$par %>%
        tibble::as_tibble_row() %>%
        return()
    }
    new_hp_i <- sapply(list_ID,
                       floop,
                       simplify = FALSE,
                       USE.NAMES = TRUE) %>%
      tibble::enframe(name = "ID") %>%
      tidyr::unnest(cols = .data$value)
  }

  list(
    "hp_0" = new_hp_0,
    "hp_i" = new_hp_i
  ) %>%
    return()
}
