#' E-Step of the EM algorithm
#'
#' Expectation step of the EM algorithm to compute the parameters of the
#' hyper-posterior Gaussian distribution of the mean process in Magma.
#'
#' @param db A tibble or data frame. Columns required: ID, Input, Output.
#'    Additional columns for covariates can be specified.
#' @param m_0 A vector, corresponding to the prior mean of the mean GP
#'    (\eqn{\mu_0}).
#' @param kern_0 A kernel function, associated with the mean GP (\eqn{\mu_{0}}).
#' @param kern_i A kernel function, associated with the individual GPs
#'    (\eqn{\Psi_{i}).
#' @param hp_0 A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_0}.
#' @param hp_i A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#' numerical issues when inverting nearly singular matrices.
#'
#' @return A named list, containing the elements \code{mean}, a tibble
#' containing the Input and associated Output of the hyper-posterior mean
#' parameter, and \code{cov}, the hyper-posterior covariance matrix.
#'
#' @examples
#' db <- simu_db(N = 10, common_input = T)
#' m_0 <- rep(0, 10)
#' hp_0 <- hp()
#' hp_i <- hp("SE", list_ID = unique(db$ID))
#' e_step(db, m_0, "SE", "SE", hp_0, hp_i, 0.001)
#'
#' db_async <- simu_db(N = 10, common_input = F)
#' m_0_async <- rep(0, db_async$Input %>% unique() %>% length())
#' e_step(db_async, m_0_async, "SE", "SE", hp_0, hp_i, 0.001)
e_step <- function(db, m_0, kern_0, kern_i, hp_0, hp_i, pen_diag) {
  ## Define the union of all reference Inputs in the dataset
  all_t <- unique(db$Input) %>% sort()
  ## Compute all the inverse covariance matrices
  inv_0 <- kern_to_inv(all_t, kern_0, hp_0, pen_diag)
  list_inv_i <- list_kern_to_inv(db, kern_i, hp_i, pen_diag)
  ## Create a named list of Output values for all individuals
  list_output_i <- base::split(db$Output, list(db$ID))

  ## Update the posterior inverse covariance ##
  new_inv <- inv_0
  for (inv_i in list_inv_i)
  {
    ## Collect the input's common indices between mean and individual processes
    common_times <- intersect(row.names(inv_i), row.names(new_inv))
    ## Sum the common inverse covariance's terms
    new_inv[common_times, common_times] <- new_inv[common_times, common_times] +
      inv_i[common_times, common_times]
  }
  ##############################################

  ## Update the posterior mean ##
  weighted_0 <- inv_0 %*% m_0
  for (i in names(list_inv_i))
  {
    ## Compute the weighted mean for the i-th individual
    weighted_i <- list_inv_i[[i]] %*% list_output_i[[i]]
    ## Collect the input's common indices between mean and individual processes
    common_times <- intersect(row.names(weighted_i), row.names(weighted_0))
    ## Sum the common weighted mean's terms
    weighted_0[common_times, ] <- weighted_0[common_times, ] +
      weighted_i[common_times, ]
  }

  ## Fast or slow matrix inversion if nearly singular
  new_cov <- tryCatch(solve(new_inv), error = function(e) {
    MASS::ginv(new_inv)
  })
  ## Compute the updated mean parameter
  new_mean <- new_cov %*% weighted_0 %>% as.vector()
  ##############################################

  list(
    "mean" = tibble::tibble("Input" = all_t, "Output" = new_mean),
    "cov" = new_cov
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
#' @param old_hp_0 A tibble or data frame, containing the hyper-parameters
#'    from the previous M-step (or initialisation) associated with the mean GP
#'    (\eqn{\mu_{0}}).
#' @param old_hp_i A tibble or data frame, containing the hyper-parameters
#'    from the previous M-step (or initialisation) associated with the
#'    individual GPs (\eqn{f_{i}}).
#' @param mean A tibble, coming out of the E step, containing the Input and
#'    associated Output of the hyper-posterior mean parameter.
#' @param cov A matrix, coming out of the E step, being the hyper-posterior
#'    covariance parameter.
#' @param kern_0 A kernel function, associated with the mean GP (\eqn{\mu_{0}}).
#' @param kern_i A kernel function, associated with the individual GPs
#'    (\eqn{f_{i}).
#' @param m_0 A vector, corresponding to the prior mean of the mean GP
#'    (\eqn{\mu_0}).
#' @param common_hp A logical value, indicating whether the set of
#'    hyper-parameters is assumed to be common to all indiviuals.
#'
#' @return A named list, containing the elements \code{hp_0}, a tibble
#'    containing the hyper-parameters associated with the mean GP,
#'    (\eqn{\mu_{0}}), \code{hp_i}, a tibble containing the hyper-parameters
#'    associated with the individual GPs (\eqn{\mu_{0}}).
#'
#' @examples
#' db <- simu_db(N = 10, common_input = T)
#' m_0 <- rep(0, 10)
#' hp_0 <- hp()
#' hp_i <- hp("SE", list_ID = unique(db$ID))
#' post <- e_step(db, m_0, "SE", "SE", hp_0, hp_i, 0.001)
#'
#' m_step(db, m_0, "SE", "SE", hp_0, hp_i, post$mean, post$cov, F, 0.001)
m_step <- function(db, m_0, kern_0, kern_i, old_hp_0, old_hp_i,
                   mean, cov, common_hp, pen_diag) {

  list_ID <- unique(db$ID)
  list_hp_0 <- old_hp_0 %>% names()
  list_hp_i <- old_hp_i %>%
    dplyr::select(-ID) %>%
    names()

  new_hp_0 <- optimr::opm(
    par = old_hp_0,
    fn = logL_GP_mod,
    gr = gr_GP_mod,
    db = mean,
    mean = m_0,
    kern = kern_0,
    new_cov = cov,
    pen_diag = pen_diag,
    method = "L-BFGS-B",
    control = list(kkt = FALSE)
  ) %>%
    dplyr::select(list_hp_0) %>%
    tibble::as_tibble()

  if (common_hp) {
    new_hp_i <- optimr::opm(
      par = old_hp_i %>% select(- .data$ID) %>% slice(1),
      fn = logL_GP_mod_common_hp,
      gr = gr_GP_mod_common_hp,
      db = db,
      mean = mean,
      kern = kern_i,
      new_cov = cov,
      method = "L-BFGS-B",
      control = list(kkt = F)
    ) %>%
      dplyr::select(list_hp_i) %>%
      tibble::as_tibble() %>%
      dplyr::uncount(weights = length(list_ID)) %>%
      dplyr::mutate(ID = list_ID, .before = 1)
  }
  else {
    floop <- function(i) {
      ## Extract the i-th specific inputs
      input_i <- db %>%
        dplyr::filter(.data$ID == i) %>%
        dplyr::pull(.data$Input)
      ## Extract the mean values associated with the i-th specific inputs
      mean_i = mean %>%
        dplyr::filter(.data$Input %in% input_i) %>%
        dplyr::pull(.data$Output)
      ## Extract the covariance values associated with the i-th specific inputs
      cov_i = cov[as.character(input_i), as.character(input_i)]
      ## Extract the data associated with the i-th individual
      db_i = db %>%
        dplyr::filter(.data$ID == i) %>%
        dplyr::select(- .data$ID)
      ## Extract the hyper-parameters associated with the i-th individual
      par_i = old_hp_i %>%
        filter(.data$ID == i) %>%
        select(- .data$ID)

      optimr::opm(
        par = par_i,
        fn = logL_GP_mod,
        gr = gr_GP_mod,
        db = db_i,
        mean = mean_i,
        kern = kern_i,
        new_cov = cov_i,
        pen_diag = 0,
        method = "L-BFGS-B",
        control = list(kkt = F)
      ) %>%
        dplyr::select(list_hp_i) %>%
        tibble::as_tibble()
      return()
    }
    new_hp_i <- sapply(list_ID, floop, simplify = FALSE, USE.NAMES = TRUE) %>%
      enframe(name = 'ID', value = hp) %>%
      unnest(cols = hp)
  }

  list(
    "hp_0" = new_hp_0,
    "hp_i" = new_hp_i
  ) %>%
    return()
}
