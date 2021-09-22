#' Variation of the E-Step of the EM algorithm
#'
#' Expectation step of the variation of the EM algorithm to compute the parameters of the
#' hyper-posterior Gaussian distribution of the mean process in Magma.
#'
#' @param db A tibble or data frame. Columns required: ID, Input, Output.
#'    Additional columns for covariates can be specified.
#' @param kern_k A kernel function, associated with the mean GP.
#' @param kern_i A kernel function, associated with the individual GPs.
#' @param hp_k A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_k}.
#' @param hp_i A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#' @param m_k prior means of the mu_k processes.
#' @param old_hp_mixture values au hp_mixture from previous iterations.
#'
#' @return A named list, containing the elements \code{mean}, a tibble
#' containing the Input and associated Output of the hyper-posterior mean
#' parameter, \code{cov}, the hyper-posterior covariance matrix,
#' and \code{hp_mixture}, the probability to belong to a cluster for an individual.
#' @export
#'
#' @examples
#' \donttest{
#' k = seq_len(3)
#' m_k <- c("K1" = 0, "K2" = 0, "K3" = 0)
#'
#' db <- simu_db(N = 10, common_input = TRUE)
#' hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
#' hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))
#'
#' old_hp_mixture = MagmaClustR:::ini_hp_mixture(db = db, k = length(k), nstart = 50)
#' prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
#' hp_k[['prop_mixture']] = sapply( prop_mixture_1, function(x) x %>% unlist() %>% mean() )
#'
#' MagmaClustR:::e_step_VEM(db, m_k, "SE", "SE", hp_k, hp_i, old_hp_mixture ,0.001)
#'
#' }
#'
e_step_VEM = function(db,
                      m_k,
                      kern_k,
                      kern_i,
                      hp_k,
                      hp_i,
                      old_hp_mixture,
                      pen_diag) {

  prop_mixture_k = hp_k$prop_mixture
  all_input = unique(db$Input) %>% sort()
  t_clust = tibble::tibble('ID' = rep(names(m_k), each = length(all_input)),
                           'Input' = rep(all_input, length(m_k)))

  ## Compute all the inverse covariance matrices
  list_inv_k = list_kern_to_inv(t_clust, kern_k, hp_k, pen_diag)
  list_inv_i = list_kern_to_inv(db, kern_i, hp_i, pen_diag)

  ## Create a named list of Output values for all individuals
  list_output_i <- base::split(db$Output, list(db$ID))

  ## Update each mu_k parameters for each cluster ##
  floop = function(k)
  {
    new_inv = list_inv_k[[k]]
    hp_mixture = old_hp_mixture[k]
    for(x in list_inv_i %>% names())
    {
      inv_i = list_inv_i[[x]]
      ## Collect the input's common indices between mean and individual processes
      common_times = intersect(row.names(inv_i), row.names(new_inv))
      ## Sum the common inverse covariance's terms
      new_inv[common_times, common_times] = new_inv[common_times, common_times] +
        as.double(hp_mixture[x,]) * inv_i[common_times, common_times]
    }

    ## Fast or slow matrix inversion if nearly singular
    tryCatch(new_inv %>% chol() %>% chol2inv(),
             error = function(e){MASS::ginv(new_inv)}) %>%
      `rownames<-`(all_input) %>%
      `colnames<-`(all_input) %>%
      return()
  }
  cov_k = sapply(names(m_k), floop, simplify = FALSE, USE.NAMES = TRUE)

  ##############################################

  ## Update the posterior mean for each cluster ##

  floop2 = function(k)
  {
    prior_mean = m_k[[k]]; prior_inv = list_inv_k[[k]]
    hp_mixture <- old_hp_mixture[k]

    if(length(prior_mean) == 1){prior_mean = rep(prior_mean, ncol(prior_inv))}
    weighted_mean = prior_inv %*% prior_mean

    for(i in list_inv_i %>% names())
    {
      ## Compute the weighted mean for the i-th individual
      weighted_i = as.double(hp_mixture[i,]) * list_inv_i[[i]] %*% list_output_i[[i]]
      ## Collect the input's common indices between mean and individual processes
      common_times = intersect(row.names(weighted_i), row.names(weighted_mean))
      ## Sum the common weighted mean's terms
      weighted_mean[common_times,] = weighted_mean[common_times,] +
        weighted_i[common_times,]
    }

    ## Compute the updated mean parameter
    new_mean = cov_k[[k]] %*% weighted_mean %>% as.vector()
    tibble::tibble('Input' = all_input, 'Output' = new_mean) %>% return()
  }
  mean_k = sapply(names(m_k), floop2, simplify = FALSE, USE.NAMES = TRUE)

  ## Update hp_mixture
  c_i = 0
  c_k = 0
  mat_elbo = matrix(NA, nrow = length(names(m_k)), ncol = length(unique(db$ID)) )

  for(i in unique(db$ID))
  { c_i = c_i + 1
  input_i = db %>% dplyr::filter(.data$ID == i) %>% dplyr::pull(.data$Input)
  for(k in names(m_k))
  { c_k = c_k + 1

    ## Extract the i-th specific hyper-parameters.
    hp_i_i = hp_i %>% dplyr::filter(.data$ID == i)
    ## Extract the data associated with the i-th individual
    db_i = db %>% dplyr::filter(.data$ID == i) %>% dplyr::select(-.data$ID)
    ## Extract the mean values associated with the i-th specific inputs
    mean_k_i = mean_k[[k]] %>%
      dplyr::filter(.data$Input %in% input_i) %>%
      dplyr::pull(.data$Output)
    ## Extract the covariance values associated with the i-th specific inputs
    cov_k_i = cov_k[[k]][as.character(input_i), as.character(input_i)]

    mat_elbo[c_k,c_i] = - logL_GP_mod(hp_i_i, db_i, mean_k_i , kern_i, cov_k_i, pen_diag)
    if(is.na(mat_elbo[c_k,c_i])){print(i)}
  }
  c_k = 0
  }

  ## We need to use the 'log-sum-exp' trick: exp(x - max(x)) / sum exp(x - max(x)) to remain numerically stable
  mat_L = mat_elbo %>% apply(2,function(x) exp(x - max(x)))

  hp_mixture <- (prop_mixture_k * mat_L) %>% apply(2,function(x) x / sum(x)) %>%
    `rownames<-`(names(m_k)) %>%
    t %>%
    tibble::as_tibble()

  hp_mixture <- tibble::tibble('ID' = unique(db$ID)) %>%
    dplyr::mutate(hp_mixture)

  list('mean' = mean_k,
       'cov' = cov_k,
       'hp_mixture' = hp_mixture) %>%
    return()

}


#' Variation of the M-Step of the EM algorithm
#'
#' Maximization step of the variation of the EM algorithm to compute hyper-parameters of all the
#' kernels involved in Magma.
#'
#' @param db A tibble or data frame. Columns required: ID, Input, Output.
#'    Additional columns for covariates can be specified.
#' @param list_mu_param List of parameters of the K mean GPs.
#' @param kern_k kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0)
#' @param kern_i kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i)
#' @param m_k prior value of the mean parameter of the mean GPs (mu_k). Length = 1 or nrow(unique(db$Input))
#' @param common_hp_k boolean indicating whether hp are common among mean GPs (for each mu_k)
#' @param common_hp_i boolean indicating whether hp are common among individual GPs (for each y_i)
#' @param old_hp_i A named vector, tibble or data frame, containing the hyper-parameters
#'    from the previous M-step (or initialisation) associated with the
#'    individual GPs.
#' @param old_hp_k A named vector, tibble or data frame, containing the hyper-parameters
#'    from the previous M-step (or initialisation) associated with the clusters.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#' numerical issues when inverting nearly singular matrices.
#'
#' @return A named list, containing the elements \code{hp_k}, a tibble
#'    containing the hyper-parameters associated with each cluster,
#'    \code{hp_i}, a tibble containing the hyper-parameters
#'    associated with the individual GPs, and \code{prop_mixture_k},
#'    a tibble containing the hyper-parameters associated with each individual,
#'    indicating in which cluster it belongs.
#'
#' @export
#'
#' @examples
#' \donttest{
#'
#' ## Common inputs across individuals & cluster and differents HPs across individuals & Cluster
#' k = seq_len(2)
#' m_k <- c("K1" = 0, "K2" = 0)
#'
#' db <- simu_db(N = 10, common_input = FALSE)
#' hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
#' hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))
#'
#' old_hp_mixture = MagmaClustR:::ini_hp_mixture(db = db, k = length(k), nstart = 50)
#' prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
#' hp_k[['prop_mixture']] = sapply( prop_mixture_1, function(x) x %>% unlist() %>% mean() )
#'
#' post = MagmaClustR:::e_step_VEM(db, m_k, "SE", "SE", hp_k, hp_i, old_hp_mixture ,0.001)
#'
#' MagmaClustR:::m_step_VEM(db, hp_k, hp_i, post, "SE", "SE", m_k, FALSE, FALSE, 2)
#'
#'
#' ## Different inputs across individuals & cluster and common HPs
#' k = seq_len(4)
#' m_k <- c("K1" = 0, "K2" = 0, "K3" = 0, "K4" = 0)
#' db <- simu_db(N = 10, common_input = FALSE)
#' hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k), common_hp = TRUE)
#' hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID), common_hp = TRUE)
#'
#' old_hp_mixture = MagmaClustR:::ini_hp_mixture(db = db, k = length(k), nstart = 50)
#' prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
#' hp_k[['prop_mixture']] = sapply( prop_mixture_1, function(x) x %>% unlist() %>% mean() )
#'
#' post = MagmaClustR:::e_step_VEM(db, m_k, "SE", "SE", hp_k, hp_i, old_hp_mixture ,0.001)
#'
#' MagmaClustR:::m_step_VEM(db, hp_k, hp_i, post, "SE", "SE", m_k, TRUE, TRUE, 0.1)
#'
#'
#' ## Different kernels
#' kernels <- c("SE", "LIN", "PERIO", "SE + PERIO", "SE * LIN + PERIO" )
#' k = seq_len(2)
#' m_k <- c("K1" = 0, "K2" = 0)
#'
#' for(i in kernels) {
#'
#' ## Common inputs across individuals & cluster and differents HPs across individuals & Cluster
#'
#' db <- simu_db(N = 10, common_input = TRUE)
#' hp_k <- MagmaClustR:::hp(i, list_ID = names(m_k))
#' hp_i <- MagmaClustR:::hp(i, list_ID = unique(db$ID))
#'
#' old_hp_mixture = MagmaClustR:::ini_hp_mixture(db = db, k = length(k), nstart = 50)
#' prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
#' hp_k[['prop_mixture']] = sapply( prop_mixture_1, function(x) x %>% unlist() %>% mean() )
#'
#' post = MagmaClustR:::e_step_VEM(db, m_k, i, i, hp_k, hp_i, old_hp_mixture ,0.001)
#'
#' MagmaClustR:::m_step_VEM(db, hp_k, hp_i, post, i, i, m_k, FALSE, FALSE, 25) -> a
#' paste("kernel =",i ,
#' "Common inputs across individuals & cluster and differents HPs across individuals & Cluster") %>% print()
#' print(a)
#' print(" ")
#'
#'
#' ## if Error in svd(X) appear, increase the pen_diag
#'
#'}
m_step_VEM = function(db, old_hp_k, old_hp_i, list_mu_param, kern_k, kern_i, m_k, common_hp_k, common_hp_i, pen_diag)
{
  list_ID_k = names(m_k)
  list_ID_i = unique(db$ID)

  list_hp_i <- old_hp_i %>%
    dplyr::select(-.data$ID) %>%
    names()

  list_hp_k <- old_hp_k %>%
    dplyr::select(-.data$ID) %>%
    dplyr::select(-.data$prop_mixture) %>%
    names()

  ## Detect whether the kernel_0 provides derivatives for its hyper-parameters
  if (kern_k %>% is.function()) {
    if (!("deriv" %in% methods::formalArgs(kern_k))) {
      gr_GP_mod <- NULL
      gr_GP_mod_common_hp_k <- NULL
    }
  }

  ## Detect whether the kernel_i provides derivatives for its hyper-parameters
  if (kern_i %>% is.function()) {
    if (!("deriv" %in% methods::formalArgs(kern_i))) {
      gr_clust_multi_GP_common_hp_i <- NULL
      gr_clust_multi_GP <- NULL
    }
  }


  ## Check whether hyper-parameters are common to all individuals
  if(common_hp_i)
  {
    ## Extract the hyper-parameters associated with the i-th individual
    par_i <-  old_hp_i %>%
      dplyr::select(-.data$ID) %>%
      dplyr::slice(1)

    ## Optimise hyper-parameters of the individual processes
    new_theta_i <- optimr::opm(
      par = par_i,
      fn = elbo_clust_multi_GP_common_hp_i,
      gr = gr_clust_multi_GP_common_hp_i,
      db = db,
      mu_k_param = list_mu_param,
      kern = kern_i,
      pen_diag = pen_diag,
      method = "L-BFGS-B",
      control = list(kkt = F)
      ) %>%
        dplyr::select(list_hp_i) %>%
        tibble::as_tibble() %>%
        tidyr::uncount(weights = length(list_ID_i)) %>%
        dplyr::mutate('ID' = list_ID_i, .before = 1)

  }
  else {
    loop2 = function(i) {
      ## Extract the hyper-parameters associated with the i-th individual
      par_i <-  old_hp_i %>%
        dplyr::filter(.data$ID == i) %>%
        dplyr::select(-.data$ID)
      ## Extract the data associated with the i-th individual
      db_i <- db %>% dplyr::filter(.data$ID == i)

      ## Optimise hyper-parameters of the individual processes
      optimr::opm(
        par = par_i,
        fn = elbo_clust_multi_GP,
        gr = gr_clust_multi_GP,
        db = db_i,
        pen_diag = pen_diag,
        mu_k_param = list_mu_param,
        kern = kern_i,
        method = "L-BFGS-B",
        control = list(kkt = F)
      ) %>%
        dplyr::select(list_hp_i) %>%
        tibble::as_tibble() %>%
        return()
    }
    new_theta_i = sapply(list_ID_i, loop2, simplify = FALSE, USE.NAMES = TRUE) %>%
      tibble::enframe(name = "ID") %>%
      tidyr::unnest(cols = .data$value)
  }

  ## Check whether hyper-parameters are common to all cluster
  if(common_hp_k)
  {
    ## Extract the hyper-parameters associated with the k-th cluster
    par_k <- old_hp_k %>%
      dplyr::select(-.data$ID) %>%
      dplyr::slice(1) %>%
      dplyr::select(-.data$prop_mixture)

    ## Optimise hyper-parameters of the processes of each cluster
    new_theta_k <- optimr::opm(
      par = par_k,
      fn = elbo_GP_mod_common_hp_k,
      gr = gr_GP_mod_common_hp_k,
      db = list_mu_param$mean,
      mean = m_k,
      kern = kern_k,
      post_cov = list_mu_param$cov,
      pen_diag = pen_diag,
      method = "L-BFGS-B",
      control = list(kkt = F)
      ) %>%
        dplyr::select(list_hp_k) %>%
        tibble::as_tibble() %>%
        tidyr::uncount(weights = length(list_ID_k)) %>%
        dplyr::mutate('ID' = list_ID_k, .before = 1) %>%
        dplyr::mutate("prop_mixture" = pen_diag)


  }
  else
  {
    loop = function(k) {
      ## Extract the hyper-parameters associated with the k-th cluster
      par_k <- old_hp_k %>%
        dplyr::filter(.data$ID == k) %>%
        dplyr::select(-.data$ID) %>%
        dplyr::select(-.data$prop_mixture)
      ## Extract the data associated with the k-th cluster
      db_k <- list_mu_param$mean[[k]]
      ## Extract the mean values associated with the k-th specific inputs
      mean_k <- m_k[[k]]
      ## Extract the covariance values associated with the k-th specific inputs
      post_cov_k <- list_mu_param$cov[[k]]

      ## Optimise hyper-parameters of the processes of each cluster
      c(optimr::opm(
        par = par_k,
        logL_GP_mod,
        gr = gr_GP_mod,
        db = db_k,
        mean = mean_k,
        kern = kern_k,
        post_cov = post_cov_k,
        pen_diag = pen_diag,
        method = "L-BFGS-B",
        control = list(kkt = FALSE)
      ) %>%
        dplyr::select(list_hp_k) %>%
        tibble::as_tibble(),
        pen_diag) %>%
        return()
    }

    new_theta_k = sapply(list_ID_k, loop, simplify = FALSE, USE.NAMES = TRUE) %>%
      tibble::enframe(name = "ID") %>%
      tidyr::unnest_auto(.data$value) %>%
      dplyr::rename_at(dplyr::vars(names(.) %>%
      utils::tail(1)),dplyr::funs(paste('prop_mixture')) )

  }

  prop_mixture_k_1 <- list_mu_param$hp_mixture %>% dplyr::select(-.data$ID)

  prop_mixture_k = sapply( prop_mixture_k_1, function(x) x %>% unlist() %>% mean() ) %>%
    t %>%
    tibble::as_tibble()

  list(
    "hp_k" = new_theta_k,
    "hp_i" = new_theta_i,
    'prop_mixture_k' = prop_mixture_k
  ) %>%
    return()
}
