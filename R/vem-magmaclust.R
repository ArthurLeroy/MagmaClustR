#' E-Step of the VEM algorithm
#'
#' Expectation step of the Variational EM algorithm used to compute
#' the parameters of the hyper-posteriors distributions
#' for the mean processes and mixture variables involved in MagmaClust.
#'
#' @param db A tibble or data frame. Columns required: ID, Input_ID, Input,
#'    Output_ID, Output.
#' @param m_k A named list of vectors, corresponding to the prior mean
#'    parameters of the K mean GPs.
#' @param kern_k A kernel function, associated with the K mean GPs.
#' @param kern_t A kernel function, associated with the M individual GPs.
#' @param weight_inv_k A number, indicating the weight that the user wants to
#'  attribute to the inverse prior covariances inv_k.
#' @param hp_k A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_k}.
#' @param hp_t A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern_i}.
#' @param old_mixture A list of mixture values from the previous iteration.
#' @param iter A number, indicating the current iteration of the VEM algorithm.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#'    numerical issues when inverting nearly singular matrices.
#'
#' @return A named list, containing the elements \code{mean}, a tibble
#'    containing the Input and associated Output of the hyper-posterior mean
#'    parameters, \code{cov}, the hyper-posterior covariance matrices,
#'    and \code{mixture}, the probabilities to belong to each cluster for each
#'    individual.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
ve_step <- function(db,
                    m_k,
                    kern_k,
                    kern_t,
                    weight_inv_k,
                    hp_k,
                    hp_t,
                    old_mixture,
                    iter,
                    pen_diag) {

  list_ID_task <- unique(db$Task_ID)
  list_output_ID <-  db$Output_ID %>% unique()

  # Get the union of all unique input points from the training data
  all_inputs <- db %>%
    dplyr::select(-c(Task_ID, Output)) %>%
    unique() %>%
    tidyr::separate(Reference,
                    into = c("Output_ID_temp", "Input_temp"),
                    sep = ";",
                    remove = FALSE) %>%
    dplyr::mutate(Input_temp_numeric = as.numeric(Input_temp)) %>%
    dplyr::arrange(Output_ID_temp, Input_temp_numeric) %>%
    dplyr::select(c(Input_1, Reference, Output_ID))

  all_input <- all_inputs %>%
    dplyr::pull(Reference)

  ## Sort the database according to Reference
  ## WARNING: ARRRANGE !!!!!
  # db <- db %>% dplyr::arrange(Reference, .by_group = TRUE)

  prop_mixture_k <- hp_k %>%
    dplyr::select(c(Cluster_ID, Output_ID, prop_mixture))

  ## Format a sequence of inputs for all clusters
  t_clust <- tidyr::expand_grid("Cluster_ID" = names(m_k),
                                all_inputs
  )

  list_inv_k <- list()
  list_inv_t <- list()

  # Loop over clusters
  for(k in t_clust$Cluster_ID %>% unique){
    # Subset t_clust and hp_k on k cluster
    t_clust_k <- t_clust %>%
      dplyr::filter(Cluster_ID == k)

    hp_k_subset <- hp_k %>%
      dplyr::filter(Cluster_ID == k)

    # Compute the covariance matrix of the mean process of the
    # k cluster
    cov_k <- suppressWarnings(kern_to_cov(input = t_clust_k,
                         kern = kern_k,
                         hp = hp_k_subset %>%
                           dplyr::select(-c(Cluster_ID, prop_mixture))
                         )
    )

    references <- rownames(cov_k)
    inv_k <- cov_k %>% chol_inv_jitter(pen_diag = pen_diag)

    # Re-apply the stored names to the inverted matrix
    dimnames(inv_k) <- list(references, references)
    inv_k <- weight_inv_k * inv_k
    list_inv_k[[k]] <- inv_k
  }


  # Loop over tasks
  for (t in list_ID_task) {
    # Isolate the data and HPs for the current task
    db_t <- db %>% dplyr::filter(Task_ID == t) %>%
      dplyr::select(-c(Output, Task_ID))
    hp_t_indiv <- hp_t %>% dplyr::filter(Task_ID == t)

    if(length(list_output_ID) > 1 && !(kern_t %>% is.character())){
      # MO case with dependent outputs
      # Call kern_to_cov directly.
      # It will handle the multi-output structure and the noise addition internally.
      # 'kern_t' is expected to be the 'convolution_kernel' function.
      K_task_t <- kern_to_cov(
        input = db_t,
        kern = kern_t,
        hp = hp_t_indiv
      )
    } else if (length(list_output_ID) > 1 && kern_t %>% is.character()){
      # MO case with independent outputs
      hp_t_indiv <- hp_t_indiv %>% dplyr::select(-c(Task_ID, Output_ID))

      K_task_t <- kern_to_cov(
        input = db_t,
        kern = kern_t,
        hp = hp_t_indiv
      )
    } else{
      # Extract all_inputs to call kern_to_cov() on the single output case
      all_inputs_t <- db %>%
        dplyr::filter(Task_ID == t) %>%
        dplyr::select(-c(Task_ID, Output, Output_ID)) %>%
        unique()

      K_task_t <- kern_to_cov(
        input = all_inputs_t,
        kern = kern_t,
        hp = hp_t_indiv
      )
    }

    # Store the correct row/column names before they are lost during inversion
    task_references <- rownames(K_task_t)

    # Invert the covariance matrix (this strips the names)
    K_inv_t <- K_task_t %>% chol_inv_jitter(pen_diag = pen_diag)

    # Re-apply the stored names to the inverted matrix
    dimnames(K_inv_t) <- list(task_references, task_references)

    # Add the inverted matrix to the list
    # The rownames are already correctly set by kern_to_cov
    list_inv_t[[t]] <- K_inv_t
  }

  ## Create a named list of Output values for all individuals
  list_output_t <- base::split(db$Output, list(db$Task_ID))

  ## Update each mu_k parameters for each cluster ##
  floop <- function(k) {
    post_inv <- list_inv_k[[k]]
    tau_k <- old_mixture %>% dplyr::select(Task_ID, k)
    for (t in list_inv_t %>% names())
    {
      ## Extract the corresponding mixture probability
      tau_t_k <- tau_k %>%
        dplyr::filter(Task_ID == t) %>%
        dplyr::pull(k)

      inv_t <- list_inv_t[[t]]
      ## Collect input's common indices between mean and individual processes
      co_input <- intersect(row.names(inv_t), row.names(post_inv))
      ## Sum the common inverse covariance's terms
      post_inv[co_input, co_input] <- post_inv[co_input, co_input] +
        tau_t_k * inv_t[co_input, co_input]
    }

    post_inv %>%
      chol_inv_jitter(pen_diag = pen_diag) %>%
      `rownames<-`(all_input) %>%
      `colnames<-`(all_input) %>%
      return()
  }
  cov_k <- sapply(tidyselect::all_of(names(m_k)),
                  floop,
                  simplify = FALSE,
                  USE.NAMES = TRUE)

  ## Update the posterior mean for each cluster ##
  floop2 <- function(k) {
    prior_mean <- m_k[[k]]
    prior_inv <- list_inv_k[[k]]
    tau_k <- old_mixture %>% dplyr::select(Task_ID, k)

    weighted_mean <- prior_inv %*% prior_mean

    for (t in list_inv_t %>% names())
    {
      ## Extract the corresponding mixture probability
      tau_t_k <- tau_k %>%
        dplyr::filter(Task_ID == t) %>%
        dplyr::pull(k)
      ## Compute the weighted mean for the i-th individual
      weighted_t <- tau_t_k * list_inv_t[[t]] %*% list_output_t[[t]]
      ## Collect input's common indices between mean and individual processes
      co_input <- intersect(row.names(weighted_t), row.names(weighted_mean))
      ## Sum the common weighted mean's terms
      weighted_mean[co_input, ] <- weighted_mean[co_input, ] +
        weighted_t[co_input, ]
    }

    ## Compute the updated mean parameter
    new_mean <- cov_k[[k]] %*% weighted_mean %>% as.vector()
    tibble::tibble(all_inputs,
                   "Output" = new_mean) %>% return()
  }
  mean_k <- sapply(names(m_k), floop2, simplify = FALSE, USE.NAMES = TRUE)

  ## Update mixture (skip first iteration to avoid bad HP initialisation issues)
  if(iter == 1){
    mixture <- old_mixture
  } else{
    mixture <- update_mixture(
      db = db,
      mean_k,
      cov_k,
      hp = hp_t,
      kern = kern_t,
      prop_mixture_k,
      pen_diag
    )
  }

  list(
    "mean" = mean_k,
    "cov" = cov_k,
    "mixture" = mixture
  ) %>%
    return()
}


#' M-Step of the VEM algorithm
#'
#' Maximization step of the Variational EM algorithm used to compute
#' hyper-parameters of all the kernels involved in MagmaClust.
#'
#' @param db A tibble or data frame. Columns required: Task_ID, Input_ID, Input,
#'    Output_ID, Output.
#' @param list_mu_param List of parameters of the K mean GPs.
#' @param kern_k A kernel used to compute the covariance matrix of the mean GP
#'    at corresponding timestamps.
#' @param kern_t A kernel used to compute the covariance matrix of individuals
#'    GP at corresponding timestamps.
#' @param m_k A named list of prior mean parameters for the K mean GPs.
#'    Length = 1 or length(unique(db$Output_ID) or
#'    nb_cluster*length(unique(db$Output_ID).
#' @param shared_hp_clusts A boolean indicating whether hyper-parameters are common
#'    among the mean GPs.
#' @param shared_hp_tasks A boolean indicating whether hyper-parameters are common
#'    among the individual GPs.
#' @param old_hp_t A named vector, tibble or data frame, containing the
#'    hyper-parameters from the previous  M-step (or initialisation) associated
#'    with the individual GPs.
#' @param old_hp_k A named vector, tibble or data frame, containing the
#'    hyper-parameters from the previous M-step (or initialisation) associated
#'    with the mean GPs.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#' numerical issues when inverting nearly singular matrices.
#'
#' @return A named list, containing the elements \code{hp_k}, a tibble
#'    containing the hyper-parameters associated with each cluster,
#'    \code{hp_t}, a tibble containing the hyper-parameters
#'    associated with the individual GPs, and \code{prop_mixture_k},
#'    a tibble containing the hyper-parameters associated with each task,
#'    indicating the probabilities to belong to each cluster.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
vm_step <- function(db,
                    old_hp_k,
                    old_hp_t,
                    list_mu_param,
                    kern_k,
                    kern_t,
                    m_k,
                    pen_diag,
                    shared_hp_tasks) {

  # browser()
  list_ID_k <- names(m_k)
  list_ID_t <- unique(db$Task_ID)
  output_ids_vector <- unique(db$Output_ID)

  list_hp_t <- old_hp_t %>%
    dplyr::select(-Task_ID) %>%
    names()

  list_hp_k <- old_hp_k %>%
    dplyr::select(-Cluster_ID) %>%
    dplyr::select(-prop_mixture) %>%
    names()

  ## Detect whether kernel_k provides derivatives for its hyper-parameters
  if (kern_k %>% is.function()) {
    if (!("deriv" %in% methods::formalArgs(kern_k))) {
      gr_GP_mod <- NULL
      gr_GP_mod_common_hp_k <- NULL
    }
  }

  ## Detect whether kernel_t provides derivatives for its hyper-parameters
  if (kern_t %>% is.function()) {
    if (!("deriv" %in% methods::formalArgs(kern_t))) {
      gr_clust_multi_GP_common_hp_t <- NULL
      gr_clust_multi_GP <- NULL
    }
  }

  ## Check whether hyper-parameters are common to all individuals
  if (shared_hp_tasks) {
    if (length(db$Output_ID %>% unique()) > 1){
      # Prepare parameters for optim() in MO case
      hp_per_output <- old_hp_t %>%
        dplyr::group_by(Output_ID) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Task_ID, -l_u_t) %>%
        tidyr::pivot_longer(cols = -Output_ID,
                            names_to = "hp_name",
                            values_to = "value") %>%
        dplyr::mutate(specific_name = paste(hp_name, Output_ID, sep = "_")) %>%
        dplyr::select(specific_name, value) %>%
        tibble::deframe()

      shared_hp_l_u_t <- old_hp_t$l_u_t[1]
      names(shared_hp_l_u_t) <- "l_u_t"
      par_t <- c(hp_per_output, shared_hp_l_u_t)
      hp_col_names <- names(par_t)
    } else {
      # Prepare parameters for optim() in single output case
      par_t <- old_hp_t %>%
        dplyr::select(-c(Task_ID, Output_ID)) %>%
        dplyr::slice(1)
      hp_col_names <- names(par)
    }

    ## Optimise hyper-parameters of the individual processes
    new_hp_t <- stats::optim(
      par = par_t,
      fn = elbo_clust_multi_GP_shared_hp_tasks,
      gr = gr_clust_multi_GP_shared_hp_tasks,
      db = db,
      hyperpost = list_mu_param,
      kern = kern_t,
      pen_diag = pen_diag,
      method = "L-BFGS-B",
      hp_col_names = hp_col_names,
      output_ids = output_ids_vector,
      control = list(factr = 1e13, maxit = 25)
    )$par %>%
      tibble::as_tibble_row()

    # Reshape results
    if (length(db$Output_ID %>% unique()) > 1){
      # MO case
      new_hp_t <- new_hp_t %>%
        tidyr::pivot_longer(
          cols = -dplyr::any_of("l_u_t"),
          names_to = c(".value", "Output_ID"),
          names_pattern = "(.+)_(\\d+)$"
        ) %>%
        tidyr::crossing(Task_ID = list_ID_t, .)
    } else {
      # Single output case
      new_hp_t$Output_ID <- "1"
      new_hp_t <- new_hp_t %>%
        dplyr::mutate(Task_ID = list(list_ID_t)) %>%
        tidyr::unnest(Task_ID)
      new_hp_t <- new_hp_t
    }
  } else {
    cat("HPs are task-specific.\n")

    floop <- function(t) {
      cat(paste0("Optimizing for task ", t, "...\n"))
      # Filter data for the current task
      db_t <- db %>% dplyr::filter(Task_ID == t) %>%
        dplyr::select(-Task_ID)

      input_t <- db_t %>% dplyr::pull(Reference)

      post_mean_t <- post_mean %>%
        dplyr::filter(Reference %in% input_t) %>%
        dplyr::pull(Output)

      post_cov_t <- post_cov[as.character(input_t), as.character(input_t)]

      old_hp_t_task <- old_hp_t %>% dplyr::filter(Task_ID == t)

      # Prepare parameters for optim() for this specific task
      if (length(db$Output_ID %>% unique()) > 1){
        hp_per_output <- old_hp_t_task %>%
          dplyr::select(-Task_ID, -l_u_t) %>%
          tidyr::pivot_longer(cols = -Output_ID,
                              names_to = "hp_name",
                              values_to = "value") %>%
          dplyr::mutate(specific_name = paste(hp_name, Output_ID, sep = "_")) %>%
          dplyr::select(specific_name, value) %>%
          tibble::deframe()

        shared_hp <- old_hp_t_task$l_u_t[1]
        names(shared_hp) <- "l_u_t"

        par_t <- c(hp_per_output, shared_hp)
        hp_col_names <- names(par_t)
      } else {
        ## Extract the hyper-parameters associated with the i-th individual
        par_t <- old_hp_t %>%
          dplyr::filter(Task_ID == t) %>%
          dplyr::select(-c(Task_ID, Output_ID))
        hp_col_names <- names(par_t)
      }

      ## Optimise hyper-parameters of the individual processes
      stats::optim(
        par = par_t,
        fn = elbo_clust_multi_GP,
        gr = gr_clust_multi_GP,
        db = db_t,
        pen_diag = pen_diag,
        hyperpost = list_mu_param,
        kern = kern_t,
        method = "L-BFGS-B",
        hp_col_names = hp_col_names,
        output_ids = output_ids,
        control = list(factr = 1e13, maxit = 25)
      )$par %>%
        tibble::as_tibble_row()
    }
    new_hp_t <- sapply(list_ID_t, loop2, simplify = FALSE, USE.NAMES = TRUE) %>%
      tibble::enframe(name = "Task_ID") %>%
      tidyr::unnest(cols = value)
  }

  ## Re-attribute each set of task specific HPs to the Cluster_ID in which
  ## the task belongs
  # new_hp_t$Cluster_ID <- old_hp_t$Cluster_ID

  ## Compute the prop mixture of each cluster
  prop_mixture <- list_mu_param$mixture %>%
    dplyr::select(-Task_ID) %>%
    colMeans()

  list(
    "hp_k" = old_hp_k,
    "hp_t" = new_hp_t
  ) %>%
    return()
}

#' Update the mixture probabilities for each individual and each cluster
#'
#' @param db A tibble or data frame. Columns required: \code{ID},
#'    \code{Input}, \code{Output}. Additional columns for covariates can be
#'    specified.
#' @param mean_k A list of the K hyper-posterior mean parameters.
#' @param cov_k A list of the K hyper-posterior covariance matrices.
#' @param hp A named vector, tibble or data frame of hyper-parameters
#'    associated with \code{kern}, the individual process' kernel. The
#'    columns/elements should be named according to the hyper-parameters
#'    that are used in \code{kern}.
#' @param kern A kernel function, defining the covariance structure of
#'    the individual GPs.
#' @param prop_mixture A tibble containing the hyper-parameters associated
#'    with each individual, indicating in which cluster it belongs.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#' numerical issues when inverting nearly singular matrices.
#'
#' @return Compute the hyper-posterior multinomial distributions by updating
#'    mixture probabilities.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
update_mixture <- function(db,
                           mean_k,
                           cov_k,
                           hp,
                           kern,
                           prop_mixture,
                           pen_diag) {

  # browser()
  # 1. Initialisation des identifiants
  ID_t <- unique(db$Task_ID)
  ID_k <- names(mean_k) # c("K1", "K2", ...)
  output_vector_ids <- unique(db$Output_ID)

  # 2. Extraction du vecteur des proportions (Optimisation)
  # Au lieu de le chercher dans la boucle, on le prépare une fois pour toutes.
  # prop_mixture est un tibble: Cluster_ID, Output_ID, prop_mixture
  # Comme la prop est la même pour tous les Output_ID d'un cluster, on prend la première occurence.

  if(prop_mixture %>% is_tibble()){
    vec_prop_map <- prop_mixture %>%
      dplyr::select(Cluster_ID, prop_mixture) %>%
      dplyr::distinct() %>%
      tibble::deframe() # Crée un vecteur nommé c("K1" = 0.5, "K2" = 0.5)

    # On s'assure que l'ordre correspond à ID_k
    vec_prop <- vec_prop_map[ID_k]
  } else {
    vec_prop <- prop_mixture
  }

  # Matrice pour stocker les log-vraisemblances p(y | Z=k)
  # Lignes = Clusters, Colonnes = Tâches
  mat_logL <- matrix(NA, nrow = length(ID_k), ncol = length(ID_t))

  # 3. Boucle sur les tâches
  # (On peut utiliser seq_along pour éviter les compteurs manuels c_t)
  for (i_t in seq_along(ID_t)) {

    t <- ID_t[i_t]

    ## Extract data for task t
    # db_t contient TOUS les outputs pour la tâche t (c'est crucial pour le MO)
    db_t <- db %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::select(-Task_ID)

    # input_t doit contenir les références de TOUS les points de TOUS les outputs
    # pour construire la bonne matrice de covariance jointe
    input_t <- db_t %>% dplyr::pull(Reference)

    ## Extract hyper-parameters for task t
    hp_t <- hp %>% dplyr::filter(Task_ID == t)
    hp_col_names <- colnames(hp_t)

    # 4. Boucle sur les clusters
    for (i_k in seq_along(ID_k)) {

      k <- ID_k[i_k]

      ## Extract Mean Process for Cluster k
      # On filtre sur 'input_t' pour récupérer les moyennes correspondant exactement
      # aux points observés de la tâche (Output 1 ET Output 2...)
      # browser()
      mean_k_t <- mean_k[[k]] %>%
        dplyr::filter(Reference %in% input_t) %>%
        # Assurons-nous que l'ordre est le même que dans db_t
        # dplyr::arrange(match(Reference, input_t)) %>%
        dplyr::pull(Output)

      ## Extract Covariance Process for Cluster k
      # C'est ici que la covariance CROISÉE est récupérée si 'cov_k' est bien construite
      cov_k_t <- cov_k[[k]][as.character(input_t), as.character(input_t)]

      ## Calcul de la Log-Vraisemblance Multi-Output
      # logL_GP_mod doit être capable de gérer des vecteurs/matrices MO
      mat_logL[i_k, i_t] <- - logL_GP_mod(
        hp = hp_t,
        db = db_t,
        mean = mean_k_t,
        kern = kern,
        post_cov = cov_k_t,
        pen_diag = pen_diag,
        hp_col_names = hp_col_names,
        output_ids = output_vector_ids
      )
    }
  }

  # 5. Calcul des Responsabilités (Softmax avec Log-Sum-Exp trick)
  # mat_logL contient log p(y | k)
  # On veut : tau_ik = pi_k * p(y|k) / sum(...)
  # Soit en log : log(tau) = log(pi_k) + log p(y|k) - log(sum(...))

  # Ajout du log-prior (log pi_k) à la log-vraisemblance
  log_numerator <- mat_logL + log(vec_prop)

  # Application du Log-Sum-Exp par colonne (par tâche) pour normaliser
  # T_ik = exp( Num_ik - max_j(Num_ij) ) / sum_j( exp( Num_ij - max ) )

  mat_tau <- apply(log_numerator, 2, function(col_scores) {
    max_score <- max(col_scores)
    exp_scores <- exp(col_scores - max_score)
    return(exp_scores / sum(exp_scores))
  })

  # 6. Mise en forme du résultat
  mat_tau %>%
    t() %>% # Transpose pour avoir Tâches x Clusters
    round(5) %>%
    tibble::as_tibble() %>%
    `colnames<-`(ID_k) %>% # Nomme les colonnes K1, K2...
    dplyr::mutate("Task_ID" = ID_t, .before = 1) %>%
    return()
}
