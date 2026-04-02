# ==========================================================================
# generate_random_config.R
# Fonctions utilitaires partagées par tous les scripts de benchmark NeurIPS
# ==========================================================================

#' Génère une configuration aléatoire respectant les contraintes MOMT
#' Lengthscales identiques entre outputs et clusters
#' @param n_out Nombre d'outputs
#' @param n_clust Nombre de clusters (par défaut 3)
generate_random_config <- function(n_out, n_clust = 3) {
  hp_mp_list <- list()
  jitter_leng <- 0.2
  jitter_var  <- 0.5
  global_l0_min <- log(1 / 1000)
  global_l0_max <- log(1 / 100)

  # Lengthscale unique pour les processus moyens
  safe_max_l0 <- global_l0_max - jitter_leng - 0.05
  center_l0 <- runif(1, global_l0_min, safe_max_l0)
  l0_val <- runif(1, center_l0 - jitter_leng, center_l0 + jitter_leng)
  l0_val <- max(min(l0_val, global_l0_max), global_l0_min)
  lu0_val <- l0_val

  for (o in 1:n_out) {
    center_S0 <- runif(1, log(1), log(50))
    for (k in 1:n_clust) {
      S0_val <- runif(1, center_S0 - jitter_var, center_S0 + jitter_var)
      S0_val <- max(min(S0_val, log(20)), log(1))
      hp_mp_list[[length(hp_mp_list) + 1]] <- tibble(
        temp_cluster_id = k, output_id = o,
        l0_min = l0_val, l0_max = l0_val,
        S0_min = S0_val, S0_max = S0_val,
        lu0_min = lu0_val, lu0_max = lu0_val
      )
    }
  }

  hp_mp_config <- bind_rows(hp_mp_list) %>%
    dplyr::arrange(temp_cluster_id, output_id) %>%
    dplyr::select(-temp_cluster_id)

  # Paramètres tâches : lengthscale unique
  base_task_params <- list()
  lt_val  <- runif(1, log(1 / 1000), log(1 / 100))
  lut_val <- lt_val

  for (o in 1:n_out) {
    St_val <- runif(1, log(0.4), log(1))
    base_task_params[[length(base_task_params) + 1]] <- tibble(
      output_id = o,
      lt_min = lt_val, lt_max = lt_val,
      St_min = St_val, St_max = St_val,
      noise_min = -3, noise_max = -3,
      lu_min = lut_val, lu_max = lut_val
    )
  }
  base_task_df <- bind_rows(base_task_params)

  hp_task_list_full <- list()
  for (k in 1:n_clust) hp_task_list_full[[k]] <- base_task_df
  hp_task_config <- bind_rows(hp_task_list_full)

  return(list(mp = hp_mp_config, task = hp_task_config))
}

#' LogL_GP (définie en dehors du package pour éviter les problèmes de scope)
logL_GP_outside_package <- function(hp, db, mean, kern, post_cov, pen_diag, hp_col_names) {
  if (!(hp %>% tibble::is_tibble()) && length(db$Output_ID %>% unique()) > 1) {
    hp_tibble <- reconstruct_hp(par_vector = hp, hp_names = hp_col_names,
                                output_ids = db$Output_ID %>% unique())
  } else if (!(hp %>% tibble::is_tibble()) && length(db$Output_ID %>% unique()) == 1) {
    hp_tibble <- hp %>% t() %>% tibble::as_tibble() %>% stats::setNames(hp_col_names)
  } else {
    hp_tibble <- hp
  }
  inputs <- db %>% dplyr::select(-Output)
  if (length(db$Output_ID %>% unique()) > 1) {
    cov <- kern_to_cov(inputs, kern, hp_tibble) + post_cov
    inv <- cov %>% chol_inv_jitter(pen_diag = pen_diag)
  } else {
    all_inputs <- db %>% dplyr::select(-c(Output, Output_ID)) %>% unique() %>%
      dplyr::arrange(Reference)
    inv <- kern_to_inv(input = all_inputs, kern = kern, hp = hp_tibble, pen_diag = pen_diag)
  }
  (-dmnorm(db$Output, mean, inv, log = TRUE)) %>% return()
}

#' Scaling des données : Output / max(abs(Output)) par (Output_ID, Cluster_ID)
#' Retourne la liste (data_scaled, scale_factors)
scale_data <- function(data_obs) {
  scale_factors <- data_obs %>%
    dplyr::group_by(Output_ID, Cluster_ID) %>%
    dplyr::summarise(scale_factor = max(abs(Output)), .groups = "drop")

  data_scaled <- data_obs %>%
    dplyr::left_join(scale_factors, by = c("Output_ID", "Cluster_ID")) %>%
    dplyr::mutate(Output = Output / scale_factor) %>%
    dplyr::select(-scale_factor)

  return(list(data_scaled = data_scaled, scale_factors = scale_factors))
}

#' Dé-scaling des prédictions
#' Pour une tâche de test dont on connaît le cluster (via les données simulées),
#' on retrouve les scale_factors correspondants.
#' pred_df doit contenir au moins : Input, Output_ID, Mean, Var
#' task_cluster_id : le Cluster_ID de la tâche (issu des données simulées)
unscale_predictions <- function(pred_df, scale_factors, task_cluster_id) {
  sf_task <- scale_factors %>%
    dplyr::filter(Cluster_ID == task_cluster_id)

  pred_df %>%
    dplyr::left_join(sf_task %>% dplyr::select(Output_ID, scale_factor),
                     by = "Output_ID") %>%
    dplyr::mutate(
      Mean = Mean * scale_factor,
      Var  = Var * scale_factor^2
    ) %>%
    dplyr::select(-scale_factor)
}
