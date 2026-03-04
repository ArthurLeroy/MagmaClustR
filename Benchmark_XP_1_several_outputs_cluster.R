# ==========================================================================
# Benchmark_XP_1_several_outputs_cluster.R
# Version adaptée pour le cluster IMT (SLURM) — PARAMÉTRÉ
#
# Usage :
#   Rscript Benchmark_XP_1_several_outputs_cluster.R --n_out=3 --n_train=15
#
# Ce script traite UNE SEULE combinaison (n_out, n_train) avec 100 itérations.
# Le job SLURM lance 15 instances en parallèle (5 n_out × 3 n_train).
#
# REPRISE AUTOMATIQUE : si le script est relancé, les itérations déjà
# terminées (fichiers .rds présents) sont automatiquement sautées.
# Les sauvegardes sont atomiques (écriture .tmp puis renommage).
# ==========================================================================

# --- 0. PARSING ARGUMENTS ---
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, name, default = NULL) {
  pattern <- paste0("^--", name, "=")
  match <- grep(pattern, args, value = TRUE)
  if (length(match) > 0) {
    return(as.numeric(sub(pattern, "", match[1])))
  }
  if (!is.null(default)) return(default)
  stop(paste0("Argument requis manquant : --", name))
}

N_OUT_TARGET   <- parse_arg(args, "n_out")
N_TRAIN_TARGET <- parse_arg(args, "n_train")

# --- 1. SETUP & LIBRARIES ---
username <- Sys.getenv("USER")
pkg_dir  <- file.path("/scratch", username, "MagmaClustR")
base_dir <- file.path("/scratch", username, "NeurIPS_experiments", "Experience_1")

setwd(pkg_dir)

# library(Metrics)
library(mvtnorm)
library(tidyverse)
# library(devtools)
library(matrixStats)
library(MagmaClustR)
# load_all()

convolution_kernel <- MagmaClustR:::convolution_kernel

cat(paste0("=== Benchmark XP1 MOMT ===\n"))
cat(paste0("  n_out   = ", N_OUT_TARGET, "\n"))
cat(paste0("  n_train = ", N_TRAIN_TARGET, "\n"))
cat(paste0("  Date    : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste0("  Host    : ", Sys.info()["nodename"], "\n"))
cat(paste0("  PID     : ", Sys.getpid(), "\n\n"))
flush.console()

# --- 2. FONCTIONS UTILITAIRES ---
generate_random_config <- function(n_out, n_clust = 3) {

  hp_mp_list <- list()
  jitter_leng <- 0.2
  jitter_var <- 0.5
  global_l0_min <- log(1/1000)
  global_l0_max <- log(1/100)
  global_lu0_max <- log(1/100)

  for(o in 1:n_out){
    safe_max_l0 <- global_l0_max - jitter_leng - 0.05
    center_l0 <- runif(1, global_l0_min, safe_max_l0)
    center_S0 <- runif(1, log(1), log(50))

    valid_center_lu0 <- FALSE
    while(!valid_center_lu0){
      center_lu0 <- runif(1, log(1/300), global_lu0_max)
      if(center_lu0 > center_l0 + jitter_leng) valid_center_lu0 <- TRUE
    }

    for(k in 1:n_clust){
      l0_val <- runif(1, center_l0 - jitter_leng, center_l0 + jitter_leng)
      l0_val <- max(min(l0_val, global_l0_max), global_l0_min)

      S0_min_bound <- log(1)
      S0_max_bound <- log(20)
      S0_val <- runif(1, center_S0 - jitter_var, center_S0 + jitter_var)
      S0_val <- max(min(S0_val, S0_max_bound), S0_min_bound)

      lu0_min_bound <- log(1/300)
      valid_lu0 <- FALSE
      counter <- 0
      while(!valid_lu0){
        lu0_val <- runif(1, center_lu0 - jitter_leng, center_lu0 + jitter_leng)
        lu0_val <- max(min(lu0_val, global_lu0_max), lu0_min_bound)
        if(lu0_val > l0_val) valid_lu0 <- TRUE
        counter <- counter + 1
        if(counter > 100) {
          lu0_val <- min(l0_val + 0.01, global_lu0_max)
          valid_lu0 <- TRUE
        }
      }

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

  base_task_params <- list()
  for(o in 1:n_out){
    lt_val <- runif(1, log(1/1000), log(1/100))
    St_val <- runif(1, log(0.4), log(1))
    valid_lut <- FALSE
    while(!valid_lut){
      lut_val <- runif(1, log(1/300), log(1/100))
      if(lut_val > lt_val) valid_lut <- TRUE
    }
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
  for(k in 1:n_clust) hp_task_list_full[[k]] <- base_task_df
  hp_task_config <- bind_rows(hp_task_list_full)

  return(list(mp = hp_mp_config, task = hp_task_config))
}

logL_GP_outside_package <- function(hp, db, mean, kern, post_cov, pen_diag, hp_col_names) {
  if(!(hp %>% tibble::is_tibble()) && length(db$Output_ID %>% unique()) > 1){
    hp_tibble <- reconstruct_hp(par_vector = hp, hp_names = hp_col_names, output_ids = db$Output_ID %>% unique())
  } else if (!(hp %>% tibble::is_tibble()) && length(db$Output_ID %>% unique()) == 1){
    hp_tibble <- hp %>% t() %>% tibble::as_tibble() %>% stats::setNames(hp_col_names)
  } else {
    hp_tibble <- hp
  }
  inputs <- db %>% dplyr::select(-Output)
  if(length(db$Output_ID %>% unique()) > 1){
    cov <- kern_to_cov(inputs, kern, hp_tibble) + post_cov
    inv <- cov %>% chol_inv_jitter(pen_diag = pen_diag)
  } else {
    all_inputs <- db %>% dplyr::select(-c(Output, Output_ID)) %>% unique() %>% dplyr::arrange(Reference)
    inv <- kern_to_inv(input = all_inputs, kern = kern, hp = hp_tibble, pen_diag = pen_diag)
  }
  (-dmnorm(db$Output, mean, inv, log = T)) %>% return()
}


# --- 3. CONFIGURATION ---
n_out <- N_OUT_TARGET
n_iterations <- 100
n_pred_max <- 1000
n_pred_subsets <- c(1, 10, 100, 1000)
n_train <- N_TRAIN_TARGET
n_clusters <- 3
n_points_fixed <- 30
n_tasks_total <- n_train + n_pred_max

# Sauvegarde des métadonnées
run_info <- list(
  script = "Benchmark_XP_1_several_outputs_cluster.R",
  n_out = n_out, n_train = n_train, n_iterations = n_iterations,
  started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  hostname = as.character(Sys.info()["nodename"]),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "N/A"),
  r_version = R.version.string, pid = Sys.getpid()
)
run_info_dir <- file.path(base_dir, "run_info")
if (!dir.exists(run_info_dir)) dir.create(run_info_dir, recursive = TRUE)
run_info_file <- file.path(run_info_dir,
  paste0("momt_nout", n_out, "_ntrain", n_train, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"))
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)

# Création de l'arborescence
for (n_pred in n_pred_subsets) {
  task_folder_name <- paste0("train_", n_train, "_pred_", n_pred)
  for (subdir in c("Datasets", "Models_MOMT", "Predictions_MOMT")) {
    d <- file.path(base_dir, subdir, paste0("n_out_", n_out), task_folder_name)
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
}

# --- 3b. DÉTECTION DES ITÉRATIONS DÉJÀ TERMINÉES (reprise après crash) ---
completed_iters <- c()
for (iter_check in 1:n_iterations) {
  # Une itération est considérée terminée si le fichier predictions pour n_pred=1 existe
  pred_file_check <- file.path(base_dir, "Predictions_MOMT", paste0("n_out_", n_out),
                               paste0("train_", n_train, "_pred_1"),
                               paste0("predictions_iter_", iter_check, ".rds"))
  if (file.exists(pred_file_check)) {
    completed_iters <- c(completed_iters, iter_check)
  }
}

if (length(completed_iters) > 0) {
  remaining <- setdiff(1:n_iterations, completed_iters)
  cat(paste0("  Itérations déjà terminées : ", length(completed_iters), "/", n_iterations, "\n"))
  if (length(remaining) > 0) {
    cat(paste0("  Reprise à partir de l'itération ", min(remaining), "\n\n"))
  } else {
    cat("  Toutes les itérations sont déjà terminées. Rien à faire.\n\n")
  }
  flush.console()
} else {
  cat(paste0("  Aucune itération précédente trouvée, démarrage depuis le début.\n\n"))
  flush.console()
}

# --- 4. BOUCLE PRINCIPALE (sur les itérations uniquement) ---
t_global_start <- Sys.time()

for (iter in 1:n_iterations) {
  # --- Reprise : skip si itération déjà terminée ---
  if (iter %in% completed_iters) {
    cat(paste0("[n_out=", n_out, " n_train=", n_train, "] Iter ", iter,
               "/", n_iterations, " déjà terminée, skip.\n"))
    flush.console()
    next
  }

  if(exists("trained_model")) rm(trained_model)
  if(exists("pred_res")) rm(pred_res)

  tryCatch({
    current_seed <- n_out * 10000 + n_train * 100 + iter
    set.seed(current_seed)

    cat(paste0("\n[n_out=", n_out, " n_train=", n_train, "] Iter ", iter, "/", n_iterations,
               " (seed=", current_seed, ") [", format(Sys.time(), "%H:%M:%S"), "]\n"))
    flush.console()

    # --- A. Génération Config ---
    configs <- generate_random_config(n_out, n_clust = n_clusters)
    my_priors <- rep(0, n_out * n_clusters)

    t_simu_start <- Sys.time()
    sim_results <- simulate_magmaclust_data_convol(
      nb_clusters = n_clusters,
      total_tasks = n_tasks_total,
      points_per_output_grid = rep(200, n_out),
      grid_ranges = rep(list(c(-1, 1)), n_out),
      shared_hp_clusts = FALSE,
      hp_config_mean_process = configs$mp,
      prior_means = my_priors,
      shared_hp_tasks = TRUE,
      hp_config_tasks = configs$task,
      n_points_per_task_range = c(n_points_fixed, n_points_fixed),
      shared_hp_outputs = FALSE,
      shared_grid_outputs = FALSE,
      shared_grid_clusters = TRUE
    )
    duration_simu <- as.numeric(difftime(Sys.time(), t_simu_start, units = "secs"))
    cat(paste0("  Simulation: ", round(duration_simu, 1), "s\n"))
    flush.console()

    data_obs <- sim_results$simulated_data_df

    # --- B. Split Train / Test ---
    all_task_ids <- unique(data_obs$Task_ID)
    n_tasks_available <- length(all_task_ids)

    if (n_tasks_available < n_tasks_total) {
      n_pred_max_actual <- n_tasks_available - n_train
      if (n_pred_max_actual < 1) stop("Pas assez de tâches!")
    } else {
      n_pred_max_actual <- n_pred_max
    }

    shuffled_task_ids <- as.character(sample(all_task_ids, size = n_tasks_available))
    train_task_ids <- shuffled_task_ids[1:n_train]
    test_task_ids_max <- shuffled_task_ids[(n_train + 1):(n_train + n_pred_max_actual)]

    start_1_global <- -0.7
    end_1_global   <- -0.4

    train_data_model <- data_obs %>%
      filter(Task_ID %in% train_task_ids) %>%
      filter(Output_ID != "1" | (Output_ID == "1" & !(Input >= start_1_global & Input <= end_1_global)))

    pred_tasks_data_all <- list()
    test_tasks_data_all <- list()

    for (test_task_id in test_task_ids_max) {
      task_target_O1 <- data_obs %>%
        filter(Task_ID == test_task_id, Output_ID == "1") %>% arrange(Input)

      indices_in_gap_1 <- which(task_target_O1$Input >= start_1_global & task_target_O1$Input <= end_1_global)
      n_in_gap_1 <- length(indices_in_gap_1)
      points_needed_from_end <- 10 - n_in_gap_1

      if (points_needed_from_end <= 0) {
        final_test_indices <- indices_in_gap_1[1:min(n_in_gap_1, 10)]
      } else {
        remaining_indices <- setdiff(1:nrow(task_target_O1), indices_in_gap_1)
        last_points_indices <- remaining_indices[order(task_target_O1$Input[remaining_indices], decreasing = TRUE)]
        indices_from_end <- last_points_indices[1:min(length(last_points_indices), points_needed_from_end)]
        final_test_indices <- c(indices_in_gap_1, indices_from_end)
      }

      test_data_truth <- task_target_O1[final_test_indices, ] %>% dplyr::arrange(Input)
      test_tasks_data_all[[test_task_id]] <- test_data_truth

      pred_task_input_data <- bind_rows(
        data_obs %>% filter(Task_ID == test_task_id, Output_ID != "1"),
        task_target_O1[-final_test_indices, ]
      ) %>% dplyr::arrange(Output_ID, Input)
      pred_tasks_data_all[[test_task_id]] <- pred_task_input_data
    }

    # --- C. Pré-processing : Initialisation GP Vanille ---
    hp_k_extracted_list <- list()
    max_lengthscale_observed <- -Inf
    n_unique_tasks <- length(unique(train_data_model$Task_ID))

    if (n_unique_tasks <= n_clusters) {
      unique_tasks <- unique(train_data_model$Task_ID)
      ini_mix <- tibble(Task_ID = unique_tasks)
      for (k_idx in 1:n_clusters) {
        ini_mix[[paste0("K", k_idx)]] <- ifelse(seq_along(unique_tasks) == k_idx, 1, 0)
      }
    } else {
      ini_mix <- ini_mixture(data = train_data_model, k = n_clusters)
    }

    cluster_mapping <- ini_mix %>%
      pivot_longer(cols = starts_with("K"), names_to = "Cluster_ID", values_to = "Probability") %>%
      group_by(Task_ID) %>% slice_max(Probability, n = 1, with_ties = FALSE) %>%
      ungroup() %>% select(Task_ID, Cluster_ID)

    clusters_ids <- paste0("K", 1:n_clusters)
    train_data_model <- train_data_model %>%
      select(-any_of("Cluster_ID")) %>% left_join(cluster_mapping, by = "Task_ID")

    for (k_id in clusters_ids) {
      for (o_id in 1:n_out) {
        sub_data_agg <- train_data_model %>%
          dplyr::mutate(Input = round(Input, 6)) %>%
          dplyr::filter(Cluster_ID == k_id, Output_ID == as.character(o_id)) %>%
          dplyr::group_by(Input) %>%
          dplyr::summarise(Output = mean(Output, na.rm = TRUE), .groups = "drop") %>%
          dplyr::mutate(Output_ID = as.factor("1"), Input_ID = "1")

        mean_emp <- mean(sub_data_agg$Output)
        mean_vec <- rep(mean_emp, nrow(sub_data_agg))

        best_ll <- +Inf
        best_hp_gp <- NULL
        for(seed_retry in 1:10) {
          tryCatch({
            set.seed(seed_retry * 1000 + as.numeric(o_id))
            hp_tmp <- train_gp(data = sub_data_agg, kern = "SE", prior_mean = mean_emp, ini_hp = NULL)

            sub_data_agg_format_logL <- data.frame(
              Input_1 = as.numeric(sub_data_agg$Input),
              Output_ID = as.character(sub_data_agg$Output_ID),
              Output = as.numeric(sub_data_agg$Output),
              stringsAsFactors = FALSE
            )
            sub_data_agg_format_logL$Reference <- paste0("o", sub_data_agg_format_logL$Output_ID, ";", sub_data_agg_format_logL$Input_1)

            ll_val <- tryCatch({
              logL_GP_outside_package(hp = hp_tmp, db = sub_data_agg_format_logL,
                mean = mean_vec, kern = "SE", post_cov = 0, pen_diag = 1e-10,
                hp_col_names = c("se_variance", "se_lengthscale"))
            }, error = function(e) return(+Inf))

            if(ll_val < best_ll) { best_ll <- ll_val; best_hp_gp <- hp_tmp }
          }, error = function(e) {
            cat(paste0("  [train_gp ERR] k=", k_id, " o=", o_id,
                       " seed=", seed_retry, ": ", e$message, "\n"))
          })
        }

        if(is.null(best_hp_gp)) {
          cat(paste0("  [ATTENTION] Aucun train_gp n'a réussi pour k=", k_id, " o=", o_id,
                     ". Utilisation de valeurs par défaut.\n"))
          best_hp_gp <- list(se_lengthscale = 0, se_variance = 0, noise = -3)
        }
        l_val <- best_hp_gp$se_lengthscale
        v_val <- best_hp_gp$se_variance
        if(l_val > max_lengthscale_observed) max_lengthscale_observed <- l_val

        hp_k_extracted_list[[length(hp_k_extracted_list) + 1]] <- tibble(
          Cluster_ID = k_id, Output_ID = as.factor(o_id),
          l_t = l_val, S_t = v_val, noise = best_hp_gp$noise, prior_mean = mean_emp
        )
      }
    }

    ini_hp_k_raw <- bind_rows(hp_k_extracted_list)
    ini_hp_k_raw$noise <- -3
    ini_hp_k <- ini_hp_k_raw %>% dplyr::mutate(l_u_t = max_lengthscale_observed) %>% dplyr::select(-prior_mean)

    # --- D. Initialisation Tâches ---
    hp_cluster_1 <- ini_hp_k %>%
      dplyr::filter(Cluster_ID == clusters_ids[[1]]) %>%
      dplyr::mutate(Output_Num = as.numeric(as.character(Output_ID))) %>%
      dplyr::arrange(Output_Num)

    ini_hp_t <- hp(
      kern = convolution_kernel,
      list_task_ID = unique(as.character(train_data_model$Task_ID)),
      list_output_ID = unique(as.character(train_data_model$Output_ID)),
      shared_hp_outputs = FALSE, shared_hp_tasks = TRUE, noise = TRUE,
      hp_config = tibble(
        output_id = hp_cluster_1$Output_Num,
        lt_min = hp_cluster_1$l_t, lt_max = hp_cluster_1$l_t,
        St_min = rep(0, n_out), St_max = rep(0, n_out),
        lu_min = hp_cluster_1$l_u_t, lu_max = hp_cluster_1$l_u_t,
        noise_min = -3, noise_max = -3
      )
    )

    # --- E. Entraînement MagmaClust ---
    if (n_unique_tasks > n_clusters) {
      ini_mix <- ini_mixture(data = train_data_model, k = n_clusters)
    }
    prior_means_vec <- ini_hp_k_raw %>% arrange(Cluster_ID, Output_ID) %>% pull(prior_mean)

    t_model_start <- Sys.time()
    trained_model <- train_magmaclust(
      data = train_data_model %>% dplyr::select(-Cluster_ID),
      nb_cluster = n_clusters,
      prior_mean_k = prior_means_vec,
      ini_mixture = ini_mix,
      ini_hp_k = ini_hp_k %>% dplyr::select(-noise),
      ini_hp_t = ini_hp_t,
      kern_k = convolution_kernel, kern_t = convolution_kernel,
      shared_hp_clusts = FALSE, shared_hp_tasks = TRUE,
      pen_diag = 1e-8, cv_threshold = 0.001
    )
    duration_train <- as.numeric(difftime(Sys.time(), t_model_start, units = "secs"))
    cat(paste0("  Training: ", round(duration_train, 1), "s\n"))
    flush.console()

    # --- F. Hyperpostérieur global ---
    all_pred_inputs_global <- list()
    all_test_inputs_global <- list()
    for (tid in test_task_ids_max) {
      all_pred_inputs_global[[tid]] <- pred_tasks_data_all[[tid]] %>% dplyr::select(Output_ID, Input_ID, Input)
      all_test_inputs_global[[tid]] <- test_tasks_data_all[[tid]] %>% dplyr::select(Input, Input_ID, Output_ID)
    }
    grid_hyperpost_global <- bind_rows(bind_rows(all_pred_inputs_global), bind_rows(all_test_inputs_global)) %>%
      dplyr::distinct() %>% dplyr::arrange(Output_ID, Input_ID, Input)

    t_hp_start <- Sys.time()
    hyperpost_global <- hyperposterior_clust(
      trained_model = trained_model, data = trained_model$ini_args$data,
      mixture = trained_model$hyperpost$mixture,
      hp_k = trained_model$hp_k, hp_t = trained_model$hp_t,
      grid_inputs = grid_hyperpost_global,
      kern_k = trained_model$kern_k, kern_t = trained_model$kern_t,
      prior_mean_k = trained_model$prior_mean_k
    )
    duration_hyperpost_global <- as.numeric(difftime(Sys.time(), t_hp_start, units = "secs"))
    cat(paste0("  Hyperpost: ", round(duration_hyperpost_global, 1), "s\n"))
    flush.console()

    # --- G. Boucle sur les subsets ---
    for (n_pred in n_pred_subsets) {
      if (n_pred > n_pred_max_actual) {
        cat(paste0("    n_pred=", n_pred, " SKIP\n")); next
      }
      cat(paste0("    n_pred=", n_pred, "... "))

      task_folder_name <- paste0("train_", n_train, "_pred_", n_pred)
      dir_datasets <- file.path(base_dir, "Datasets", paste0("n_out_", n_out), task_folder_name)
      dir_models_momt <- file.path(base_dir, "Models_MOMT", paste0("n_out_", n_out), task_folder_name)
      dir_predictions_momt <- file.path(base_dir, "Predictions_MOMT", paste0("n_out_", n_out), task_folder_name)

      test_task_ids <- test_task_ids_max[1:n_pred]
      pred_tasks_data <- pred_tasks_data_all[test_task_ids]
      test_tasks_data <- test_tasks_data_all[test_task_ids]

      datasets_to_save <- list(
        data_obs = data_obs, train_data = train_data_model,
        pred_tasks_data = pred_tasks_data, test_tasks_data = test_tasks_data,
        train_task_ids = train_task_ids, test_task_ids = test_task_ids,
        configs = configs, seed = current_seed, n_train = n_train, n_pred = n_pred
      )
      # Sauvegarde atomique : écriture .tmp puis renommage
      tmp_ds <- file.path(dir_datasets, paste0("datasets_iter_", iter, ".rds.tmp"))
      saveRDS(datasets_to_save, tmp_ds)
      file.rename(tmp_ds, file.path(dir_datasets, paste0("datasets_iter_", iter, ".rds")))

      model_to_save <- list(
        trained_model = trained_model, t_training = duration_train,
        t_hyperpost_global = duration_hyperpost_global,
        seed = current_seed, n_train = n_train
      )
      tmp_mod <- file.path(dir_models_momt, paste0("model_iter_", iter, ".rds.tmp"))
      saveRDS(model_to_save, tmp_mod)
      file.rename(tmp_mod, file.path(dir_models_momt, paste0("model_iter_", iter, ".rds")))

      all_predictions <- list()
      t_pred_total <- 0
      for (test_task_id in test_task_ids) {
        test_grid_inputs <- test_tasks_data[[test_task_id]] %>% dplyr::select(Input, Input_ID, Output_ID) %>% dplyr::distinct()

        t_pred_start <- Sys.time()
        pred_res <- pred_magmaclust(
          data = pred_tasks_data[[test_task_id]] %>% dplyr::select(-c(Cluster_ID)) %>% dplyr::mutate(Output_ID = as.factor(Output_ID)),
          trained_model = trained_model, grid_inputs = test_grid_inputs,
          kern = convolution_kernel, hyperpost = hyperpost_global,
          get_full_cov = TRUE, plot = FALSE
        )
        t_pred_task <- as.numeric(difftime(Sys.time(), t_pred_start, units = "secs"))
        t_pred_total <- t_pred_total + t_pred_task

        all_predictions[[test_task_id]] <- list(
          prediction = pred_res, truth = test_tasks_data[[test_task_id]],
          inputs_pred = pred_tasks_data[[test_task_id]], t_pred = t_pred_task
        )
      }

      predictions_to_save <- list(
        predictions = all_predictions, t_hyperpost_global = duration_hyperpost_global,
        t_pred_total = t_pred_total, t_training = duration_train,
        seed = current_seed, test_task_ids = test_task_ids,
        n_train = n_train, n_pred = n_pred
      )
      tmp_pred <- file.path(dir_predictions_momt, paste0("predictions_iter_", iter, ".rds.tmp"))
      saveRDS(predictions_to_save, tmp_pred)
      file.rename(tmp_pred, file.path(dir_predictions_momt, paste0("predictions_iter_", iter, ".rds")))

      cat(paste0("OK (", round(t_pred_total, 1), "s)\n"))
      flush.console()
      if(exists("all_predictions")) rm(all_predictions); gc()
    }

    if(exists("trained_model")) rm(trained_model)
    if(exists("hyperpost_global")) rm(hyperpost_global)
    if(exists("pred_res")) rm(pred_res)
    gc()

  }, error = function(e) {
    cat(paste0("\n!!! ERREUR [n_out=", n_out, ", n_train=", n_train, ", iter=", iter, "] :\n"))
    # Afficher toute la chaîne d'erreurs
    cat("  Message: ", conditionMessage(e), "\n")
    if (!is.null(e$parent)) {
      cat("  Cause:   ", conditionMessage(e$parent), "\n")
      if (!is.null(e$parent$parent)) {
        cat("  Cause 2: ", conditionMessage(e$parent$parent), "\n")
      }
    }
    cat("  Call:    ", deparse(e$call), "\n\n")
    flush.console()
  })
}

# --- Fin ---
duration_total <- as.numeric(difftime(Sys.time(), t_global_start, units = "hours"))
cat(paste0("\n=== TERMINÉ n_out=", n_out, " n_train=", n_train, " ===\n"))
cat(paste0("Durée : ", round(duration_total, 2), "h\n"))

run_info$ended_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_info$duration_hours <- round(duration_total, 2)
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)
