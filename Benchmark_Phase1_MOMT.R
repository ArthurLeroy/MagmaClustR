# ==========================================================================
# Benchmark_Phase1_MOMT.R
# Phase 1 : Test paramétrique — MOMT (MagmaClust multi-output multi-task)
#
# Usage :
#   Rscript Benchmark_Phase1_MOMT.R --param=probleme --config=default --seed=1
#
# Paramètres possibles :
#   --param   : probleme, hp_clusters, hp_tasks, inputs, n_out, n_train
#   --config  : default, variation
#   --seed    : 1..10
#
# Ce script génère les données, entraîne MOMT, et produit les prédictions.
# Les scripts MO et MT réutilisent les données générées ici.
# ==========================================================================

# --- 0. PARSING ARGUMENTS ---
args <- commandArgs(trailingOnly = TRUE)

parse_arg_str <- function(args, name, default = NULL) {
  pattern <- paste0("^--", name, "=")
  match <- grep(pattern, args, value = TRUE)
  if (length(match) > 0) return(sub(pattern, "", match[1]))
  if (!is.null(default)) return(default)
  stop(paste0("Argument requis manquant : --", name))
}

parse_arg_num <- function(args, name, default = NULL) {
  val <- parse_arg_str(args, name, if (!is.null(default)) as.character(default) else NULL)
  return(as.numeric(val))
}

PARAM  <- parse_arg_str(args, "param")
CONFIG <- parse_arg_str(args, "config")
SEED   <- parse_arg_num(args, "seed")

stopifnot(PARAM %in% c("probleme", "hp_clusters", "hp_tasks", "inputs", "n_out", "n_train"))
stopifnot(CONFIG %in% c("default", "variation"))
stopifnot(SEED >= 1 & SEED <= 10)

# --- 1. CONFIGURATION PAR DÉFAUT ---
EXPERIMENT          <- "interpolation"
SHARED_HP_CLUSTS    <- TRUE
SHARED_HP_TASKS     <- TRUE
SHARED_GRID_OUTPUTS <- FALSE
SHARED_GRID_CLUSTERS <- FALSE
N_OUT               <- 2
N_TRAIN             <- 20
N_PRED              <- 10
N_CLUSTERS          <- 3
N_POINTS            <- 30

# Surcharge selon le paramètre qui varie
if (CONFIG == "variation") {
  switch(PARAM,
    "probleme"     = { EXPERIMENT <- "forecasting" },
    "hp_clusters"  = { SHARED_HP_CLUSTS <- FALSE },
    "hp_tasks"     = { SHARED_HP_TASKS <- FALSE },
    "inputs"       = { SHARED_GRID_OUTPUTS <- TRUE; SHARED_GRID_CLUSTERS <- TRUE },
    "n_out"        = { N_OUT <- 8 },
    "n_train"      = { N_TRAIN <- 100 }
  )
}

n_tasks_total <- N_TRAIN + N_PRED

# --- 2. SETUP & LIBRARIES ---
username  <- Sys.getenv("USER")
pkg_dir   <- file.path("/scratch", username, "MagmaClustR")
base_path <- file.path("/scratch", username, "Phase1_experiments", PARAM, CONFIG)

setwd(pkg_dir)

library(mvtnorm)
library(tidyverse)
library(matrixStats)
library(MagmaClustR)

convolution_kernel <- MagmaClustR:::convolution_kernel

cat(paste0("=== Phase 1 MOMT | param=", PARAM, " config=", CONFIG, " seed=", SEED, " ===\n"))
cat(paste0("  experiment       = ", EXPERIMENT, "\n"))
cat(paste0("  shared_hp_clusts = ", SHARED_HP_CLUSTS, "\n"))
cat(paste0("  shared_hp_tasks  = ", SHARED_HP_TASKS, "\n"))
cat(paste0("  shared_grid_out  = ", SHARED_GRID_OUTPUTS, "\n"))
cat(paste0("  shared_grid_clus = ", SHARED_GRID_CLUSTERS, "\n"))
cat(paste0("  n_out            = ", N_OUT, "\n"))
cat(paste0("  n_train          = ", N_TRAIN, "\n"))
cat(paste0("  n_pred           = ", N_PRED, "\n"))
cat(paste0("  Date             : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste0("  Host             : ", Sys.info()["nodename"], "\n"))
cat(paste0("  PID              : ", Sys.getpid(), "\n\n"))

# Création des répertoires
dir_datasets         <- file.path(base_path, "Datasets")
dir_models_momt      <- file.path(base_path, "Models_MOMT")
dir_predictions_momt <- file.path(base_path, "Predictions_MOMT")
dir_run_info         <- file.path(base_path, "run_info")

for (d in c(dir_datasets, dir_models_momt, dir_predictions_momt, dir_run_info)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Sauvegarde métadonnées
run_info <- list(
  script = "Benchmark_Phase1_MOMT.R",
  param = PARAM, config = CONFIG, seed = SEED,
  experiment = EXPERIMENT, shared_hp_clusts = SHARED_HP_CLUSTS,
  shared_hp_tasks = SHARED_HP_TASKS, shared_grid_outputs = SHARED_GRID_OUTPUTS,
  shared_grid_clusters = SHARED_GRID_CLUSTERS,
  n_out = N_OUT, n_train = N_TRAIN, n_pred = N_PRED,
  started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  hostname = as.character(Sys.info()["nodename"]),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "N/A"),
  r_version = R.version.string, pid = Sys.getpid()
)
run_info_file <- file.path(dir_run_info,
  paste0("momt_seed", SEED, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"))
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)

t_global_start <- Sys.time()

# --- 3. FONCTIONS UTILITAIRES ---
generate_random_config <- function(n_out, n_clust = 3, shared_hp_tasks = TRUE) {
  hp_mp_list <- list()
  jitter_leng <- 0.2
  jitter_var  <- 0.5
  global_l0_min  <- log(1/1000)
  global_l0_max  <- log(1/100)
  global_lu0_max <- log(1/100)

  for (o in 1:n_out) {
    safe_max_l0 <- global_l0_max - jitter_leng - 0.05
    center_l0 <- runif(1, global_l0_min, safe_max_l0)
    center_S0 <- runif(1, log(1), log(50))

    valid_center_lu0 <- FALSE
    while (!valid_center_lu0) {
      center_lu0 <- runif(1, log(1/300), global_lu0_max)
      if (center_lu0 > center_l0 + jitter_leng) valid_center_lu0 <- TRUE
    }

    for (k in 1:n_clust) {
      l0_val <- runif(1, center_l0 - jitter_leng, center_l0 + jitter_leng)
      l0_val <- max(min(l0_val, global_l0_max), global_l0_min)

      S0_min_bound <- log(1)
      S0_max_bound <- log(20)
      S0_val <- runif(1, center_S0 - jitter_var, center_S0 + jitter_var)
      S0_val <- max(min(S0_val, S0_max_bound), S0_min_bound)

      lu0_min_bound <- log(1/300)
      valid_lu0 <- FALSE
      counter <- 0
      while (!valid_lu0) {
        lu0_val <- runif(1, center_lu0 - jitter_leng, center_lu0 + jitter_leng)
        lu0_val <- max(min(lu0_val, global_lu0_max), lu0_min_bound)
        if (lu0_val > l0_val) valid_lu0 <- TRUE
        counter <- counter + 1
        if (counter > 100) {
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

  # HP des tâches : si shared_hp_tasks=FALSE, on met un range pour varier
  task_jitter <- if (shared_hp_tasks) 0 else 0.2
  base_task_params <- list()
  for (o in 1:n_out) {
    lt_val <- runif(1, log(1/1000), log(1/100))
    St_val <- runif(1, log(0.4), log(1))
    valid_lut <- FALSE
    while (!valid_lut) {
      lut_val <- runif(1, log(1/300), log(1/100))
      if (lut_val > lt_val) valid_lut <- TRUE
    }
    base_task_params[[length(base_task_params) + 1]] <- tibble(
      output_id = o,
      lt_min = lt_val - task_jitter, lt_max = lt_val + task_jitter,
      St_min = St_val - task_jitter, St_max = St_val + task_jitter,
      noise_min = -3, noise_max = -3,
      lu_min = lut_val - task_jitter, lu_max = lut_val + task_jitter
    )
  }
  base_task_df <- bind_rows(base_task_params)

  hp_task_list_full <- list()
  for (k in 1:n_clust) hp_task_list_full[[k]] <- base_task_df
  hp_task_config <- bind_rows(hp_task_list_full)

  return(list(mp = hp_mp_config, task = hp_task_config))
}

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
  (-dmnorm(db$Output, mean, inv, log = T)) %>% return()
}

# --- 4. GÉNÉRATION DES DONNÉES ---
tryCatch({
  set.seed(SEED)

  configs   <- generate_random_config(N_OUT, n_clust = N_CLUSTERS, shared_hp_tasks = SHARED_HP_TASKS)
  my_priors <- rep(0, N_OUT * N_CLUSTERS)

  t_simu_start <- Sys.time()
  sim_results <- simulate_magmaclust_data_convol(
    nb_clusters           = N_CLUSTERS,
    total_tasks           = n_tasks_total,
    points_per_output_grid = rep(200, N_OUT),
    grid_ranges           = rep(list(c(-1, 1)), N_OUT),
    shared_hp_clusts      = SHARED_HP_CLUSTS,
    hp_config_mean_process = configs$mp,
    prior_means           = my_priors,
    shared_hp_tasks       = SHARED_HP_TASKS,
    hp_config_tasks       = configs$task,
    n_points_per_task_range = c(N_POINTS, N_POINTS),
    shared_hp_outputs     = FALSE,
    shared_grid_outputs   = SHARED_GRID_OUTPUTS,
    shared_grid_clusters  = SHARED_GRID_CLUSTERS
  )
  duration_simu <- as.numeric(difftime(Sys.time(), t_simu_start, units = "secs"))
  cat(paste0("  Simulation: ", round(duration_simu, 1), "s\n"))

  data_obs <- sim_results$simulated_data_df

  # --- 5. SPLIT TRAIN / TEST ---
  all_task_ids <- unique(data_obs$Task_ID)
  shuffled_task_ids <- as.character(sample(all_task_ids, size = length(all_task_ids)))
  train_task_ids    <- shuffled_task_ids[1:N_TRAIN]
  test_task_ids     <- shuffled_task_ids[(N_TRAIN + 1):(N_TRAIN + N_PRED)]

  pred_tasks_data_all <- list()
  test_tasks_data_all <- list()

  if (EXPERIMENT == "interpolation") {
    gap_start <- -0.5
    gap_end   <-  0.5

    # Données d'entraînement : gap retiré sur Output 1
    train_data_model <- data_obs %>%
      dplyr::filter(Task_ID %in% train_task_ids) %>%
      dplyr::filter(Output_ID != "1" | (Output_ID == "1" & !(Input >= gap_start & Input <= gap_end)))

    for (tid in test_task_ids) {
      task_full <- data_obs %>% dplyr::filter(Task_ID == tid)
      task_O1   <- task_full %>% dplyr::filter(Output_ID == "1") %>% dplyr::arrange(Input)

      test_tasks_data_all[[tid]] <- task_O1 %>% dplyr::filter(Input >= gap_start & Input <= gap_end)
      pred_tasks_data_all[[tid]] <- dplyr::bind_rows(
        task_full %>% dplyr::filter(Output_ID != "1"),
        task_O1 %>% dplyr::filter(!(Input >= gap_start & Input <= gap_end))
      ) %>% dplyr::arrange(Output_ID, Input)
    }

  } else {
    # Forecasting : données d'entraînement complètes
    train_data_model <- data_obs %>% dplyr::filter(Task_ID %in% train_task_ids)

    for (tid in test_task_ids) {
      task_full <- data_obs %>% dplyr::filter(Task_ID == tid)
      task_O1   <- task_full %>% dplyr::filter(Output_ID == "1") %>% dplyr::arrange(Input)

      test_tasks_data_all[[tid]] <- task_O1 %>% dplyr::filter(Input >= 0 & Input <= 1)
      pred_tasks_data_all[[tid]] <- dplyr::bind_rows(
        task_full %>% dplyr::filter(Output_ID != "1"),
        task_O1 %>% dplyr::filter(!(Input >= 0 & Input <= 1))
      ) %>% dplyr::arrange(Output_ID, Input)
    }
  }

  # --- 6. PRÉ-PROCESSING : INITIALISATION GP VANILLE ---
  hp_k_extracted_list <- list()
  max_lengthscale_observed <- -Inf
  n_unique_tasks <- length(unique(train_data_model$Task_ID))

  if (n_unique_tasks <= N_CLUSTERS) {
    unique_tasks <- unique(train_data_model$Task_ID)
    ini_mix <- tibble(Task_ID = unique_tasks)
    for (k_idx in 1:N_CLUSTERS) {
      ini_mix[[paste0("K", k_idx)]] <- ifelse(seq_along(unique_tasks) == k_idx, 1, 0)
    }
  } else {
    ini_mix <- ini_mixture(data = train_data_model, k = N_CLUSTERS)
  }

  cluster_mapping <- ini_mix %>%
    tidyr::pivot_longer(cols = starts_with("K"), names_to = "Cluster_ID", values_to = "Probability") %>%
    dplyr::group_by(Task_ID) %>%
    dplyr::slice_max(Probability, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(Task_ID, Cluster_ID)

  clusters_ids <- paste0("K", 1:N_CLUSTERS)
  train_data_model <- train_data_model %>%
    dplyr::select(-any_of("Cluster_ID")) %>%
    dplyr::left_join(cluster_mapping, by = "Task_ID")

  for (k_id in clusters_ids) {
    for (o_id in 1:N_OUT) {
      sub_data_agg <- train_data_model %>%
        dplyr::mutate(Input = round(Input, 6)) %>%
        dplyr::filter(Cluster_ID == k_id, Output_ID == as.character(o_id)) %>%
        dplyr::group_by(Input) %>%
        dplyr::summarise(Output = mean(Output, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(Output_ID = as.factor("1"), Input_ID = "1")

      mean_emp <- mean(sub_data_agg$Output)
      mean_vec <- rep(mean_emp, nrow(sub_data_agg))

      best_ll    <- +Inf
      best_hp_gp <- NULL
      for (seed_retry in 1:10) {
        tryCatch({
          set.seed(seed_retry * 1000 + as.numeric(o_id))
          hp_tmp <- train_gp(data = sub_data_agg, kern = "SE", prior_mean = mean_emp, ini_hp = NULL)

          sub_data_agg_format_logL <- data.frame(
            Input_1 = as.numeric(sub_data_agg$Input),
            Output_ID = as.character(sub_data_agg$Output_ID),
            Output = as.numeric(sub_data_agg$Output),
            stringsAsFactors = FALSE
          )
          sub_data_agg_format_logL$Reference <- paste0("o", sub_data_agg_format_logL$Output_ID,
                                                        ";", sub_data_agg_format_logL$Input_1)

          ll_val <- tryCatch({
            logL_GP_outside_package(hp = hp_tmp, db = sub_data_agg_format_logL,
              mean = mean_vec, kern = "SE", post_cov = 0, pen_diag = 1e-10,
              hp_col_names = c("se_variance", "se_lengthscale"))
          }, error = function(e) return(+Inf))

          if (ll_val < best_ll) { best_ll <- ll_val; best_hp_gp <- hp_tmp }
        }, error = function(e) {
          cat(paste0("  [train_gp ERR] k=", k_id, " o=", o_id,
                     " seed=", seed_retry, ": ", e$message, "\n"))
        })
      }

      if (is.null(best_hp_gp)) {
        cat(paste0("  [ATTENTION] Aucun train_gp n'a réussi pour k=", k_id, " o=", o_id,
                   ". Valeurs par défaut.\n"))
        best_hp_gp <- list(se_lengthscale = 0, se_variance = 0, noise = -3)
      }
      l_val <- best_hp_gp$se_lengthscale
      v_val <- best_hp_gp$se_variance
      if (l_val > max_lengthscale_observed) max_lengthscale_observed <- l_val

      hp_k_extracted_list[[length(hp_k_extracted_list) + 1]] <- tibble(
        Cluster_ID = k_id, Output_ID = as.factor(o_id),
        l_t = l_val, S_t = v_val, noise = best_hp_gp$noise, prior_mean = mean_emp
      )
    }
  }

  ini_hp_k_raw <- bind_rows(hp_k_extracted_list)
  ini_hp_k_raw$noise <- -3
  ini_hp_k <- ini_hp_k_raw %>%
    dplyr::mutate(l_u_t = max_lengthscale_observed) %>%
    dplyr::select(-prior_mean)

  # --- 7. INITIALISATION TÂCHES ---
  hp_cluster_1 <- ini_hp_k %>%
    dplyr::filter(Cluster_ID == clusters_ids[[1]]) %>%
    dplyr::mutate(Output_Num = as.numeric(as.character(Output_ID))) %>%
    dplyr::arrange(Output_Num)

  ini_hp_t <- hp(
    kern = convolution_kernel,
    list_task_ID = unique(as.character(train_data_model$Task_ID)),
    list_output_ID = unique(as.character(train_data_model$Output_ID)),
    shared_hp_outputs = FALSE, shared_hp_tasks = SHARED_HP_TASKS, noise = TRUE,
    hp_config = tibble(
      output_id = hp_cluster_1$Output_Num,
      lt_min = hp_cluster_1$l_t, lt_max = hp_cluster_1$l_t,
      St_min = rep(0, N_OUT), St_max = rep(0, N_OUT),
      lu_min = hp_cluster_1$l_u_t, lu_max = hp_cluster_1$l_u_t,
      noise_min = -3, noise_max = -3
    )
  )

  # --- 8. ENTRAÎNEMENT MAGMACLUST ---
  if (n_unique_tasks > N_CLUSTERS) {
    ini_mix <- ini_mixture(data = train_data_model, k = N_CLUSTERS)
  }
  prior_means_vec <- ini_hp_k_raw %>%
    dplyr::arrange(Cluster_ID, Output_ID) %>%
    dplyr::pull(prior_mean)

  t_model_start <- Sys.time()
  trained_model <- train_magmaclust(
    data = train_data_model %>% dplyr::select(-Cluster_ID),
    nb_cluster     = N_CLUSTERS,
    prior_mean_k   = prior_means_vec,
    ini_mixture    = ini_mix,
    ini_hp_k       = ini_hp_k %>% dplyr::select(-noise),
    ini_hp_t       = ini_hp_t,
    kern_k         = convolution_kernel,
    kern_t         = convolution_kernel,
    shared_hp_clusts = SHARED_HP_CLUSTS,
    shared_hp_tasks  = SHARED_HP_TASKS,
    pen_diag       = 1e-8,
    cv_threshold   = 0.001
  )
  duration_train <- as.numeric(difftime(Sys.time(), t_model_start, units = "secs"))
  cat(paste0("  Training: ", round(duration_train, 1), "s\n"))

  # --- 9. HYPERPOSTÉRIEUR GLOBAL ---
  all_pred_inputs_global <- list()
  all_test_inputs_global <- list()
  for (tid in test_task_ids) {
    all_pred_inputs_global[[tid]] <- pred_tasks_data_all[[tid]] %>%
      dplyr::select(Output_ID, Input_ID, Input)
    all_test_inputs_global[[tid]] <- test_tasks_data_all[[tid]] %>%
      dplyr::select(Input, Input_ID, Output_ID)
  }
  grid_hyperpost_global <- dplyr::bind_rows(
    dplyr::bind_rows(all_pred_inputs_global),
    dplyr::bind_rows(all_test_inputs_global)
  ) %>% dplyr::distinct() %>% dplyr::arrange(Output_ID, Input_ID, Input)

  t_hp_start <- Sys.time()
  hyperpost_global <- hyperposterior_clust(
    trained_model = trained_model,
    data          = trained_model$ini_args$data,
    mixture       = trained_model$hyperpost$mixture,
    hp_k          = trained_model$hp_k,
    hp_t          = trained_model$hp_t,
    grid_inputs   = grid_hyperpost_global,
    kern_k        = trained_model$kern_k,
    kern_t        = trained_model$kern_t,
    prior_mean_k  = trained_model$m_k
  )
  duration_hyperpost <- as.numeric(difftime(Sys.time(), t_hp_start, units = "secs"))
  cat(paste0("  Hyperpost: ", round(duration_hyperpost, 1), "s\n"))

  # Sauvegarde du modèle
  model_to_save <- list(
    trained_model      = trained_model,
    t_training         = duration_train,
    t_hyperpost_global = duration_hyperpost,
    seed               = SEED,
    n_train            = N_TRAIN
  )
  saveRDS(model_to_save, file.path(dir_models_momt, paste0("model_seed_", SEED, ".rds")))

  # --- 10. PRÉDICTION + MÉTRIQUES ---
  all_predictions <- list()
  t_pred_total <- 0

  for (test_task_id in test_task_ids) {
    test_grid_inputs <- test_tasks_data_all[[test_task_id]] %>%
      dplyr::select(Input, Input_ID, Output_ID) %>%
      dplyr::distinct()

    t_pred_start <- Sys.time()
    pred_res <- pred_magmaclust(
      data           = pred_tasks_data_all[[test_task_id]] %>%
                         dplyr::select(-any_of("Cluster_ID")) %>%
                         dplyr::mutate(Output_ID = as.factor(Output_ID)),
      trained_model  = trained_model,
      grid_inputs    = test_grid_inputs,
      kern           = convolution_kernel,
      hyperpost      = hyperpost_global,
      get_full_cov   = TRUE,
      plot           = FALSE
    )
    t_pred_task <- as.numeric(difftime(Sys.time(), t_pred_start, units = "secs"))
    t_pred_total <- t_pred_total + t_pred_task

    all_predictions[[test_task_id]] <- list(
      prediction  = pred_res,
      truth       = test_tasks_data_all[[test_task_id]],
      inputs_pred = pred_tasks_data_all[[test_task_id]],
      t_pred      = t_pred_task
    )
  }

  # Sauvegarde des datasets (pour MO et MT)
  datasets_to_save <- list(
    data_obs        = data_obs,
    train_data      = train_data_model,
    pred_tasks_data = pred_tasks_data_all,
    test_tasks_data = test_tasks_data_all,
    train_task_ids  = train_task_ids,
    test_task_ids   = test_task_ids,
    configs         = configs,
    seed            = SEED,
    n_train         = N_TRAIN,
    n_pred          = N_PRED
  )
  saveRDS(datasets_to_save, file.path(dir_datasets, paste0("datasets_seed_", SEED, ".rds")))

  # Sauvegarde des prédictions
  predictions_to_save <- list(
    predictions        = all_predictions,
    t_hyperpost_global = duration_hyperpost,
    t_pred_total       = t_pred_total,
    t_training         = duration_train,
    seed               = SEED,
    test_task_ids      = test_task_ids,
    n_train            = N_TRAIN,
    n_pred             = N_PRED
  )
  saveRDS(predictions_to_save, file.path(dir_predictions_momt, paste0("predictions_seed_", SEED, ".rds")))

  cat(paste0("  Prédictions OK (", round(t_pred_total, 1), "s total)\n"))

  rm(trained_model, hyperpost_global, all_predictions)
  gc()

}, error = function(e) {
  cat(paste0("\n!!! ERREUR [param=", PARAM, " config=", CONFIG, " seed=", SEED, "] :\n"))
  cat("  Message: ", conditionMessage(e), "\n")
  if (!is.null(e$parent)) cat("  Cause: ", conditionMessage(e$parent), "\n")
  cat("  Call: ", deparse(e$call), "\n\n")
})

# --- Fin ---
duration_total <- as.numeric(difftime(Sys.time(), t_global_start, units = "mins"))
cat(paste0("\n=== MOMT TERMINÉ [param=", PARAM, " config=", CONFIG, " seed=", SEED, "] ===\n"))
cat(paste0("Durée : ", round(duration_total, 2), " min\n"))

run_info$ended_at       <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_info$duration_mins  <- round(duration_total, 2)
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)
