# ==========================================================================
# Benchmark_NeurIPS_MOMT.R
# NeurIPS : Benchmark MOMT (MagmaClust multi-output multi-task)
#
# Usage :
#   Rscript Benchmark_NeurIPS_MOMT.R --n_out=2 --n_train=30 --n_pred=1 \
#     --problem=interpolation --seed=1
#
# Ce script :
#   1. Génère les données simulées
#   2. Scale les données (Output / max(|Output|) par Output_ID × Cluster_ID)
#   3. Pré-processing : clustering initial sur données brutes, HP vanille sur scalées
#   4. Entraîne MOMT (convolution_kernel, shared_hp_clusts=TRUE, shared_hp_tasks=TRUE)
#   5. Calcule l'hyperpostérieur global (une seule fois pour toutes les tâches pred)
#   6. Prédit tâche par tâche, sauvegarde prédictions + datasets
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

N_OUT     <- parse_arg_num(args, "n_out")
N_TRAIN   <- parse_arg_num(args, "n_train")
N_PRED    <- parse_arg_num(args, "n_pred")
PROBLEM   <- parse_arg_str(args, "problem")
SEED      <- parse_arg_num(args, "seed")

stopifnot(N_OUT %in% c(2, 4, 8))
stopifnot(N_TRAIN %in% c(15, 30, 100))
stopifnot(N_PRED %in% c(1, 10, 100))
stopifnot(PROBLEM %in% c("interpolation", "forecasting"))
stopifnot(SEED >= 1 & SEED <= 5)

# --- Paramètres fixes ---
N_CLUSTERS     <- 3
N_POINTS_FIXED <- 30

n_tasks_total <- N_TRAIN + N_PRED

# --- 1. SETUP & LIBRARIES ---
username  <- Sys.getenv("USER")
pkg_dir   <- file.path("/scratch", username, "MagmaClustR")
base_path <- file.path("/scratch", username, "NeurIPS_experiments",
                        paste0("out", N_OUT, "_train", N_TRAIN, "_pred", N_PRED),
                        PROBLEM)

setwd(pkg_dir)

library(mvtnorm)
library(tidyverse)
library(matrixStats)
library(MagmaClustR)

convolution_kernel <- MagmaClustR:::convolution_kernel

# Charger les fonctions utilitaires
source(file.path("/scratch", username, "NeurIPS_experiments", "utils", "generate_random_config.R"))

config_label <- paste0("out", N_OUT, "_train", N_TRAIN, "_pred", N_PRED)

cat(paste0("=== NeurIPS MOMT | ", config_label, " | ", PROBLEM, " | seed=", SEED, " ===\n"))
cat(paste0("  n_out    = ", N_OUT, "\n"))
cat(paste0("  n_train  = ", N_TRAIN, "\n"))
cat(paste0("  n_pred   = ", N_PRED, "\n"))
cat(paste0("  problem  = ", PROBLEM, "\n"))
cat(paste0("  Date     : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste0("  Host     : ", Sys.info()["nodename"], "\n"))
cat(paste0("  PID      : ", Sys.getpid(), "\n\n"))

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
  script = "Benchmark_NeurIPS_MOMT.R",
  n_out = N_OUT, n_train = N_TRAIN, n_pred = N_PRED,
  problem = PROBLEM, seed = SEED,
  n_clusters = N_CLUSTERS, n_points = N_POINTS_FIXED,
  started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  hostname = as.character(Sys.info()["nodename"]),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "N/A"),
  r_version = R.version.string, pid = Sys.getpid()
)
run_info_file <- file.path(dir_run_info,
  paste0("momt_seed", SEED, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"))
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)

t_global_start <- Sys.time()

# --- 2. GÉNÉRATION DES DONNÉES ---
tryCatch({
  set.seed(SEED)

  configs   <- generate_random_config(N_OUT, n_clust = N_CLUSTERS)
  my_priors <- rep(0, N_OUT * N_CLUSTERS)

  t_simu_start <- Sys.time()
  sim_results <- simulate_magmaclust_data_convol(
    nb_clusters            = N_CLUSTERS,
    total_tasks            = n_tasks_total,
    points_per_output_grid = rep(200, N_OUT),
    grid_ranges            = rep(list(c(-1, 1)), N_OUT),
    shared_hp_clusts       = FALSE,
    hp_config_mean_process = configs$mp,
    prior_means            = my_priors,
    shared_hp_tasks        = TRUE,
    hp_config_tasks        = configs$task,
    n_points_per_task_range = c(N_POINTS_FIXED, N_POINTS_FIXED),
    shared_hp_outputs      = FALSE,
    shared_grid_outputs    = FALSE,
    shared_grid_clusters   = TRUE
  )
  duration_simu <- as.numeric(difftime(Sys.time(), t_simu_start, units = "secs"))
  cat(paste0("  Simulation: ", round(duration_simu, 1), "s\n"))

  # Données brutes (pour ini_mixture + cluster lookup pour unscaling)
  data_raw <- sim_results$simulated_data_df

  # --- 3. SCALING ---
  scaling_result <- scale_data(data_raw)
  data_obs       <- scaling_result$data_scaled
  scale_factors  <- scaling_result$scale_factors

  # --- 4. SPLIT TRAIN / TEST ---
  all_task_ids      <- unique(data_obs$Task_ID)
  shuffled_task_ids <- as.character(sample(all_task_ids, size = length(all_task_ids)))
  train_task_ids    <- shuffled_task_ids[1:N_TRAIN]
  test_task_ids     <- shuffled_task_ids[(N_TRAIN + 1):(N_TRAIN + N_PRED)]

  pred_tasks_data_all <- list()
  test_tasks_data_all <- list()

  if (PROBLEM == "interpolation") {
    gap_start <- -0.5
    gap_end   <-  0.5

    # Entraînement : gap retiré sur Output 1 pour TOUTES les tâches train
    train_data_model <- data_obs %>%
      dplyr::filter(Task_ID %in% train_task_ids) %>%
      dplyr::filter(Output_ID != "1" |
                    (Output_ID == "1" & !(Input >= gap_start & Input <= gap_end)))

    for (tid in test_task_ids) {
      task_full <- data_obs %>% dplyr::filter(Task_ID == tid)
      task_O1   <- task_full %>% dplyr::filter(Output_ID == "1") %>% dplyr::arrange(Input)

      # Points test : ceux dans le gap sur Output 1
      test_tasks_data_all[[tid]] <- task_O1 %>%
        dplyr::filter(Input >= gap_start & Input <= gap_end)

      # Contexte : tous les autres outputs + Output 1 hors gap
      pred_tasks_data_all[[tid]] <- dplyr::bind_rows(
        task_full %>% dplyr::filter(Output_ID != "1"),
        task_O1 %>% dplyr::filter(!(Input >= gap_start & Input <= gap_end))
      ) %>% dplyr::arrange(Output_ID, Input)
    }

  } else {
    # Forecasting : entraînement complet
    train_data_model <- data_obs %>% dplyr::filter(Task_ID %in% train_task_ids)

    for (tid in test_task_ids) {
      task_full <- data_obs %>% dplyr::filter(Task_ID == tid)
      task_O1   <- task_full %>% dplyr::filter(Output_ID == "1") %>% dplyr::arrange(Input)

      # Points test : Output 1, intervalle [0, 1]
      test_tasks_data_all[[tid]] <- task_O1 %>%
        dplyr::filter(Input >= 0 & Input <= 1)

      # Contexte : tous les autres outputs + Output 1 hors [0, 1]
      pred_tasks_data_all[[tid]] <- dplyr::bind_rows(
        task_full %>% dplyr::filter(Output_ID != "1"),
        task_O1 %>% dplyr::filter(!(Input >= 0 & Input <= 1))
      ) %>% dplyr::arrange(Output_ID, Input)
    }
  }

  # --- 5. PRÉ-PROCESSING : CLUSTERING INITIAL SUR DONNÉES BRUTES ---
  n_unique_tasks <- length(unique(train_data_model$Task_ID))

  # ini_mixture sur données BRUTES (non scalées) pour un clustering plus net
  if (n_unique_tasks <= N_CLUSTERS) {
    unique_tasks <- unique(train_data_model$Task_ID)
    ini_mix <- tibble(Task_ID = unique_tasks)
    for (k_idx in 1:N_CLUSTERS) {
      ini_mix[[paste0("K", k_idx)]] <- ifelse(seq_along(unique_tasks) == k_idx, 1, 0)
    }
    cat("  (init manuelle car n_train <= n_clusters)\n")
  } else {
    # Données brutes pour ini_mixture (les mêmes tâches train, mais non scalées)
    data_raw_train <- data_raw %>% dplyr::filter(Task_ID %in% train_task_ids)
    ini_mix <- ini_mixture(data = data_raw_train, k = N_CLUSTERS)
  }

  cluster_mapping <- ini_mix %>%
    tidyr::pivot_longer(cols = starts_with("K"), names_to = "Cluster_ID",
                        values_to = "Probability") %>%
    dplyr::group_by(Task_ID) %>%
    dplyr::slice_max(Probability, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(Task_ID, Cluster_ID)

  clusters_ids <- paste0("K", 1:N_CLUSTERS)
  train_data_model <- train_data_model %>%
    dplyr::select(-any_of("Cluster_ID")) %>%
    dplyr::left_join(cluster_mapping, by = "Task_ID")

  # --- 6. INITIALISATION HP VANILLE (sur données scalées) ---
  hp_k_extracted_list <- list()
  max_lengthscale_observed <- -Inf

  for (k_id in clusters_ids) {
    for (o_id in 1:N_OUT) {
      sub_data_agg <- train_data_model %>%
        dplyr::mutate(Input = round(Input, 6)) %>%
        dplyr::filter(Cluster_ID == k_id, Output_ID == as.character(o_id)) %>%
        dplyr::group_by(Input) %>%
        dplyr::summarise(Output = mean(Output, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(Output_ID = as.factor("1"))

      mean_emp <- mean(sub_data_agg$Output)
      mean_vec <- rep(mean_emp, nrow(sub_data_agg))

      best_ll    <- +Inf
      best_hp_gp <- NULL

      for (seed_retry in 1:10) {
        tryCatch({
          set.seed(seed_retry * 1000 + as.numeric(o_id))
          hp_tmp <- suppressWarnings(suppressMessages(
            train_gp(data = sub_data_agg, kern = "SE",
                     prior_mean = mean_emp, ini_hp = NULL)
          ))

          sub_data_agg_format_logL <- data.frame(
            Input_1   = as.numeric(sub_data_agg$Input),
            Output_ID = as.character(sub_data_agg$Output_ID),
            Output    = as.numeric(sub_data_agg$Output),
            stringsAsFactors = FALSE
          )
          sub_data_agg_format_logL$Reference <- paste0(
            "o", sub_data_agg_format_logL$Output_ID, ";",
            sub_data_agg_format_logL$Input_1
          )

          ll_val <- tryCatch({
            logL_GP_outside_package(
              hp = hp_tmp, db = sub_data_agg_format_logL,
              mean = mean_vec, kern = "SE", post_cov = 0,
              pen_diag = 1e-10, hp_col_names = c("se_variance", "se_lengthscale")
            )
          }, error = function(e) return(+Inf))

          if (ll_val < best_ll) {
            best_ll    <- ll_val
            best_hp_gp <- hp_tmp
          }
        }, error = function(e) {
          cat(paste0("  [train_gp ERR] k=", k_id, " o=", o_id,
                     " retry=", seed_retry, ": ", e$message, "\n"))
        })
      }

      if (is.null(best_hp_gp)) {
        cat(paste0("  [WARN] Aucun train_gp OK pour k=", k_id, " o=", o_id, "\n"))
        best_hp_gp <- list(se_lengthscale = 0, se_variance = 0, noise = -3)
      }

      l_val <- best_hp_gp$se_lengthscale
      v_val <- best_hp_gp$se_variance
      if (l_val > max_lengthscale_observed) max_lengthscale_observed <- l_val

      hp_k_extracted_list[[length(hp_k_extracted_list) + 1]] <- tibble(
        Cluster_ID = k_id, Output_ID = as.factor(o_id),
        l_t = l_val, S_t = v_val, noise = best_hp_gp$noise,
        prior_mean = mean_emp
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
    shared_hp_outputs = FALSE,
    shared_hp_tasks = TRUE,
    noise = TRUE,
    hp_config = tibble(
      output_id = hp_cluster_1$Output_Num,
      lt_min = hp_cluster_1$l_t, lt_max = hp_cluster_1$l_t,
      St_min = rep(0, N_OUT), St_max = rep(0, N_OUT),
      lu_min = hp_cluster_1$l_u_t, lu_max = hp_cluster_1$l_u_t,
      noise_min = -3, noise_max = -3
    )
  )

  # --- 8. ENTRAÎNEMENT MAGMACLUST ---
  # Réutiliser ini_mix (calculé sur données brutes en étape 5)
  prior_means_vec <- ini_hp_k_raw %>%
    dplyr::arrange(Cluster_ID, Output_ID) %>%
    dplyr::pull(prior_mean)

  t_model_start <- Sys.time()
  trained_model <- train_magmaclust(
    data             = train_data_model %>% dplyr::select(-Cluster_ID),
    nb_cluster       = N_CLUSTERS,
    prior_mean_k     = prior_means_vec,
    ini_mixture      = ini_mix,
    ini_hp_k         = ini_hp_k %>% dplyr::select(-noise),
    ini_hp_t         = ini_hp_t,
    kern_k           = convolution_kernel,
    kern_t           = convolution_kernel,
    shared_hp_clusts = TRUE,
    shared_hp_tasks  = TRUE,
    pen_diag         = 1e-8,
    cv_threshold     = 0.001
  )
  duration_train <- as.numeric(difftime(Sys.time(), t_model_start, units = "secs"))
  cat(paste0("  Training: ", round(duration_train, 1), "s\n"))

  # --- 9. HYPERPOSTÉRIEUR GLOBAL (une seule fois) ---
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
    n_train            = N_TRAIN,
    ini_hp_k_raw       = ini_hp_k_raw,
    ini_hp_t           = ini_hp_t
  )
  saveRDS(model_to_save, file.path(dir_models_momt, paste0("model_seed_", SEED, ".rds")))

  # --- 10. PRÉDICTION TÂCHE PAR TÂCHE ---
  all_predictions <- list()
  t_pred_total <- 0

  for (test_task_id in test_task_ids) {
    test_grid_inputs <- test_tasks_data_all[[test_task_id]] %>%
      dplyr::select(Input, Input_ID, Output_ID) %>%
      dplyr::distinct()

    t_pred_start <- Sys.time()
    pred_res <- pred_magmaclust(
      data          = pred_tasks_data_all[[test_task_id]] %>%
                        dplyr::select(-any_of("Cluster_ID")) %>%
                        dplyr::mutate(Output_ID = as.factor(Output_ID)),
      trained_model = trained_model,
      grid_inputs   = test_grid_inputs,
      kern          = convolution_kernel,
      hyperpost     = hyperpost_global,
      get_full_cov  = TRUE,
      plot          = FALSE
    )
    t_pred_task <- as.numeric(difftime(Sys.time(), t_pred_start, units = "secs"))
    t_pred_total <- t_pred_total + t_pred_task

    # Récupérer le cluster de la tâche pour le dé-scaling au moment des métriques
    task_cluster <- data_raw %>%
      dplyr::filter(Task_ID == test_task_id) %>%
      dplyr::pull(Cluster_ID) %>%
      unique()

    all_predictions[[test_task_id]] <- list(
      prediction    = pred_res,
      truth         = test_tasks_data_all[[test_task_id]],
      inputs_pred   = pred_tasks_data_all[[test_task_id]],
      t_pred        = t_pred_task,
      cluster_id    = task_cluster[1]
    )
  }

  # --- 11. SAUVEGARDE DATASETS (pour MO et MT) ---
  datasets_to_save <- list(
    data_obs          = data_obs,
    data_raw          = data_raw,
    scale_factors     = scale_factors,
    train_data        = train_data_model,
    pred_tasks_data   = pred_tasks_data_all,
    test_tasks_data   = test_tasks_data_all,
    train_task_ids    = train_task_ids,
    test_task_ids     = test_task_ids,
    configs           = configs,
    seed              = SEED,
    n_out             = N_OUT,
    n_train           = N_TRAIN,
    n_pred            = N_PRED,
    problem           = PROBLEM,
    cluster_mapping   = cluster_mapping
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
  saveRDS(predictions_to_save,
          file.path(dir_predictions_momt, paste0("predictions_seed_", SEED, ".rds")))

  cat(paste0("  Prédictions OK (", round(t_pred_total, 1), "s total)\n"))

  rm(trained_model, hyperpost_global, all_predictions)
  gc()

}, error = function(e) {
  cat(paste0("\n!!! ERREUR [MOMT out=", N_OUT, " train=", N_TRAIN,
             " pred=", N_PRED, " ", PROBLEM, " seed=", SEED, "] :\n"))
  cat("  Message: ", conditionMessage(e), "\n")
  if (!is.null(e$parent)) cat("  Cause: ", conditionMessage(e$parent), "\n")
  cat("  Call: ", deparse(e$call), "\n\n")
})

# --- Fin ---
duration_total <- as.numeric(difftime(Sys.time(), t_global_start, units = "mins"))
cat(paste0("\n=== MOMT TERMINÉ [out=", N_OUT, " train=", N_TRAIN,
           " pred=", N_PRED, " ", PROBLEM, " seed=", SEED, "] ===\n"))
cat(paste0("Durée : ", round(duration_total, 2), " min\n"))

run_info$ended_at      <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_info$duration_mins <- round(duration_total, 2)
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)
