# ==========================================================================
# Benchmark_XP_1_MT_cluster.R
# Version adaptée pour le cluster IMT (SLURM) — PARAMÉTRÉ
#
# Usage :
#   Rscript Benchmark_XP_1_MT_cluster.R --n_out=3 --n_train=15
#
# MagmaClust single-output (Output_ID == "1" uniquement).
# Charge les données/modèles générés par le script MOMT.
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
username  <- Sys.getenv("USER")
pkg_dir   <- file.path("/scratch", username, "MagmaClustR")
base_path <- file.path("/scratch", username, "NeurIPS_experiments", "Experience_1")

setwd(pkg_dir)

library(tidyverse)
library(readr)
library(devtools)
devtools::load_all(pkg_dir)

n_out <- N_OUT_TARGET
n_train <- N_TRAIN_TARGET
n_iterations <- 100
n_clusters <- 3
n_pred_max <- 1000
n_pred_subsets <- c(1, 10, 100, 1000)

cat(paste0("=== Benchmark XP1 MT (Single-Output) ===\n"))
cat(paste0("  n_out   = ", n_out, "\n"))
cat(paste0("  n_train = ", n_train, "\n"))
cat(paste0("  Date    : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste0("  Host    : ", Sys.info()["nodename"], "\n"))
cat(paste0("  PID     : ", Sys.getpid(), "\n\n"))

# Sauvegarde métadonnées
run_info <- list(
  script = "Benchmark_XP_1_MT_cluster.R",
  n_out = n_out, n_train = n_train, n_iterations = n_iterations,
  started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  hostname = as.character(Sys.info()["nodename"]),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "N/A"),
  r_version = R.version.string, pid = Sys.getpid()
)
run_info_dir <- file.path(base_path, "run_info")
if (!dir.exists(run_info_dir)) dir.create(run_info_dir, recursive = TRUE)
run_info_file <- file.path(run_info_dir,
  paste0("mt_nout", n_out, "_ntrain", n_train, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"))
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)

t_global_start <- Sys.time()

# --- 2. BOUCLE SUR LES ITÉRATIONS ---
for (iter in 1:n_iterations) {
  tryCatch({
    cat(paste0("[n_out=", n_out, " n_train=", n_train, "] MT Iter ", iter, "/", n_iterations,
               " [", format(Sys.time(), "%H:%M:%S"), "]... "))

    # dir_n_out <- file.path(base_path, paste0("n_out_", n_out))
    # dir_datasets <- file.path(dir_n_out, "Datasets", paste0("train_", n_train, "_pred_", n_pred_max))
    # dir_models_momt <- file.path(dir_n_out, "Models_MOMT", paste0("train_", n_train, "_pred_", n_pred_max))
    # dir_models_mt <- file.path(dir_n_out, "Models_MT", paste0("train_", n_train, "_pred_", n_pred_max))
    dir_datasets <- file.path(base_path, "Datasets", paste0("n_out_", n_out), paste0("train_", n_train, "_pred_", n_pred_max))
    dir_models_momt <- file.path(base_path, "Models_MOMT", paste0("n_out_", n_out), paste0("train_", n_train, "_pred_", n_pred_max))
    dir_models_mt <- file.path(base_path, "Models_MT", paste0("n_out_", n_out), paste0("train_", n_train, "_pred_", n_pred_max))
    if (!dir.exists(dir_models_mt)) dir.create(dir_models_mt, recursive = TRUE)

    file_name_datasets <- file.path(dir_datasets, paste0("datasets_iter_", iter, ".rds"))
    file_name_model_momt <- file.path(dir_models_momt, paste0("model_iter_", iter, ".rds"))

    if (!file.exists(file_name_datasets)) { cat("SKIP (pas de données)\n"); next }
    if (!file.exists(file_name_model_momt)) { cat("SKIP (pas de modèle MOMT)\n"); next }

    datasets <- readRDS(file_name_datasets)
    model_momt <- readRDS(file_name_model_momt)
    trained_model_momt <- model_momt$trained_model

    data_obs <- datasets$data_obs
    train_data <- datasets$train_data
    pred_tasks_data_all <- datasets$pred_tasks_data
    test_tasks_data_all <- datasets$test_tasks_data
    test_task_ids_max <- datasets$test_task_ids

    # Filtrage Output_ID == "1"
    train_data_o1 <- train_data %>%
      dplyr::filter(Output_ID == "1") %>%
      dplyr::mutate(Output_ID = as.factor(Output_ID)) %>%
      dplyr::select(-Cluster_ID)

    # Extraction HP initiaux
    ini_hp_k <- trained_model_momt$hp_k %>%
      dplyr::filter(Output_ID == "1") %>%
      dplyr::select(-c(l_u_t)) %>%
      dplyr::mutate(se_lengthscale = l_t, se_variance = S_t) %>%
      dplyr::select(-c(l_t, S_t))

    ini_hp_t <- trained_model_momt$hp_t %>%
      dplyr::filter(Output_ID == "1") %>%
      dplyr::filter(Task_ID %in% unique(train_data$Task_ID)) %>%
      dplyr::select(-c(l_u_t, Task_ID)) %>%
      dplyr::mutate(se_lengthscale = l_t, se_variance = S_t) %>%
      dplyr::select(-c(l_t, S_t)) %>%
      dplyr::slice(1)

    indices <- seq(from = 1, to = length(trained_model_momt$ini_args$prior_mean_k), by = n_out)
    prior_mean_k_o1 <- trained_model_momt$ini_args$prior_mean_k[indices]
    ini_mixture <- trained_model_momt$ini_args$ini_mixture

    # Entraînement
    t_train_start <- Sys.time()
    trained_model_mt <- train_magmaclust(
      data = train_data_o1, ini_mixture = ini_mixture,
      ini_hp_k = ini_hp_k, nb_cluster = n_clusters,
      kern_k = "SE", ini_hp_t = ini_hp_t, kern_t = "SE",
      shared_hp_tasks = TRUE, shared_hp_clusts = FALSE, prior_mean_k = prior_mean_k_o1, pen_diag = 1e-10
    )
    duration_train <- as.numeric(difftime(Sys.time(), t_train_start, units = "secs"))

    # Sauvegarde modèle MT
    model_mt_to_save <- list(
      trained_model = trained_model_mt, t_training = duration_train,
      seed = datasets$seed, n_train = n_train
    )
    saveRDS(model_mt_to_save, file.path(dir_models_mt, paste0("model_iter_", iter, ".rds")))

    # --- Hyperpostérieur global (1 seul calcul, comme MOMT) ---
    n_tasks_total <- length(test_task_ids_max)
    save_points <- sort(unique(pmin(n_pred_subsets, n_tasks_total)))
    n_tasks_needed <- max(save_points)

    # Construire la grille globale : tous les inputs (pred + test) de toutes les tâches
    all_pred_inputs_global <- list()
    all_test_inputs_global <- list()
    for (tid in test_task_ids_max[1:n_tasks_needed]) {
      pred_task_o1 <- pred_tasks_data_all[[tid]] %>% dplyr::filter(Output_ID == "1")
      test_task_o1 <- test_tasks_data_all[[tid]] %>% dplyr::filter(Output_ID == "1")
      all_pred_inputs_global[[tid]] <- pred_task_o1 %>% dplyr::select(Output_ID, Input_ID, Input)
      all_test_inputs_global[[tid]] <- test_task_o1 %>% dplyr::select(Input, Input_ID, Output_ID)
    }
    grid_hyperpost_global <- bind_rows(bind_rows(all_pred_inputs_global), bind_rows(all_test_inputs_global)) %>%
      dplyr::distinct() %>% dplyr::arrange(Output_ID, Input_ID, Input)

    t_hp_start <- Sys.time()
    hyperpost_global <- hyperposterior_clust(
      trained_model = trained_model_mt, data = trained_model_mt$ini_args$data,
      mixture = trained_model_mt$hyperpost$mixture,
      hp_k = trained_model_mt$hp_k, hp_t = trained_model_mt$hp_t,
      grid_inputs = grid_hyperpost_global,
      kern_k = trained_model_mt$kern_k, kern_t = trained_model_mt$kern_t,
      prior_mean_k = trained_model_mt$m_k
    )
    duration_hyperpost_global <- as.numeric(difftime(Sys.time(), t_hp_start, units = "secs"))
    cat(paste0("  Hyperpost: ", round(duration_hyperpost_global, 1), "s\n"))

    # --- Prédiction tâche par tâche (réutilisation de l'hyperpost global) ---
    all_predictions <- list()

    for (idx in 1:n_tasks_needed) {
      test_task_id <- test_task_ids_max[idx]

      pred_task_input_data <- pred_tasks_data_all[[test_task_id]] %>% dplyr::filter(Output_ID == "1")
      test_data_truth <- test_tasks_data_all[[test_task_id]] %>% dplyr::filter(Output_ID == "1")

      test_grid_inputs <- test_data_truth %>%
        dplyr::select(Input, Input_ID, Output_ID) %>%
        dplyr::distinct() %>% dplyr::mutate(Output_ID = as.factor(Output_ID))

      pred_data_task <- pred_task_input_data %>%
        dplyr::select(-Task_ID, -Cluster_ID) %>%
        dplyr::mutate(Output_ID = as.factor(Output_ID))

      t_pred_start <- Sys.time()
      pred_res <- pred_magmaclust(
        data = pred_data_task, trained_model = trained_model_mt,
        grid_inputs = test_grid_inputs, kern = "SE",
        hyperpost = hyperpost_global, get_full_cov = TRUE, plot = FALSE
      )
      t_pred_task <- as.numeric(difftime(Sys.time(), t_pred_start, units = "secs"))

      all_predictions[[test_task_id]] <- list(
        prediction = pred_res, truth = test_data_truth,
        inputs_pred = pred_task_input_data,
        t_pred = t_pred_task
      )

      # Sauvegarde immédiate dès qu'un subset est complet
      if (idx %in% save_points) {
        for (n_pred in n_pred_subsets[pmin(n_pred_subsets, n_tasks_total) == idx]) {
          test_task_ids <- test_task_ids_max[1:idx]

          dir_predictions_mt <- file.path(base_path, "Predictions_MT", paste0("n_out_", n_out), paste0("train_", n_train, "_pred_", n_pred))
          if (!dir.exists(dir_predictions_mt)) dir.create(dir_predictions_mt, recursive = TRUE)

          subset_predictions <- all_predictions[test_task_ids]
          t_pred_total <- sum(sapply(subset_predictions, function(x) x$t_pred))

          predictions_to_save <- list(
            predictions = subset_predictions, t_hyperpost_global = duration_hyperpost_global,
            t_pred_total = t_pred_total, t_training = duration_train,
            seed = datasets$seed, test_task_ids = test_task_ids,
            n_train = n_train, n_pred = n_pred
          )
          saveRDS(predictions_to_save, file.path(dir_predictions_mt, paste0("predictions_iter_", iter, ".rds")))
          cat(paste0("[saved pred_", n_pred, "] "))
        }
      }
    }

    cat(paste0("OK (Train: ", round(duration_train, 1), "s, Hyp: ", round(duration_hyperpost_global, 1), "s)\n"))

    if(exists("trained_model_mt")) rm(trained_model_mt)
    if(exists("hyperpost_global")) rm(hyperpost_global)
    if(exists("all_predictions")) rm(all_predictions)
    gc()

  }, error = function(e) {
    cat(paste0("\n!!! ERREUR [n_out=", n_out, ", n_train=", n_train, ", iter=", iter, "] : ", e$message, "\n"))
  })
}

# --- Fin ---
duration_total <- as.numeric(difftime(Sys.time(), t_global_start, units = "hours"))
cat(paste0("\n=== MT TERMINÉ n_out=", n_out, " n_train=", n_train, " ===\n"))
cat(paste0("Durée : ", round(duration_total, 2), "h\n"))

run_info$ended_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_info$duration_hours <- round(duration_total, 2)
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)
