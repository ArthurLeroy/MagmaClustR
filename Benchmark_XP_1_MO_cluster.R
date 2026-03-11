# ==========================================================================
# Benchmark_XP_1_MO_cluster.R
# Version adaptée pour le cluster IMT (SLURM) — PARAMÉTRÉ
#
# Usage :
#   Rscript Benchmark_XP_1_MO_cluster.R --n_out=3 --n_train=15
#
# Multi-Output GP sans clustering (train_gp + pred_gp avec convolution_kernel).
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
library(MagmaClustR)
# library(devtools)
# devtools::load_all()

convolution_kernel <- MagmaClustR:::convolution_kernel

n_out <- N_OUT_TARGET
n_train <- N_TRAIN_TARGET
n_iterations <- 100
n_pred_max <- 1000
n_pred_subsets <- c(1, 10, 100, 1000)

cat(paste0("=== Benchmark XP1 MO (Multi-Output GP) ===\n"))
cat(paste0("  n_out   = ", n_out, "\n"))
cat(paste0("  n_train = ", n_train, "\n"))
cat(paste0("  Date    : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste0("  Host    : ", Sys.info()["nodename"], "\n"))
cat(paste0("  PID     : ", Sys.getpid(), "\n\n"))

# Sauvegarde métadonnées
run_info <- list(
  script = "Benchmark_XP_1_MO_cluster.R",
  n_out = n_out, n_train = n_train, n_iterations = n_iterations,
  started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  hostname = as.character(Sys.info()["nodename"]),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "N/A"),
  r_version = R.version.string, pid = Sys.getpid()
)
run_info_dir <- file.path(base_path, "run_info")
if (!dir.exists(run_info_dir)) dir.create(run_info_dir, recursive = TRUE)
run_info_file <- file.path(run_info_dir,
  paste0("mo_nout", n_out, "_ntrain", n_train, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"))
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)

t_global_start <- Sys.time()

# --- 2. BOUCLE SUR LES ITÉRATIONS ---
for (iter in 1:n_iterations) {
  tryCatch({
    cat(paste0("[n_out=", n_out, " n_train=", n_train, "] MO Iter ", iter, "/", n_iterations,
               " [", format(Sys.time(), "%H:%M:%S"), "]... "))

    # dir_n_out <- file.path(base_path, paste0("n_out_", n_out))
    # dir_datasets <- file.path(dir_n_out, "Datasets", paste0("train_", n_train, "_pred_", n_pred_max))
    # dir_models_momt <- file.path(dir_n_out, "Models_MOMT", paste0("train_", n_train, "_pred_", n_pred_max))

    dir_datasets <- file.path(base_path, "Datasets", paste0("n_out_", n_out), paste0("train_", n_train, "_pred_", n_pred_max))
    dir_models_momt <- file.path(base_path, "Models_MOMT", paste0("n_out_", n_out), paste0("train_", n_train, "_pred_", n_pred_max))

    file_name_datasets <- file.path(dir_datasets, paste0("datasets_iter_", iter, ".rds"))
    file_name_model_momt <- file.path(dir_models_momt, paste0("model_iter_", iter, ".rds"))

    if (!file.exists(file_name_datasets)) { cat("SKIP (pas de données)\n"); next }
    if (!file.exists(file_name_model_momt)) { cat("SKIP (pas de modèle MOMT)\n"); next }

    datasets <- readRDS(file_name_datasets)
    model_momt <- readRDS(file_name_model_momt)
    trained_model_momt <- model_momt$trained_model

    pred_tasks_data_all <- datasets$pred_tasks_data
    test_tasks_data_all <- datasets$test_tasks_data
    test_task_ids_max <- datasets$test_task_ids

    ini_hp_template <- trained_model_momt$hp_t %>%
      dplyr::slice(1:n_out) %>% dplyr::select(-Task_ID)

    # Boucle sur les subsets
    for (n_pred in n_pred_subsets) {
      n_pred_actual <- min(n_pred, length(test_task_ids_max))
      if (n_pred_actual == 0) next

      test_task_ids <- test_task_ids_max[1:n_pred_actual]

      # dir_predictions_mo <- file.path(dir_n_out, "Predictions_MO", paste0("train_", n_train, "_pred_", n_pred))
      # dir_models_mo <- file.path(dir_n_out, "Models_MO", paste0("train_", n_train, "_pred_", n_pred))
      dir_predictions_mo <- file.path(base_path, "Predictions_MO", paste0("n_out_", n_out), paste0("train_", n_train, "_pred_", n_pred))
      dir_models_mo <- file.path(base_path, "Models_MO", paste0("n_out_", n_out), paste0("train_", n_train, "_pred_", n_pred))
      if (!dir.exists(dir_predictions_mo)) dir.create(dir_predictions_mo, recursive = TRUE)
      if (!dir.exists(dir_models_mo)) dir.create(dir_models_mo, recursive = TRUE)

      all_predictions <- list()
      all_models <- list()
      t_train_total <- 0
      t_pred_total <- 0

      for (test_task_id in test_task_ids) {
        pred_task_input_data <- pred_tasks_data_all[[test_task_id]]
        test_data_truth <- test_tasks_data_all[[test_task_id]]

        train_data_task <- pred_task_input_data %>%
          dplyr::select(-Task_ID, -Cluster_ID) %>%
          dplyr::mutate(Output_ID = as.factor(Output_ID))

        mean_emp_per_output <- train_data_task %>%
          group_by(Output_ID) %>%
          summarise(Empirical_mean = mean(Output, na.rm = TRUE), .groups = 'drop')

        test_grid_inputs <- test_data_truth %>%
          dplyr::select(Input, Input_ID, Output_ID) %>%
          dplyr::distinct() %>% dplyr::mutate(Output_ID = as.factor(Output_ID))

        # Entraînement
        t_train_start <- Sys.time()
        hp_optim <- train_gp(
          data = train_data_task, ini_hp = ini_hp_template,
          kern = convolution_kernel,
          prior_mean = mean_emp_per_output$Empirical_mean %>% as.vector()
        )
        t_train_task <- as.numeric(difftime(Sys.time(), t_train_start, units = "secs"))
        t_train_total <- t_train_total + t_train_task

        # Prédiction
        t_pred_start <- Sys.time()
        pred_res <- pred_gp(
          data = train_data_task, kern = convolution_kernel,
          hp = hp_optim,
          mean = mean_emp_per_output$Empirical_mean %>% as.vector(),
          grid_inputs = test_grid_inputs, get_full_cov = TRUE, plot = FALSE
        )
        t_pred_task <- as.numeric(difftime(Sys.time(), t_pred_start, units = "secs"))
        t_pred_total <- t_pred_total + t_pred_task

        all_predictions[[test_task_id]] <- list(
          prediction = pred_res, truth = test_data_truth,
          inputs_pred = pred_task_input_data,
          t_train = t_train_task, t_pred = t_pred_task
        )
        all_models[[test_task_id]] <- list(
          hp_optim = hp_optim,
          prior_mean = mean_emp_per_output$Empirical_mean
        )
      }

      # Sauvegarde
      predictions_to_save <- list(
        predictions = all_predictions, t_train_total = t_train_total,
        t_pred_total = t_pred_total, seed = datasets$seed,
        test_task_ids = test_task_ids, n_train = n_train, n_pred = n_pred
      )
      models_to_save <- list(
        models = all_models, t_train_total = t_train_total,
        seed = datasets$seed, n_train = n_train, n_pred = n_pred
      )

      saveRDS(predictions_to_save, file.path(dir_predictions_mo, paste0("predictions_iter_", iter, ".rds")))
      saveRDS(models_to_save, file.path(dir_models_mo, paste0("models_iter_", iter, ".rds")))
    }

    cat(paste0("OK\n"))
    if(exists("all_predictions")) rm(all_predictions)
    if(exists("all_models")) rm(all_models)
    gc()

  }, error = function(e) {
    cat(paste0("\n!!! ERREUR [n_out=", n_out, ", n_train=", n_train, ", iter=", iter, "] : ", e$message, "\n"))
  })
}

# --- Fin ---
duration_total <- as.numeric(difftime(Sys.time(), t_global_start, units = "hours"))
cat(paste0("\n=== MO TERMINÉ n_out=", n_out, " n_train=", n_train, " ===\n"))
cat(paste0("Durée : ", round(duration_total, 2), "h\n"))

run_info$ended_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_info$duration_hours <- round(duration_total, 2)
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)
