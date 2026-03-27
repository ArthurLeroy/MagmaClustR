# ==========================================================================
# Benchmark_Phase1_NC_MO.R
# Phase 1 NC : Test paramétrique — MO (Multi-Output GP sans clustering)
#
# Usage :
#   Rscript Benchmark_Phase1_NC_MO.R --param=probleme --config=default --seed=1
#
# Charge les données générées par le script NC_MOMT.
# Entraîne un GP multi-output (train_gp + pred_gp avec convolution_kernel)
# pour chaque tâche de test.
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

stopifnot(PARAM %in% c("probleme", "hp_clusters", "hp_tasks", "inputs", "n_train"))
stopifnot(CONFIG %in% c("default", "variation"))
stopifnot(SEED >= 1 & SEED <= 10)

# --- 1. CONFIGURATION PAR DÉFAUT ---
N_OUT   <- 2
N_TRAIN <- 20

if (CONFIG == "variation") {
  switch(PARAM,
    "n_train" = { N_TRAIN <- 100 }
  )
}

N_PRED <- 10

# --- 2. SETUP & LIBRARIES ---
username  <- Sys.getenv("USER")
pkg_dir   <- file.path("/scratch", username, "MagmaClustR")
base_path <- file.path("/scratch", username, "Phase1_NC_experiments", PARAM, CONFIG)

setwd(pkg_dir)

library(tidyverse)
library(readr)
library(MagmaClustR)

convolution_kernel <- MagmaClustR:::convolution_kernel

cat(paste0("=== Phase 1 NC MO | param=", PARAM, " config=", CONFIG, " seed=", SEED, " ===\n"))
cat(paste0("  n_out   = ", N_OUT, "\n"))
cat(paste0("  n_train = ", N_TRAIN, "\n"))
cat(paste0("  n_pred  = ", N_PRED, "\n"))
cat(paste0("  Date    : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste0("  Host    : ", Sys.info()["nodename"], "\n"))
cat(paste0("  PID     : ", Sys.getpid(), "\n\n"))

# Création des répertoires
dir_models_mo      <- file.path(base_path, "Models_MO")
dir_predictions_mo <- file.path(base_path, "Predictions_MO")
dir_run_info       <- file.path(base_path, "run_info")

for (d in c(dir_models_mo, dir_predictions_mo, dir_run_info)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Sauvegarde métadonnées
run_info <- list(
  script = "Benchmark_Phase1_NC_MO.R",
  param = PARAM, config = CONFIG, seed = SEED,
  n_out = N_OUT, n_train = N_TRAIN, n_pred = N_PRED,
  started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  hostname = as.character(Sys.info()["nodename"]),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "N/A"),
  r_version = R.version.string, pid = Sys.getpid()
)
run_info_file <- file.path(dir_run_info,
  paste0("mo_seed", SEED, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"))
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)

t_global_start <- Sys.time()

# --- 3. CHARGEMENT DES DONNÉES ---
tryCatch({
  dir_datasets    <- file.path(base_path, "Datasets")
  dir_models_momt <- file.path(base_path, "Models_MOMT")

  file_datasets <- file.path(dir_datasets, paste0("datasets_seed_", SEED, ".rds"))
  file_model_momt <- file.path(dir_models_momt, paste0("model_seed_", SEED, ".rds"))

  if (!file.exists(file_datasets)) stop("Fichier datasets introuvable : ", file_datasets)
  if (!file.exists(file_model_momt)) stop("Fichier modèle MOMT NC introuvable : ", file_model_momt)

  datasets         <- readRDS(file_datasets)
  model_momt       <- readRDS(file_model_momt)
  trained_model_momt <- model_momt$trained_model

  pred_tasks_data_all <- datasets$pred_tasks_data
  test_tasks_data_all <- datasets$test_tasks_data
  test_task_ids       <- datasets$test_task_ids

  # HP template depuis le modèle MOMT NC
  ini_hp_template <- trained_model_momt$hp_t %>%
    dplyr::slice(1:N_OUT) %>%
    dplyr::select(-any_of("Task_ID"))

  # --- 4. ENTRAÎNEMENT + PRÉDICTION PAR TÂCHE ---
  all_predictions <- list()
  all_models      <- list()

  for (test_task_id in test_task_ids) {
    cat(paste0("  MO tâche ", test_task_id, "... "))

    pred_task_input_data <- pred_tasks_data_all[[test_task_id]]
    test_data_truth      <- test_tasks_data_all[[test_task_id]]

    train_data_task <- pred_task_input_data %>%
      dplyr::select(-Task_ID) %>%
      dplyr::mutate(Output_ID = as.factor(Output_ID))

    mean_emp_per_output <- train_data_task %>%
      dplyr::group_by(Output_ID) %>%
      dplyr::summarise(Empirical_mean = mean(Output, na.rm = TRUE), .groups = "drop")

    test_grid_inputs <- test_data_truth %>%
      dplyr::select(Input, Input_ID, Output_ID) %>%
      dplyr::distinct() %>%
      dplyr::mutate(Output_ID = as.factor(Output_ID))

    # Entraînement
    t_train_start <- Sys.time()
    hp_optim <- train_gp(
      data       = train_data_task,
      ini_hp     = ini_hp_template,
      kern       = convolution_kernel,
      prior_mean = mean_emp_per_output$Empirical_mean %>% as.vector()
    )
    t_train_task <- as.numeric(difftime(Sys.time(), t_train_start, units = "secs"))

    # Prédiction
    t_pred_start <- Sys.time()
    pred_res <- pred_gp(
      data        = train_data_task,
      kern        = convolution_kernel,
      hp          = hp_optim,
      mean        = mean_emp_per_output$Empirical_mean %>% as.vector(),
      grid_inputs = test_grid_inputs,
      get_full_cov = TRUE,
      plot        = FALSE
    )
    t_pred_task <- as.numeric(difftime(Sys.time(), t_pred_start, units = "secs"))

    all_predictions[[test_task_id]] <- list(
      prediction  = pred_res,
      truth       = test_data_truth,
      inputs_pred = pred_task_input_data,
      t_train     = t_train_task,
      t_pred      = t_pred_task
    )
    all_models[[test_task_id]] <- list(
      hp_optim   = hp_optim,
      prior_mean = mean_emp_per_output$Empirical_mean
    )

    cat(paste0("OK (train=", round(t_train_task, 1), "s, pred=", round(t_pred_task, 1), "s)\n"))
  }

  # --- 5. SAUVEGARDE ---
  t_train_total <- sum(sapply(all_predictions, function(x) x$t_train))
  t_pred_total  <- sum(sapply(all_predictions, function(x) x$t_pred))

  predictions_to_save <- list(
    predictions   = all_predictions,
    t_train_total = t_train_total,
    t_pred_total  = t_pred_total,
    seed          = SEED,
    test_task_ids = test_task_ids,
    n_train       = N_TRAIN,
    n_pred        = N_PRED
  )

  models_to_save <- list(
    models        = all_models,
    t_train_total = t_train_total,
    seed          = SEED,
    n_train       = N_TRAIN,
    n_pred        = N_PRED
  )

  saveRDS(predictions_to_save, file.path(dir_predictions_mo, paste0("predictions_seed_", SEED, ".rds")))
  saveRDS(models_to_save, file.path(dir_models_mo, paste0("models_seed_", SEED, ".rds")))

  cat(paste0("\n  Total: train=", round(t_train_total, 1), "s, pred=", round(t_pred_total, 1), "s\n"))

  rm(all_predictions, all_models)
  gc()

}, error = function(e) {
  cat(paste0("\n!!! ERREUR [NC MO param=", PARAM, " config=", CONFIG, " seed=", SEED, "] : ", e$message, "\n"))
})

# --- Fin ---
duration_total <- as.numeric(difftime(Sys.time(), t_global_start, units = "mins"))
cat(paste0("\n=== NC MO TERMINÉ [param=", PARAM, " config=", CONFIG, " seed=", SEED, "] ===\n"))
cat(paste0("Durée : ", round(duration_total, 2), " min\n"))

run_info$ended_at      <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_info$duration_mins <- round(duration_total, 2)
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)
