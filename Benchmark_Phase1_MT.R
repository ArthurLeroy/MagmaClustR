# ==========================================================================
# Benchmark_Phase1_MT.R
# Phase 1 : Test paramétrique — MT (MagmaClust single-output, Output_ID=="1")
#
# Usage :
#   Rscript Benchmark_Phase1_MT.R --param=probleme --config=default --seed=1
#
# Charge les données/modèle MOMT générés par Benchmark_Phase1_MOMT.R.
# Entraîne MagmaClust sur l'Output 1 uniquement (single-output, kernel SE).
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
SHARED_HP_CLUSTS <- TRUE
SHARED_HP_TASKS  <- TRUE
N_OUT            <- 2
N_TRAIN          <- 20

if (CONFIG == "variation") {
  switch(PARAM,
    "hp_clusters" = { SHARED_HP_CLUSTS <- FALSE },
    "hp_tasks"    = { SHARED_HP_TASKS <- FALSE },
    "n_out"       = { N_OUT <- 8 },
    "n_train"     = { N_TRAIN <- 100 }
  )
}

N_PRED     <- 10
N_CLUSTERS <- 3

# --- 2. SETUP & LIBRARIES ---
username  <- Sys.getenv("USER")
pkg_dir   <- file.path("/scratch", username, "MagmaClustR")
base_path <- file.path("/scratch", username, "Phase1_experiments", PARAM, CONFIG)

setwd(pkg_dir)

library(tidyverse)
library(readr)
library(MagmaClustR)

cat(paste0("=== Phase 1 MT | param=", PARAM, " config=", CONFIG, " seed=", SEED, " ===\n"))
cat(paste0("  shared_hp_clusts = ", SHARED_HP_CLUSTS, "\n"))
cat(paste0("  shared_hp_tasks  = ", SHARED_HP_TASKS, "\n"))
cat(paste0("  n_out            = ", N_OUT, "\n"))
cat(paste0("  n_train          = ", N_TRAIN, "\n"))
cat(paste0("  n_pred           = ", N_PRED, "\n"))
cat(paste0("  Date             : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste0("  Host             : ", Sys.info()["nodename"], "\n"))
cat(paste0("  PID              : ", Sys.getpid(), "\n\n"))

# Création des répertoires
dir_models_mt      <- file.path(base_path, "Models_MT")
dir_predictions_mt <- file.path(base_path, "Predictions_MT")
dir_run_info       <- file.path(base_path, "run_info")

for (d in c(dir_models_mt, dir_predictions_mt, dir_run_info)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Sauvegarde métadonnées
run_info <- list(
  script = "Benchmark_Phase1_MT.R",
  param = PARAM, config = CONFIG, seed = SEED,
  shared_hp_clusts = SHARED_HP_CLUSTS, shared_hp_tasks = SHARED_HP_TASKS,
  n_out = N_OUT, n_train = N_TRAIN, n_pred = N_PRED,
  started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  hostname = as.character(Sys.info()["nodename"]),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "N/A"),
  r_version = R.version.string, pid = Sys.getpid()
)
run_info_file <- file.path(dir_run_info,
  paste0("mt_seed", SEED, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"))
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)

t_global_start <- Sys.time()

# --- 3. CHARGEMENT DES DONNÉES ---
tryCatch({
  dir_datasets    <- file.path(base_path, "Datasets")
  dir_models_momt <- file.path(base_path, "Models_MOMT")

  file_datasets   <- file.path(dir_datasets, paste0("datasets_seed_", SEED, ".rds"))
  file_model_momt <- file.path(dir_models_momt, paste0("model_seed_", SEED, ".rds"))

  if (!file.exists(file_datasets)) stop("Fichier datasets introuvable : ", file_datasets)
  if (!file.exists(file_model_momt)) stop("Fichier modèle MOMT introuvable : ", file_model_momt)

  datasets           <- readRDS(file_datasets)
  model_momt         <- readRDS(file_model_momt)
  trained_model_momt <- model_momt$trained_model

  train_data          <- datasets$train_data
  pred_tasks_data_all <- datasets$pred_tasks_data
  test_tasks_data_all <- datasets$test_tasks_data
  test_task_ids       <- datasets$test_task_ids

  # --- 4. FILTRAGE OUTPUT_ID == "1" ---
  train_data_o1 <- train_data %>%
    dplyr::filter(Output_ID == "1") %>%
    dplyr::mutate(Output_ID = as.factor(Output_ID)) %>%
    dplyr::select(-any_of("Cluster_ID"))

  # Extraction HP initiaux depuis le modèle MOMT
  ini_hp_k <- trained_model_momt$hp_k %>%
    dplyr::filter(Output_ID == "1") %>%
    dplyr::select(-any_of("l_u_t")) %>%
    dplyr::mutate(se_lengthscale = l_t, se_variance = S_t) %>%
    dplyr::select(-c(l_t, S_t))

  ini_hp_t <- trained_model_momt$hp_t %>%
    dplyr::filter(Output_ID == "1") %>%
    dplyr::filter(Task_ID %in% unique(train_data$Task_ID)) %>%
    dplyr::select(-any_of(c("l_u_t", "Task_ID"))) %>%
    dplyr::mutate(se_lengthscale = l_t, se_variance = S_t) %>%
    dplyr::select(-c(l_t, S_t)) %>%
    dplyr::slice(1)

  # Prior mean pour Output 1 seulement
  indices <- seq(from = 1, to = length(trained_model_momt$ini_args$prior_mean_k), by = N_OUT)
  prior_mean_k_o1 <- trained_model_momt$ini_args$prior_mean_k[indices]
  ini_mixture     <- trained_model_momt$ini_args$ini_mixture

  # --- 5. ENTRAÎNEMENT MAGMACLUST SINGLE-OUTPUT ---
  t_train_start <- Sys.time()
  trained_model_mt <- train_magmaclust(
    data             = train_data_o1,
    ini_mixture      = ini_mixture,
    ini_hp_k         = ini_hp_k,
    nb_cluster       = N_CLUSTERS,
    kern_k           = "SE",
    ini_hp_t         = ini_hp_t,
    kern_t           = "SE",
    shared_hp_tasks  = SHARED_HP_TASKS,
    shared_hp_clusts = SHARED_HP_CLUSTS,
    prior_mean_k     = prior_mean_k_o1,
    pen_diag         = 1e-10
  )
  duration_train <- as.numeric(difftime(Sys.time(), t_train_start, units = "secs"))
  cat(paste0("  Training: ", round(duration_train, 1), "s\n"))

  # Sauvegarde modèle MT
  model_mt_to_save <- list(
    trained_model = trained_model_mt,
    t_training    = duration_train,
    seed          = SEED,
    n_train       = N_TRAIN
  )
  saveRDS(model_mt_to_save, file.path(dir_models_mt, paste0("model_seed_", SEED, ".rds")))

  # --- 6. HYPERPOSTÉRIEUR GLOBAL ---
  all_pred_inputs_global <- list()
  all_test_inputs_global <- list()
  for (tid in test_task_ids) {
    pred_task_o1 <- pred_tasks_data_all[[tid]] %>% dplyr::filter(Output_ID == "1")
    test_task_o1 <- test_tasks_data_all[[tid]] %>% dplyr::filter(Output_ID == "1")
    all_pred_inputs_global[[tid]] <- pred_task_o1 %>% dplyr::select(Output_ID, Input_ID, Input)
    all_test_inputs_global[[tid]] <- test_task_o1 %>% dplyr::select(Input, Input_ID, Output_ID)
  }
  grid_hyperpost_global <- dplyr::bind_rows(
    dplyr::bind_rows(all_pred_inputs_global),
    dplyr::bind_rows(all_test_inputs_global)
  ) %>% dplyr::distinct() %>% dplyr::arrange(Output_ID, Input_ID, Input)

  t_hp_start <- Sys.time()
  hyperpost_global <- hyperposterior_clust(
    trained_model = trained_model_mt,
    data          = trained_model_mt$ini_args$data,
    mixture       = trained_model_mt$hyperpost$mixture,
    hp_k          = trained_model_mt$hp_k,
    hp_t          = trained_model_mt$hp_t,
    grid_inputs   = grid_hyperpost_global,
    kern_k        = trained_model_mt$kern_k,
    kern_t        = trained_model_mt$kern_t,
    prior_mean_k  = trained_model_mt$m_k
  )
  duration_hyperpost <- as.numeric(difftime(Sys.time(), t_hp_start, units = "secs"))
  cat(paste0("  Hyperpost: ", round(duration_hyperpost, 1), "s\n"))

  # --- 7. PRÉDICTION TÂCHE PAR TÂCHE ---
  all_predictions <- list()
  t_pred_total    <- 0

  for (test_task_id in test_task_ids) {
    pred_task_input_data <- pred_tasks_data_all[[test_task_id]] %>%
      dplyr::filter(Output_ID == "1")
    test_data_truth <- test_tasks_data_all[[test_task_id]] %>%
      dplyr::filter(Output_ID == "1")

    test_grid_inputs <- test_data_truth %>%
      dplyr::select(Input, Input_ID, Output_ID) %>%
      dplyr::distinct() %>%
      dplyr::mutate(Output_ID = as.factor(Output_ID))

    pred_data_task <- pred_task_input_data %>%
      dplyr::select(-Task_ID, -any_of("Cluster_ID")) %>%
      dplyr::mutate(Output_ID = as.factor(Output_ID))

    t_pred_start <- Sys.time()
    pred_res <- pred_magmaclust(
      data          = pred_data_task,
      trained_model = trained_model_mt,
      grid_inputs   = test_grid_inputs,
      kern          = "SE",
      hyperpost     = hyperpost_global,
      get_full_cov  = TRUE,
      plot          = FALSE
    )
    t_pred_task <- as.numeric(difftime(Sys.time(), t_pred_start, units = "secs"))
    t_pred_total <- t_pred_total + t_pred_task

    all_predictions[[test_task_id]] <- list(
      prediction  = pred_res,
      truth       = test_data_truth,
      inputs_pred = pred_task_input_data,
      t_pred      = t_pred_task
    )
  }

  # --- 8. SAUVEGARDE ---
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
  saveRDS(predictions_to_save, file.path(dir_predictions_mt, paste0("predictions_seed_", SEED, ".rds")))

  cat(paste0("  Prédictions OK (", round(t_pred_total, 1), "s total)\n"))

  rm(trained_model_mt, hyperpost_global, all_predictions)
  gc()

}, error = function(e) {
  cat(paste0("\n!!! ERREUR [MT param=", PARAM, " config=", CONFIG, " seed=", SEED, "] : ", e$message, "\n"))
  if (!is.null(e$parent)) cat("  Cause: ", conditionMessage(e$parent), "\n")
})

# --- Fin ---
duration_total <- as.numeric(difftime(Sys.time(), t_global_start, units = "mins"))
cat(paste0("\n=== MT TERMINÉ [param=", PARAM, " config=", CONFIG, " seed=", SEED, "] ===\n"))
cat(paste0("Durée : ", round(duration_total, 2), " min\n"))

run_info$ended_at      <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_info$duration_mins <- round(duration_total, 2)
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)
