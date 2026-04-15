# ==========================================================================
# Benchmark_NeurIPS_MO_nclust_forecasting.R
# NeurIPS : Benchmark MO (FORECASTING ONLY) avec paramètre n_clust variable
# Corrige le bug de la version originale : on enlève les observations sur
# l'intervalle de forecasting [0,1] pour TOUS les outputs (pas seulement Output 1).
# MO n'utilise pas de clustering → code identique quel que soit n_clust.
# La seule différence est le chemin des données (qui inclut n_clust).
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
N_CLUST   <- parse_arg_num(args, "n_clust", default = 1)
PROBLEM   <- parse_arg_str(args, "problem")
SEED      <- parse_arg_num(args, "seed")

stopifnot(N_OUT %in% c(2, 4, 8))
stopifnot(N_TRAIN %in% c(15, 30, 100))
stopifnot(N_PRED %in% c(1, 10, 100))
stopifnot(N_CLUST %in% c(1, 2, 3, 4))
stopifnot(PROBLEM == "forecasting")
stopifnot(SEED >= 1 & SEED <= 50)

# --- 1. SETUP & LIBRARIES ---
username  <- Sys.getenv("USER")
pkg_dir   <- file.path("/scratch", username, "MagmaClustR")
base_path <- file.path("/scratch", username, "NeurIPS_experiments",
                        paste0("out", N_OUT, "_train", N_TRAIN, "_pred", N_PRED,
                               "_clust", N_CLUST),
                        PROBLEM)

setwd(pkg_dir)

library(tidyverse)
library(readr)
library(MagmaClustR, lib.loc = "/scratch/agrenoui/R_lib_v2")

convolution_kernel <- MagmaClustR:::convolution_kernel

source(file.path("/scratch", username, "NeurIPS_experiments", "utils", "generate_random_config.R"))

config_label <- paste0("out", N_OUT, "_train", N_TRAIN, "_pred", N_PRED,
                        "_clust", N_CLUST)
n_total_outputs <- N_OUT * (N_TRAIN + N_PRED)

cat(paste0("=== NeurIPS MO | ", config_label, " | ", PROBLEM, " | seed=", SEED, " ===\n"))
cat(paste0("  n_out    = ", N_OUT, "\n"))
cat(paste0("  n_train  = ", N_TRAIN, "\n"))
cat(paste0("  n_pred   = ", N_PRED, "\n"))
cat(paste0("  n_clust  = ", N_CLUST, "\n"))
cat(paste0("  Total MO Output_IDs = ", n_total_outputs, "\n"))
cat(paste0("  problem  = ", PROBLEM, "\n"))
cat(paste0("  Date     : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste0("  Host     : ", Sys.info()["nodename"], "\n"))
cat(paste0("  PID      : ", Sys.getpid(), "\n\n"))

# Création des répertoires
dir_models_mo      <- file.path(base_path, "Models_MO")
dir_predictions_mo <- file.path(base_path, "Predictions_MO")
dir_run_info       <- file.path(base_path, "run_info")

for (d in c(dir_models_mo, dir_predictions_mo, dir_run_info)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

run_info <- list(
  script = "Benchmark_NeurIPS_MO_nclust_forecasting.R",
  n_out = N_OUT, n_train = N_TRAIN, n_pred = N_PRED, n_clust = N_CLUST,
  n_total_outputs = n_total_outputs,
  problem = PROBLEM, seed = SEED,
  started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  hostname = as.character(Sys.info()["nodename"]),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "N/A"),
  r_version = R.version.string, pid = Sys.getpid()
)
run_info_file <- file.path(dir_run_info,
  paste0("mo_seed", SEED, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"))
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)

t_global_start <- Sys.time()

# --- 2. CHARGEMENT DES DONNÉES ---
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

  data_obs            <- datasets$data_obs
  scale_factors       <- datasets$scale_factors
  train_task_ids      <- datasets$train_task_ids
  test_task_ids       <- datasets$test_task_ids
  pred_tasks_data_all <- datasets$pred_tasks_data
  test_tasks_data_all <- datasets$test_tasks_data

  # --- 3. REFORMATAGE POUR MO ---
  all_task_ids_ordered <- c(train_task_ids, test_task_ids)
  output_ids_orig <- as.character(1:N_OUT)

  mo_id_mapping <- tidyr::expand_grid(
    Task_ID = all_task_ids_ordered,
    orig_Output_ID = output_ids_orig
  ) %>%
    dplyr::mutate(mo_Output_ID = as.character(dplyr::row_number()))

  cat("  Mapping MO Output_IDs (premières lignes) :\n")
  print(head(mo_id_mapping, 10))

  # 3a. Données d'entraînement
  train_data_mo_parts <- list()
  for (tid in train_task_ids) {
    task_data <- data_obs %>%
      dplyr::filter(Task_ID == tid) %>%
      dplyr::select(-any_of(c("Task_ID", "Cluster_ID")))

    task_mapping <- mo_id_mapping %>%
      dplyr::filter(Task_ID == tid) %>%
      dplyr::select(orig_Output_ID, mo_Output_ID)
    task_data <- task_data %>%
      dplyr::left_join(task_mapping, by = c("Output_ID" = "orig_Output_ID")) %>%
      dplyr::mutate(Output_ID = mo_Output_ID) %>%
      dplyr::select(-mo_Output_ID)

    train_data_mo_parts[[tid]] <- task_data
  }

  # 3b. Données de contexte des tâches pred
  pred_context_mo_parts <- list()
  for (tid in test_task_ids) {
    task_mapping <- mo_id_mapping %>%
      dplyr::filter(Task_ID == tid) %>%
      dplyr::select(orig_Output_ID, mo_Output_ID)
    pred_data <- pred_tasks_data_all[[tid]] %>%
      dplyr::select(-any_of(c("Task_ID", "Cluster_ID"))) %>%
      dplyr::left_join(task_mapping, by = c("Output_ID" = "orig_Output_ID")) %>%
      dplyr::mutate(Output_ID = mo_Output_ID) %>%
      dplyr::select(-mo_Output_ID)

    pred_context_mo_parts[[tid]] <- pred_data
  }

  # 3c. Assemblage complet
  all_mo_data <- dplyr::bind_rows(
    dplyr::bind_rows(train_data_mo_parts),
    dplyr::bind_rows(pred_context_mo_parts)
  )

  cat(paste0("  MO data: ", nrow(all_mo_data), " observations, ",
             length(unique(all_mo_data$Output_ID)), " Output_IDs\n"))

  # 3d. Grid inputs pour la prédiction (TOUS les outputs en [0,1])
  test_grid_parts <- list()
  for (tid in test_task_ids) {
    for (oid in output_ids_orig) {
      mo_oid <- mo_id_mapping %>%
        dplyr::filter(Task_ID == tid, orig_Output_ID == oid) %>%
        dplyr::pull(mo_Output_ID)
      test_data_oid <- test_tasks_data_all[[tid]] %>%
        dplyr::filter(Output_ID == oid) %>%
        dplyr::mutate(Output_ID = mo_oid) %>%
        dplyr::select(Input, Input_ID, Output_ID) %>%
        dplyr::distinct()
      if (nrow(test_data_oid) > 0) {
        test_grid_parts[[paste0(tid, "_", oid)]] <- test_data_oid
      }
    }
  }
  test_grid_inputs <- dplyr::bind_rows(test_grid_parts)

  # --- 4. INITIALISATION HP ---
  ini_hp_momt <- trained_model_momt$hp_t %>%
    dplyr::slice(1:N_OUT) %>%
    dplyr::select(-Task_ID)

  all_output_ids_mo <- as.character(sort(as.numeric(unique(c(
    all_mo_data$Output_ID, test_grid_inputs$Output_ID)))))

  ini_hp_mo <- hp(
    kern = convolution_kernel,
    list_task_ID = "1",
    list_output_ID = all_output_ids_mo,
    shared_hp_outputs = FALSE,
    shared_hp_tasks = TRUE,
    noise = TRUE,
    hp_config = tibble(
      output_id = all_output_ids_mo,
      lt_min = rep(ini_hp_momt$l_t[1], length(all_output_ids_mo)),
      lt_max = rep(ini_hp_momt$l_t[1], length(all_output_ids_mo)),
      St_min = rep(0, length(all_output_ids_mo)),
      St_max = rep(0, length(all_output_ids_mo)),
      lu_min = rep(ini_hp_momt$l_u_t[1], length(all_output_ids_mo)),
      lu_max = rep(ini_hp_momt$l_u_t[1], length(all_output_ids_mo)),
      noise_min = -3, noise_max = -3
    )
  ) %>% dplyr::select(-Task_ID)

  prior_means_mo <- NULL

  # --- 5. ENTRAÎNEMENT ---
  t_train_start <- Sys.time()
  hp_optim <- train_gp(
    data       = all_mo_data,
    ini_hp     = ini_hp_mo,
    kern       = convolution_kernel,
    prior_mean = prior_means_mo
  )
  duration_train <- as.numeric(difftime(Sys.time(), t_train_start, units = "secs"))
  cat(paste0("  Training: ", round(duration_train, 1), "s\n"))

  # --- 6. PRÉDICTION ---
  t_pred_start <- Sys.time()
  pred_res <- pred_gp(
    data         = all_mo_data,
    kern         = convolution_kernel,
    hp           = hp_optim,
    mean         = prior_means_mo,
    grid_inputs  = test_grid_inputs,
    get_full_cov = TRUE,
    plot         = FALSE
  )
  duration_pred <- as.numeric(difftime(Sys.time(), t_pred_start, units = "secs"))
  cat(paste0("  Prediction: ", round(duration_pred, 1), "s\n"))

  # --- 7. RÉORGANISER LES PRÉDICTIONS PAR TÂCHE (tous les outputs) ---
  all_predictions <- list()
  pred_df <- pred_res$pred_gp

  for (test_task_id in test_task_ids) {
    task_pred_by_output <- list()
    for (oid in output_ids_orig) {
      mo_output_id <- mo_id_mapping %>%
        dplyr::filter(Task_ID == test_task_id, orig_Output_ID == oid) %>%
        dplyr::pull(mo_Output_ID)

      task_pred <- pred_df$pred %>%
        dplyr::filter(Output_ID == mo_output_id)

      if (nrow(task_pred) > 0) {
        task_pred_by_output[[oid]] <- task_pred
      }
    }

    task_cluster <- datasets$data_raw %>%
      dplyr::filter(Task_ID == test_task_id) %>%
      dplyr::pull(Cluster_ID) %>%
      unique()

    all_predictions[[test_task_id]] <- list(
      predictions_by_output = task_pred_by_output,
      truth       = test_tasks_data_all[[test_task_id]],
      cluster_id  = task_cluster[1]
    )
  }

  # --- 8. SAUVEGARDE ---
  models_to_save <- list(
    hp_optim       = hp_optim,
    prior_means    = prior_means_mo,
    mo_id_mapping  = mo_id_mapping,
    t_training     = duration_train,
    seed           = SEED,
    n_train        = N_TRAIN,
    n_pred         = N_PRED,
    n_clust        = N_CLUST
  )

  predictions_to_save <- list(
    predictions    = all_predictions,
    pred_gp_full   = pred_res,
    mo_id_mapping  = mo_id_mapping,
    t_train_total  = duration_train,
    t_pred_total   = duration_pred,
    seed           = SEED,
    n_clust        = N_CLUST,
    test_task_ids  = test_task_ids,
    n_train        = N_TRAIN,
    n_pred         = N_PRED
  )

  saveRDS(models_to_save,
          file.path(dir_models_mo, paste0("models_seed_", SEED, ".rds")))
  saveRDS(predictions_to_save,
          file.path(dir_predictions_mo, paste0("predictions_seed_", SEED, ".rds")))

  cat(paste0("\n  Total: train=", round(duration_train, 1),
             "s, pred=", round(duration_pred, 1), "s\n"))

  rm(all_predictions, pred_res)
  gc()

}, error = function(e) {
  cat(paste0("\n!!! ERREUR [MO out=", N_OUT, " train=", N_TRAIN,
             " pred=", N_PRED, " clust=", N_CLUST, " ", PROBLEM,
             " seed=", SEED, "] : ", e$message, "\n"))
  if (!is.null(e$parent)) cat("  Cause: ", conditionMessage(e$parent), "\n")
})

# --- Fin ---
duration_total <- as.numeric(difftime(Sys.time(), t_global_start, units = "mins"))
cat(paste0("\n=== MO TERMINÉ [out=", N_OUT, " train=", N_TRAIN,
           " pred=", N_PRED, " clust=", N_CLUST, " ", PROBLEM,
           " seed=", SEED, "] ===\n"))
cat(paste0("Durée : ", round(duration_total, 2), " min\n"))

run_info$ended_at      <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_info$duration_mins <- round(duration_total, 2)
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)
