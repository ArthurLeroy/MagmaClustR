# ==========================================================================
# Rerun_Predictions_NeurIPS_MT_nclust_interpolation.R
# NeurIPS : relance uniquement les prédictions MT (interpolation)
# à partir des datasets et modèles déjà sauvegardés sur le serveur.
#
# Usage :
#   Rscript Rerun_Predictions_NeurIPS_MT_nclust_interpolation.R \
#     --n_out=2 --n_train=30 --n_pred=1 --n_clust=1 \
#     --problem=interpolation --seed=1
# ==========================================================================

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
  as.numeric(val)
}

N_OUT   <- parse_arg_num(args, "n_out")
N_TRAIN <- parse_arg_num(args, "n_train")
N_PRED  <- parse_arg_num(args, "n_pred")
N_CLUST <- parse_arg_num(args, "n_clust", default = 1)
PROBLEM <- parse_arg_str(args, "problem")
SEED    <- parse_arg_num(args, "seed")

stopifnot(N_OUT %in% c(2, 4, 8))
stopifnot(N_TRAIN %in% c(15, 30, 100))
stopifnot(N_PRED %in% c(1, 10, 100))
stopifnot(N_CLUST %in% c(1, 2, 3, 4))
stopifnot(PROBLEM == "interpolation")
stopifnot(SEED >= 1 & SEED <= 50)

username  <- Sys.getenv("USER")
pkg_dir   <- file.path("/scratch", username, "MagmaClustR")
base_root <- file.path("/scratch", username, "NeurIPS_experiments_interpolation")

setwd(pkg_dir)

library(tidyverse)
library(readr)
library(MagmaClustR, lib.loc = "/scratch/agrenoui/R_lib_v2")

config_label <- paste0("out", N_OUT, "_train", N_TRAIN, "_pred", N_PRED,
                       "_clust", N_CLUST)

dir_datasets       <- file.path(base_root, "Datasets", config_label)
dir_models_mt      <- file.path(base_root, "Models_MT", config_label)
dir_predictions_mt <- file.path(base_root, "Predictions_MT", config_label)
dir_run_info       <- file.path(base_root, "run_info", config_label)

for (d in c(dir_predictions_mt, dir_run_info)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

run_info <- list(
  script = "Rerun_Predictions_NeurIPS_MT_nclust_interpolation.R",
  n_out = N_OUT,
  n_train = N_TRAIN,
  n_pred = N_PRED,
  n_clust = N_CLUST,
  problem = PROBLEM,
  seed = SEED,
  started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  hostname = as.character(Sys.info()["nodename"]),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "N/A"),
  r_version = R.version.string,
  pid = Sys.getpid()
)
run_info_file <- file.path(
  dir_run_info,
  paste0("rerun_pred_mt_seed", SEED, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json")
)
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)

t_global_start <- Sys.time()

tryCatch({
  file_datasets <- file.path(dir_datasets, paste0("datasets_seed_", SEED, ".rds"))
  file_model_mt <- file.path(dir_models_mt, paste0("model_seed_", SEED, ".rds"))

  if (!file.exists(file_datasets)) stop("Fichier datasets introuvable : ", file_datasets)
  if (!file.exists(file_model_mt)) stop("Fichier modèle MT introuvable : ", file_model_mt)

  datasets <- readRDS(file_datasets)
  model_mt <- readRDS(file_model_mt)

  trained_model_mt <- model_mt$trained_model
  data_raw         <- datasets$data_raw
  train_task_ids   <- datasets$train_task_ids
  test_task_ids    <- datasets$test_task_ids

  if (is.null(trained_model_mt)) stop("Objet trained_model absent du fichier modèle MT.")

  pred_tasks_data_raw <- list()
  test_tasks_data_raw <- list()
  gap_start <- -0.5
  gap_end   <-  0.5

  train_data_raw <- data_raw %>%
    dplyr::filter(Task_ID %in% train_task_ids) %>%
    dplyr::filter(Output_ID != "1" |
                  (Output_ID == "1" & !(Input >= gap_start & Input <= gap_end)))

  for (tid in test_task_ids) {
    task_full <- data_raw %>% dplyr::filter(Task_ID == tid)
    task_O1 <- task_full %>% dplyr::filter(Output_ID == "1") %>% dplyr::arrange(Input)

    test_tasks_data_raw[[tid]] <- task_O1 %>%
      dplyr::filter(Input >= gap_start & Input <= gap_end)

    pred_tasks_data_raw[[tid]] <- dplyr::bind_rows(
      task_full %>% dplyr::filter(Output_ID != "1"),
      task_O1 %>% dplyr::filter(!(Input >= gap_start & Input <= gap_end))
    ) %>%
      dplyr::arrange(Output_ID, Input)
  }

  output_ids <- unique(data_raw$Output_ID)

  all_pred_inputs_global <- list()
  all_test_inputs_global <- list()
  for (tid in test_task_ids) {
    for (oid in output_ids) {
      mt_task_id <- paste0("T", tid, "_O", oid)
      pred_data_oid <- pred_tasks_data_raw[[tid]] %>%
        dplyr::filter(Output_ID == oid)
      if (nrow(pred_data_oid) > 0) {
        all_pred_inputs_global[[mt_task_id]] <- pred_data_oid %>%
          dplyr::mutate(Output_ID = "1") %>%
          dplyr::select(Output_ID, Input_ID, Input)
      }
      test_data_oid <- test_tasks_data_raw[[tid]] %>%
        dplyr::filter(Output_ID == oid)
      if (nrow(test_data_oid) > 0) {
        all_test_inputs_global[[mt_task_id]] <- test_data_oid %>%
          dplyr::mutate(Output_ID = "1") %>%
          dplyr::select(Input, Input_ID, Output_ID)
      }
    }
  }

  grid_hyperpost_global <- dplyr::bind_rows(
    dplyr::bind_rows(all_pred_inputs_global),
    dplyr::bind_rows(all_test_inputs_global)
  ) %>% dplyr::distinct() %>% dplyr::arrange(Output_ID, Input_ID, Input)

  t_hp_start <- Sys.time()
  if (N_CLUST == 1) {
    hyperpost_global <- hyperposterior(
      trained_model = trained_model_mt,
      data          = trained_model_mt$ini_args$data,
      hp_0          = trained_model_mt$hp_0,
      hp_t          = trained_model_mt$hp_t,
      grid_inputs   = grid_hyperpost_global,
      kern_0        = trained_model_mt$kern_0,
      kern_t        = trained_model_mt$kern_t,
      prior_mean    = trained_model_mt$m_0
    )
  } else {
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
  }
  duration_hyperpost <- as.numeric(difftime(Sys.time(), t_hp_start, units = "secs"))
  cat(paste0("  Hyperpost: ", round(duration_hyperpost, 1), "s\n"))

  all_predictions <- list()
  t_pred_total <- 0

  for (test_task_id in test_task_ids) {
    task_predictions <- list()

    for (oid in output_ids) {
      pred_data_oid <- pred_tasks_data_raw[[test_task_id]] %>%
        dplyr::filter(Output_ID == oid) %>%
        dplyr::select(-any_of(c("Task_ID", "Cluster_ID", "Dataframe"))) %>%
        dplyr::mutate(Output_ID = as.factor("1"))

      if (nrow(pred_data_oid) == 0) next

      test_data_oid <- test_tasks_data_raw[[test_task_id]] %>%
        dplyr::filter(Output_ID == oid)
      if (nrow(test_data_oid) == 0) next

      test_grid_inputs <- test_data_oid %>%
        dplyr::mutate(Output_ID = as.factor("1")) %>%
        dplyr::select(Input, Input_ID, Output_ID) %>%
        dplyr::distinct()

      t_pred_start <- Sys.time()
      if (N_CLUST == 1) {
        pred_res <- pred_magma(
          data          = pred_data_oid,
          trained_model = trained_model_mt,
          grid_inputs   = test_grid_inputs,
          kern          = "SE",
          hyperpost     = hyperpost_global,
          get_full_cov  = TRUE,
          plot          = FALSE
        )
      } else {
        pred_res <- pred_magmaclust(
          data          = pred_data_oid,
          trained_model = trained_model_mt,
          grid_inputs   = test_grid_inputs,
          kern          = "SE",
          hyperpost     = hyperpost_global,
          get_full_cov  = TRUE,
          plot          = FALSE
        )
      }
      t_pred_task  <- as.numeric(difftime(Sys.time(), t_pred_start, units = "secs"))
      t_pred_total <- t_pred_total + t_pred_task

      task_predictions[[oid]] <- list(
        prediction = pred_res,
        t_pred     = t_pred_task
      )
    }

    all_predictions[[test_task_id]] <- list(
      predictions_by_output = task_predictions,
      truth                 = test_tasks_data_raw[[test_task_id]],
      inputs_pred           = pred_tasks_data_raw[[test_task_id]],
      t_pred                = sum(sapply(task_predictions, function(x) x$t_pred))
    )
  }

  predictions_to_save <- list(
    predictions        = all_predictions,
    t_hyperpost_global = duration_hyperpost,
    t_pred_total       = t_pred_total,
    t_training         = model_mt$t_training,
    seed               = SEED,
    n_clust            = N_CLUST,
    test_task_ids      = test_task_ids,
    n_train            = N_TRAIN,
    n_pred             = N_PRED
  )

  saveRDS(
    predictions_to_save,
    file.path(dir_predictions_mt, paste0("predictions_seed_", SEED, ".rds"))
  )

  cat(paste0("  Prédictions OK (", round(t_pred_total, 1), "s total)\n"))
}, error = function(e) {
  cat(paste0("\n!!! ERREUR [RERUN MT out=", N_OUT,
             " train=", N_TRAIN,
             " pred=", N_PRED,
             " clust=", N_CLUST,
             " ", PROBLEM,
             " seed=", SEED,
             "] : ", e$message, "\n"))
  if (!is.null(e$parent)) cat("  Cause: ", conditionMessage(e$parent), "\n")
})

duration_total <- as.numeric(difftime(Sys.time(), t_global_start, units = "mins"))
cat(paste0("\n=== RERUN MT TERMINÉ [out=", N_OUT,
           " train=", N_TRAIN,
           " pred=", N_PRED,
           " clust=", N_CLUST,
           " ", PROBLEM,
           " seed=", SEED,
           "] ===\n"))
cat(paste0("Durée : ", round(duration_total, 2), " min\n"))

run_info$ended_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_info$duration_mins <- round(duration_total, 2)
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)