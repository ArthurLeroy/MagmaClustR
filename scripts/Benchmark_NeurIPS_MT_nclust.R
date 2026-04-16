# ==========================================================================
# Benchmark_NeurIPS_MT_nclust.R
# NeurIPS : Benchmark MT avec paramètre n_clust variable (1, 2, 3, 4)
#
# Usage :
#   Rscript Benchmark_NeurIPS_MT_nclust.R --n_out=2 --n_train=30 --n_pred=1 \
#     --n_clust=1 --problem=interpolation --seed=1
#
# Charge les données générées par MOMT_nclust.
# Réorganise : chaque (Task_ID, Output_ID) de MOMT → une "tâche" MT avec
#   Output_ID="1", kernel SE.
#
# Branchement :
#   - n_clust=1  : train_magma(), hyperposterior(), pred_magma()
#   - n_clust>=2 : train_magmaclust(), hyperposterior_clust(), pred_magmaclust()
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
stopifnot(PROBLEM %in% c("interpolation", "forecasting"))
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

source(file.path("/scratch", username, "NeurIPS_experiments", "utils", "generate_random_config.R"))

config_label <- paste0("out", N_OUT, "_train", N_TRAIN, "_pred", N_PRED,
                        "_clust", N_CLUST)

cat(paste0("=== NeurIPS MT | ", config_label, " | ", PROBLEM, " | seed=", SEED, " ===\n"))
cat(paste0("  n_out    = ", N_OUT, " (→ MT tasks = n_tasks × n_out)\n"))
cat(paste0("  n_train  = ", N_TRAIN, " (→ MT train tasks = ", N_TRAIN * N_OUT, ")\n"))
cat(paste0("  n_pred   = ", N_PRED, " (→ MT pred tasks = ", N_PRED * N_OUT, ")\n"))
cat(paste0("  n_clust  = ", N_CLUST, "\n"))
cat(paste0("  problem  = ", PROBLEM, "\n"))
cat(paste0("  Date     : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste0("  Host     : ", Sys.info()["nodename"], "\n"))
cat(paste0("  PID      : ", Sys.getpid(), "\n\n"))

# Création des répertoires
dir_models_mt      <- file.path(base_path, "Models_MT")
dir_predictions_mt <- file.path(base_path, "Predictions_MT")
dir_run_info       <- file.path(base_path, "run_info")

for (d in c(dir_models_mt, dir_predictions_mt, dir_run_info)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

run_info <- list(
  script = "Benchmark_NeurIPS_MT_nclust.R",
  n_out = N_OUT, n_train = N_TRAIN, n_pred = N_PRED, n_clust = N_CLUST,
  problem = PROBLEM, seed = SEED,
  started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  hostname = as.character(Sys.info()["nodename"]),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "N/A"),
  r_version = R.version.string, pid = Sys.getpid()
)
run_info_file <- file.path(dir_run_info,
  paste0("mt_seed", SEED, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"))
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

  data_raw       <- datasets$data_raw
  train_task_ids <- datasets$train_task_ids
  test_task_ids  <- datasets$test_task_ids

  # Reconstruire les splits pred/test à partir des données brutes
  pred_tasks_data_raw <- list()
  test_tasks_data_raw <- list()

  if (PROBLEM == "interpolation") {
    gap_start <- -0.5; gap_end <- 0.5
    train_data_raw <- data_raw %>%
      dplyr::filter(Task_ID %in% train_task_ids) %>%
      dplyr::filter(Output_ID != "1" |
                    (Output_ID == "1" & !(Input >= gap_start & Input <= gap_end)))

    for (tid in test_task_ids) {
      task_full <- data_raw %>% dplyr::filter(Task_ID == tid)
      task_O1   <- task_full %>% dplyr::filter(Output_ID == "1") %>% dplyr::arrange(Input)
      test_tasks_data_raw[[tid]] <- task_O1 %>%
        dplyr::filter(Input >= gap_start & Input <= gap_end)
      pred_tasks_data_raw[[tid]] <- dplyr::bind_rows(
        task_full %>% dplyr::filter(Output_ID != "1"),
        task_O1 %>% dplyr::filter(!(Input >= gap_start & Input <= gap_end))
      ) %>% dplyr::arrange(Output_ID, Input)
    }
  } else {
    train_data_raw <- data_raw %>% dplyr::filter(Task_ID %in% train_task_ids)

    for (tid in test_task_ids) {
      task_full <- data_raw %>% dplyr::filter(Task_ID == tid)
      task_O1   <- task_full %>% dplyr::filter(Output_ID == "1") %>% dplyr::arrange(Input)
      test_tasks_data_raw[[tid]] <- task_O1 %>%
        dplyr::filter(Input >= 0 & Input <= 1)
      pred_tasks_data_raw[[tid]] <- dplyr::bind_rows(
        task_full %>% dplyr::filter(Output_ID != "1"),
        task_O1 %>% dplyr::filter(!(Input >= 0 & Input <= 1))
      ) %>% dplyr::arrange(Output_ID, Input)
    }
  }

  # --- 3. REFORMATAGE POUR MT ---
  train_data_mt <- train_data_raw %>%
    dplyr::select(-any_of(c("Cluster_ID", "Dataframe"))) %>%
    dplyr::mutate(
      Task_ID   = paste0("T", Task_ID, "_O", Output_ID),
      Output_ID = "1"
    )

  cat(paste0("  MT train tasks: ", length(unique(train_data_mt$Task_ID)), "\n"))

  output_ids <- unique(data_raw$Output_ID)

  # ==================================================================
  # BRANCHEMENT SELON N_CLUST
  # ==================================================================

  if (N_CLUST == 1) {
    # ================================================================
    # N_CLUST = 1 : Magma (sans clustering) pour MT
    # ================================================================

    # --- HP initiaux extraits du modèle MOMT (convolution → SE) ---
    ini_hp_0_mt <- trained_model_momt$hp_t %>%
      dplyr::filter(Output_ID == "1") %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::select(-any_of(c("l_u_t", "Task_ID"))) %>%
      dplyr::mutate(se_lengthscale = l_t, se_variance = S_t) %>%
      dplyr::select(-c(l_t, S_t))

    ini_hp_t_mt <- ini_hp_0_mt

    prior_means_mt <- mean(train_data_mt$Output, na.rm = TRUE)

    # --- Entraînement train_magma ---
    t_train_start <- Sys.time()
    trained_model_mt <- train_magma(
      data             = train_data_mt,
      ini_hp_0         = ini_hp_0_mt %>% dplyr::select(-noise),
      kern_0           = "SE",
      ini_hp_t         = ini_hp_t_mt,
      kern_t           = "SE",
      shared_hp_tasks  = TRUE,
      prior_mean       = prior_means_mt,
      pen_diag         = 1e-10
    )
    duration_train <- as.numeric(difftime(Sys.time(), t_train_start, units = "secs"))
    cat(paste0("  Training (Magma MT): ", round(duration_train, 1), "s\n"))

    # Sauvegarde modèle MT
    model_mt_to_save <- list(
      trained_model = trained_model_mt,
      t_training    = duration_train,
      seed          = SEED,
      n_train       = N_TRAIN,
      n_clust       = N_CLUST
    )
    saveRDS(model_mt_to_save, file.path(dir_models_mt, paste0("model_seed_", SEED, ".rds")))

    # --- Hyperpostérieur global ---
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
        if (oid == "1") {
          test_data_oid <- test_tasks_data_raw[[tid]]
          if (nrow(test_data_oid) > 0) {
            all_test_inputs_global[[mt_task_id]] <- test_data_oid %>%
              dplyr::mutate(Output_ID = "1") %>%
              dplyr::select(Input, Input_ID, Output_ID)
          }
        }
      }
    }

    grid_hyperpost_global <- dplyr::bind_rows(
      dplyr::bind_rows(all_pred_inputs_global),
      dplyr::bind_rows(all_test_inputs_global)
    ) %>% dplyr::distinct() %>% dplyr::arrange(Output_ID, Input_ID, Input)

    t_hp_start <- Sys.time()
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
    duration_hyperpost <- as.numeric(difftime(Sys.time(), t_hp_start, units = "secs"))
    cat(paste0("  Hyperpost: ", round(duration_hyperpost, 1), "s\n"))

    # --- Prédictions par (tâche, output) ---
    all_predictions <- list()
    t_pred_total    <- 0

    for (test_task_id in test_task_ids) {
      task_predictions <- list()

      for (oid in output_ids) {
        if (oid != "1") next

        pred_data_oid <- pred_tasks_data_raw[[test_task_id]] %>%
          dplyr::filter(Output_ID == oid) %>%
          dplyr::select(-any_of(c("Task_ID", "Cluster_ID", "Dataframe"))) %>%
          dplyr::mutate(Output_ID = as.factor("1"))

        if (nrow(pred_data_oid) == 0) next

        test_data_oid <- test_tasks_data_raw[[test_task_id]]
        if (nrow(test_data_oid) == 0) next

        test_grid_inputs <- test_data_oid %>%
          dplyr::mutate(Output_ID = as.factor("1")) %>%
          dplyr::select(Input, Input_ID, Output_ID) %>%
          dplyr::distinct()

        t_pred_start <- Sys.time()
        pred_res <- pred_magma(
          data          = pred_data_oid,
          trained_model = trained_model_mt,
          grid_inputs   = test_grid_inputs,
          kern          = "SE",
          hyperpost     = hyperpost_global,
          get_full_cov  = TRUE,
          plot          = FALSE
        )
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
        t_pred = sum(sapply(task_predictions, function(x) x$t_pred))
      )
    }

  } else {
    # ================================================================
    # N_CLUST >= 2 : MagmaClust pour MT
    # ================================================================

    # En MT, on crée N_CLUST * N_OUT clusters (un par (cluster MOMT, output d'origine))
    N_CLUST_MT <- N_CLUST * N_OUT
    output_ids_sorted <- sort(output_ids)
    cat(paste0("  N_CLUST_MT = ", N_CLUST_MT, " (", N_CLUST, " × ", N_OUT, ")\n"))

    # --- HP initiaux extraits du modèle MOMT (étendus à k*n_out clusters) ---
    ini_hp_k <- trained_model_momt$hp_k %>%
      dplyr::select(-any_of("l_u_t")) %>%
      dplyr::mutate(
        se_lengthscale = l_t,
        se_variance    = S_t,
          Cluster_ID = (as.integer(as.character(Cluster_ID)) - 1) * length(output_ids_sorted) +
                     match(as.character(Output_ID), output_ids_sorted),
        Output_ID  = as.factor("1")
      ) %>%
      dplyr::select(-c(l_t, S_t)) %>%
      dplyr::arrange(Cluster_ID)

    ini_hp_t <- trained_model_momt$hp_t %>%
      dplyr::filter(Output_ID == "1") %>%
      dplyr::filter(Task_ID %in% unique(datasets$train_data$Task_ID)) %>%
      dplyr::select(-any_of(c("l_u_t", "Task_ID"))) %>%
      dplyr::mutate(se_lengthscale = l_t, se_variance = S_t) %>%
      dplyr::select(-c(l_t, S_t)) %>%
      dplyr::slice(1)

    # Prior means par cluster MT
    cluster_mapping_momt <- datasets$cluster_mapping %>%
        dplyr::mutate(
          Task_ID = as.character(Task_ID),
          Cluster_ID = as.integer(as.character(Cluster_ID))
        )

    train_data_mt_with_cluster <- train_data_mt %>%
      dplyr::mutate(
        Task_ID_orig   = sub("^T(.+)_O.*$", "\\1", Task_ID),
        Output_ID_orig = sub("^T.+_O(.+)$", "\\1", Task_ID)
      ) %>%
      dplyr::left_join(cluster_mapping_momt, by = c("Task_ID_orig" = "Task_ID")) %>%
      dplyr::mutate(
          Cluster_ID = (as.integer(as.character(Cluster_ID)) - 1) * length(output_ids_sorted) +
                     match(Output_ID_orig, output_ids_sorted)
      )

    prior_means_mt <- train_data_mt_with_cluster %>%
      dplyr::group_by(Cluster_ID) %>%
      dplyr::summarise(prior_mean = mean(Output, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(Cluster_ID) %>%
      dplyr::pull(prior_mean)

    # ini_mixture : étendre de N_TRAIN (k clusters) à N_TRAIN × N_OUT (k*n_out clusters)
    ini_mixture_momt <- trained_model_momt$ini_args$ini_mixture

    ini_mixture_mt_rows <- list()
    for (r in 1:nrow(ini_mixture_momt)) {
      orig_tid <- as.character(ini_mixture_momt$Task_ID[r])
      for (oid in output_ids) {
        oid_idx <- match(oid, output_ids_sorted)
        new_row <- tibble(Task_ID = paste0("T", orig_tid, "_O", oid))
        for (c_momt in 1:N_CLUST) {
          p_momt <- ini_mixture_momt[[paste0("Cluster_", c_momt)]][r]
          for (o_idx in seq_along(output_ids_sorted)) {
            c_mt <- (c_momt - 1) * length(output_ids_sorted) + o_idx
            new_row[[paste0("Cluster_", c_mt)]] <- if (o_idx == oid_idx) p_momt else 0
          }
        }
        ini_mixture_mt_rows[[length(ini_mixture_mt_rows) + 1]] <- new_row
      }
    }
    ini_mixture_mt <- dplyr::bind_rows(ini_mixture_mt_rows)

    # --- Entraînement MagmaClust MT ---
    t_train_start <- Sys.time()
    trained_model_mt <- train_magmaclust(
      data             = train_data_mt,
      ini_mixture      = ini_mixture_mt,
      ini_hp_k         = ini_hp_k %>% dplyr::select(-noise),
      nb_cluster       = N_CLUST_MT,
      kern_k           = "SE",
      ini_hp_t         = ini_hp_t,
      kern_t           = "SE",
      shared_hp_tasks  = TRUE,
      shared_hp_clusts = TRUE,
      prior_mean_k     = prior_means_mt,
      pen_diag         = 1e-10
    )
    duration_train <- as.numeric(difftime(Sys.time(), t_train_start, units = "secs"))
    cat(paste0("  Training (MagmaClust MT): ", round(duration_train, 1), "s\n"))

    # Sauvegarde modèle MT
    model_mt_to_save <- list(
      trained_model = trained_model_mt,
      t_training    = duration_train,
      seed          = SEED,
      n_train       = N_TRAIN,
      n_clust       = N_CLUST
    )
    saveRDS(model_mt_to_save, file.path(dir_models_mt, paste0("model_seed_", SEED, ".rds")))

    # --- Hyperpostérieur global ---
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
        if (oid == "1") {
          test_data_oid <- test_tasks_data_raw[[tid]]
          if (nrow(test_data_oid) > 0) {
            all_test_inputs_global[[mt_task_id]] <- test_data_oid %>%
              dplyr::mutate(Output_ID = "1") %>%
              dplyr::select(Input, Input_ID, Output_ID)
          }
        }
      }
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

    # --- Prédictions par (tâche, output) ---
    all_predictions <- list()
    t_pred_total    <- 0

    for (test_task_id in test_task_ids) {
      task_predictions <- list()

      for (oid in output_ids) {
        if (oid != "1") next

        pred_data_oid <- pred_tasks_data_raw[[test_task_id]] %>%
          dplyr::filter(Output_ID == oid) %>%
          dplyr::select(-any_of(c("Task_ID", "Cluster_ID", "Dataframe"))) %>%
          dplyr::mutate(Output_ID = as.factor("1"))

        if (nrow(pred_data_oid) == 0) next

        test_data_oid <- test_tasks_data_raw[[test_task_id]]
        if (nrow(test_data_oid) == 0) next

        test_grid_inputs <- test_data_oid %>%
          dplyr::mutate(Output_ID = as.factor("1")) %>%
          dplyr::select(Input, Input_ID, Output_ID) %>%
          dplyr::distinct()

        t_pred_start <- Sys.time()
        pred_res <- pred_magmaclust(
          data          = pred_data_oid,
          trained_model = trained_model_mt,
          grid_inputs   = test_grid_inputs,
          kern          = "SE",
          hyperpost     = hyperpost_global,
          get_full_cov  = TRUE,
          plot          = FALSE
        )
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
        t_pred = sum(sapply(task_predictions, function(x) x$t_pred))
      )
    }
  }

  # --- SAUVEGARDE ---
  predictions_to_save <- list(
    predictions        = all_predictions,
    t_hyperpost_global = duration_hyperpost,
    t_pred_total       = t_pred_total,
    t_training         = duration_train,
    seed               = SEED,
    n_clust            = N_CLUST,
    test_task_ids      = test_task_ids,
    n_train            = N_TRAIN,
    n_pred             = N_PRED
  )
  saveRDS(predictions_to_save,
          file.path(dir_predictions_mt, paste0("predictions_seed_", SEED, ".rds")))

  cat(paste0("  Prédictions OK (", round(t_pred_total, 1), "s total)\n"))

  rm(trained_model_mt, hyperpost_global, all_predictions)
  gc()

}, error = function(e) {
  cat(paste0("\n!!! ERREUR [MT out=", N_OUT, " train=", N_TRAIN,
             " pred=", N_PRED, " clust=", N_CLUST, " ", PROBLEM,
             " seed=", SEED, "] : ", e$message, "\n"))
  if (!is.null(e$parent)) cat("  Cause: ", conditionMessage(e$parent), "\n")
})

# --- Fin ---
duration_total <- as.numeric(difftime(Sys.time(), t_global_start, units = "mins"))
cat(paste0("\n=== MT TERMINÉ [out=", N_OUT, " train=", N_TRAIN,
           " pred=", N_PRED, " clust=", N_CLUST, " ", PROBLEM,
           " seed=", SEED, "] ===\n"))
cat(paste0("Durée : ", round(duration_total, 2), " min\n"))

run_info$ended_at      <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_info$duration_mins <- round(duration_total, 2)
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)
