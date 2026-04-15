# ==========================================================================
# Benchmark_NeurIPS_MT_nclust_forecasting.R
# NeurIPS : Benchmark MT (FORECASTING ONLY) avec paramètre n_clust variable
# Corrige le bug de la version originale : on enlève les observations sur
# l'intervalle de forecasting [0,1] pour TOUS les outputs (pas seulement Output 1).
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
stopifnot(PROBLEM == "forecasting")
stopifnot(SEED >= 1 & SEED <= 5)

# --- 1. SETUP & LIBRARIES ---
username  <- Sys.getenv("USER")
pkg_dir   <- file.path("/scratch", username, "MagmaClustR")
base_root <- file.path("/scratch", username, "NeurIPS_experiments_forecasting")

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
dir_models_mt      <- file.path(base_root, "Models_MT", config_label)
dir_predictions_mt <- file.path(base_root, "Predictions_MT", config_label)
dir_run_info       <- file.path(base_root, "run_info", config_label)

for (d in c(dir_models_mt, dir_predictions_mt, dir_run_info)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

run_info <- list(
  script = "Benchmark_NeurIPS_MT_nclust_forecasting.R",
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
  dir_datasets    <- file.path(base_root, "Datasets", config_label)
  dir_models_momt <- file.path(base_root, "Models_MOMT", config_label)

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

  # --- Forecasting : on enlève [0,1] pour TOUS les outputs des tâches test ---
  train_data_raw <- data_raw %>% dplyr::filter(Task_ID %in% train_task_ids)

  for (tid in test_task_ids) {
    task_full <- data_raw %>% dplyr::filter(Task_ID == tid)
    test_tasks_data_raw[[tid]] <- task_full %>%
      dplyr::filter(Input >= 0 & Input <= 1) %>%
      dplyr::arrange(Output_ID, Input)
    pred_tasks_data_raw[[tid]] <- task_full %>%
      dplyr::filter(!(Input >= 0 & Input <= 1)) %>%
      dplyr::arrange(Output_ID, Input)
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

    # --- HP initialisation via GP vanille (SE) ---
    data_for_init <- train_data_mt %>%
      dplyr::mutate(Input = round(Input, 6))

    best_task_id_init <- data_for_init %>%
      dplyr::count(Task_ID) %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      dplyr::slice(1) %>%
      dplyr::pull(Task_ID)

    sub_data_agg <- data_for_init %>%
      dplyr::filter(Task_ID == best_task_id_init)

    mean_emp <- mean(sub_data_agg$Output, na.rm = TRUE)
    mean_vec <- rep(mean_emp, nrow(sub_data_agg))

    best_ll    <- +Inf
    best_hp_gp <- NULL

    for (seed_retry in 1:10) {
      tryCatch({
        set.seed(seed_retry * 1000)
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
            pen_diag = 1e-10,
            hp_col_names = c("se_variance", "se_lengthscale")
          )
        }, error = function(e) return(+Inf))

        if (ll_val < best_ll) {
          best_ll    <- ll_val
          best_hp_gp <- hp_tmp
        }
      }, error = function(e) {
        cat(paste0("  [train_gp ERR] retry=", seed_retry, ": ", e$message, "\n"))
      })
    }

    if (is.null(best_hp_gp)) {
      best_hp_gp <- list(se_lengthscale = 0, se_variance = 0, noise = -3)
    }

    ini_hp_0_mt <- tibble(
      Output_ID      = as.factor("1"),
      se_lengthscale = best_hp_gp$se_lengthscale,
      se_variance    = best_hp_gp$se_variance,
      noise          = -3
    )

    ini_hp_t_mt <- ini_hp_0_mt

    prior_means_mt <- mean_emp

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

    # --- HP initialisation via GP vanille (SE) par cluster ---
    cluster_mapping_momt <- datasets$cluster_mapping %>%
      dplyr::mutate(Task_ID = as.character(Task_ID))

    train_data_mt_with_cluster <- train_data_mt %>%
      dplyr::mutate(Task_ID_orig = sub("^T(.+)_O.*$", "\\1", Task_ID)) %>%
      dplyr::left_join(cluster_mapping_momt, by = c("Task_ID_orig" = "Task_ID"))

    clusters_ids <- sort(unique(cluster_mapping_momt$Cluster_ID))
    hp_k_extracted_list <- list()

    for (k_id in clusters_ids) {
      data_for_init_k <- train_data_mt_with_cluster %>%
        dplyr::filter(Cluster_ID == k_id) %>%
        dplyr::mutate(Input = round(Input, 6))

      best_task_id_k <- data_for_init_k %>%
        dplyr::count(Task_ID) %>%
        dplyr::arrange(dplyr::desc(n)) %>%
        dplyr::slice(1) %>%
        dplyr::pull(Task_ID)

      sub_data_agg <- data_for_init_k %>%
        dplyr::filter(Task_ID == best_task_id_k)

      mean_emp_k <- mean(sub_data_agg$Output, na.rm = TRUE)
      mean_vec_k <- rep(mean_emp_k, nrow(sub_data_agg))

      best_ll    <- +Inf
      best_hp_gp <- NULL

      for (seed_retry in 1:10) {
        tryCatch({
          set.seed(seed_retry * 1000 + which(clusters_ids == k_id))
          hp_tmp <- suppressWarnings(suppressMessages(
            train_gp(data = sub_data_agg, kern = "SE",
                     prior_mean = mean_emp_k, ini_hp = NULL)
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
              mean = mean_vec_k, kern = "SE", post_cov = 0,
              pen_diag = 1e-10,
              hp_col_names = c("se_variance", "se_lengthscale")
            )
          }, error = function(e) return(+Inf))

          if (ll_val < best_ll) {
            best_ll    <- ll_val
            best_hp_gp <- hp_tmp
          }
        }, error = function(e) {
          cat(paste0("  [train_gp ERR] k=", k_id,
                     " retry=", seed_retry, ": ", e$message, "\n"))
        })
      }

      if (is.null(best_hp_gp)) {
        cat(paste0("  [WARN] Aucun train_gp OK pour k=", k_id, "\n"))
        best_hp_gp <- list(se_lengthscale = 0, se_variance = 0, noise = -3)
      }

      hp_k_extracted_list[[length(hp_k_extracted_list) + 1]] <- tibble(
        Cluster_ID     = k_id,
        Output_ID      = as.factor("1"),
        se_lengthscale = best_hp_gp$se_lengthscale,
        se_variance    = best_hp_gp$se_variance,
        noise          = -3
      )
    }

    ini_hp_k <- bind_rows(hp_k_extracted_list)

    ini_hp_t <- ini_hp_k %>%
      dplyr::filter(Cluster_ID == clusters_ids[[1]]) %>%
      dplyr::select(-Cluster_ID)

    prior_means_mt <- train_data_mt_with_cluster %>%
      dplyr::group_by(Cluster_ID) %>%
      dplyr::summarise(prior_mean = mean(Output, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(Cluster_ID) %>%
      dplyr::pull(prior_mean)

    # ini_mixture : étendre de N_TRAIN à N_TRAIN × N_OUT
    ini_mixture_momt <- trained_model_momt$ini_args$ini_mixture

    ini_mixture_mt_rows <- list()
    for (r in 1:nrow(ini_mixture_momt)) {
      orig_tid <- as.character(ini_mixture_momt$Task_ID[r])
      for (oid in output_ids) {
        new_row <- ini_mixture_momt[r, ]
        new_row$Task_ID <- paste0("T", orig_tid, "_O", oid)
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
      nb_cluster       = N_CLUST,
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
