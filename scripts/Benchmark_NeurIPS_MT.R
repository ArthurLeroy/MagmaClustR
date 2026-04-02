# ==========================================================================
# Benchmark_NeurIPS_MT.R
# NeurIPS : Benchmark MT (MagmaClust single-output, chaque (Task,Output) → tâche)
#
# Usage :
#   Rscript Benchmark_NeurIPS_MT.R --n_out=2 --n_train=30 --n_pred=1 \
#     --problem=interpolation --seed=1
#
# Charge les données générées par MOMT (données BRUTES, non scalées).
# Réorganise : chaque (Task_ID, Output_ID) de MOMT → une "tâche" MT avec
#   Output_ID="1", kernel SE.
# Nb tâches train MT = N_TRAIN × N_OUT
# Nb tâches pred MT  = N_PRED × N_OUT
# Pour la prédiction : un pred_magmaclust() par (tâche pred, output).
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

N_CLUSTERS <- 3
N_CLUSTERS_MT <- N_OUT * N_CLUSTERS

# --- 1. SETUP & LIBRARIES ---
username  <- Sys.getenv("USER")
pkg_dir   <- file.path("/scratch", username, "MagmaClustR")
base_path <- file.path("/scratch", username, "NeurIPS_experiments",
                        paste0("out", N_OUT, "_train", N_TRAIN, "_pred", N_PRED),
                        PROBLEM)

setwd(pkg_dir)

library(tidyverse)
library(readr)
library(MagmaClustR)

source(file.path("/scratch", username, "NeurIPS_experiments", "utils", "generate_random_config.R"))

config_label <- paste0("out", N_OUT, "_train", N_TRAIN, "_pred", N_PRED)

cat(paste0("=== NeurIPS MT | ", config_label, " | ", PROBLEM, " | seed=", SEED, " ===\n"))
cat(paste0("  n_out    = ", N_OUT, " (→ MT tasks = n_tasks × n_out)\n"))
cat(paste0("  n_train  = ", N_TRAIN, " (→ MT train tasks = ", N_TRAIN * N_OUT, ")\n"))
cat(paste0("  n_pred   = ", N_PRED, " (→ MT pred tasks = ", N_PRED * N_OUT, ")\n"))
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
  script = "Benchmark_NeurIPS_MT.R",
  n_out = N_OUT, n_train = N_TRAIN, n_pred = N_PRED,
  problem = PROBLEM, seed = SEED,
  n_clusters = N_CLUSTERS,
  n_clusters_mt = N_CLUSTERS_MT,
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

  # Données BRUTES (non scalées) pour MT
  data_raw            <- datasets$data_raw
  train_task_ids      <- datasets$train_task_ids
  test_task_ids       <- datasets$test_task_ids

  # Reconstruire pred/test à partir des données brutes
  # On refait le split identique à MOMT mais sur data_raw
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
  # Chaque (Task_ID, Output_ID) → une "tâche" MT avec Task_ID = "T{task}_O{out}"
  # et Output_ID = "1"

  # 3a. Données d'entraînement
  train_data_mt <- train_data_raw %>%
    dplyr::select(-any_of("Cluster_ID")) %>%
    dplyr::mutate(
      Task_ID   = paste0("T", Task_ID, "_O", Output_ID),
      Output_ID = "1"
    )

  cat(paste0("  MT train tasks: ", length(unique(train_data_mt$Task_ID)), "\n"))

  # 3b. ini_mixture fraîche pour MT (N_CLUSTERS_MT = N_OUT × N_CLUSTERS)
  # En MT, outputs et tâches ont le même statut → plus de sens de réutiliser
  # l'ini_mixture MOMT. On en calcule une nouvelle sur les données MT.
  ini_mixture_mt <- ini_mixture(data = train_data_mt, k = N_CLUSTERS_MT)

  cat(paste0("  ini_mixture MT : ", nrow(ini_mixture_mt), " tâches × ",
             N_CLUSTERS_MT, " clusters\n"))

  # Assignation dure : dans quel cluster MT chaque tâche tombe-t-elle ?
  mt_cluster_assignment <- ini_mixture_mt %>%
    tidyr::pivot_longer(cols = starts_with("K"), names_to = "MT_Cluster",
                        values_to = "Probability") %>%
    dplyr::group_by(Task_ID) %>%
    dplyr::slice_max(Probability, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(Task_ID, MT_Cluster)

  # 3c. Mapping MT cluster → (MOMT Cluster, Output_ID original)
  # Pour chaque tâche MT "T{tid}_O{oid}", on connaît le cluster MOMT de tid
  cluster_mapping_momt <- datasets$cluster_mapping %>%
    dplyr::mutate(Task_ID = as.character(Task_ID))

  mt_task_origins <- mt_cluster_assignment %>%
    dplyr::mutate(
      Task_ID_orig = sub("^T(.+)_O.*$", "\\1", Task_ID),
      orig_Output_ID = sub("^T.+_O(.+)$", "\\1", Task_ID)
    ) %>%
    dplyr::left_join(cluster_mapping_momt, by = c("Task_ID_orig" = "Task_ID"))

  # Cluster MT dominant → (MOMT_Cluster, Output_ID) le plus fréquent
  mt_to_momt_mapping <- mt_task_origins %>%
    dplyr::count(MT_Cluster, Cluster_ID, orig_Output_ID) %>%
    dplyr::group_by(MT_Cluster) %>%
    dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(MT_Cluster, Cluster_ID, orig_Output_ID)

  cat("  Mapping MT cluster → MOMT (Cluster, Output) :\n")
  print(mt_to_momt_mapping)

  # 3d. Extraction HP initiaux depuis le modèle MOMT via le mapping
  momt_hp_k <- trained_model_momt$hp_k %>%
    dplyr::mutate(Output_ID = as.character(Output_ID),
                  Cluster_ID = as.character(Cluster_ID))

  ini_hp_k <- mt_to_momt_mapping %>%
    dplyr::left_join(momt_hp_k,
                     by = c("Cluster_ID" = "Cluster_ID",
                            "orig_Output_ID" = "Output_ID")) %>%
    dplyr::transmute(
      Cluster_ID = MT_Cluster,
      Output_ID = "1",
      se_lengthscale = l_t,
      se_variance = S_t,
      noise = -3
    )

  cat("  ini_hp_k MT :\n"); print(ini_hp_k)

  ini_hp_t <- trained_model_momt$hp_t %>%
    dplyr::filter(Output_ID == "1") %>%
    dplyr::filter(Task_ID %in% unique(datasets$train_data$Task_ID)) %>%
    dplyr::select(-any_of(c("l_u_t", "Task_ID"))) %>%
    dplyr::mutate(se_lengthscale = l_t, se_variance = S_t) %>%
    dplyr::select(-c(l_t, S_t)) %>%
    dplyr::slice(1)

  # 3e. Prior means MT : moyenne empirique par cluster MT (sur données brutes)
  train_data_mt_with_cluster <- train_data_mt %>%
    dplyr::left_join(mt_cluster_assignment, by = "Task_ID")

  prior_means_mt <- train_data_mt_with_cluster %>%
    dplyr::group_by(MT_Cluster) %>%
    dplyr::summarise(prior_mean = mean(Output, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(MT_Cluster) %>%
    dplyr::pull(prior_mean)

  cat(paste0("  Prior means MT (", length(prior_means_mt), " clusters): ",
             paste(round(prior_means_mt, 4), collapse = ", "), "\n"))

  # --- 4. ENTRAÎNEMENT MAGMACLUST SINGLE-OUTPUT ---
  t_train_start <- Sys.time()
  trained_model_mt <- train_magmaclust(
    data             = train_data_mt,
    ini_mixture      = ini_mixture_mt,
    ini_hp_k         = ini_hp_k,
    nb_cluster       = N_CLUSTERS_MT,
    kern_k           = "SE",
    ini_hp_t         = ini_hp_t,
    kern_t           = "SE",
    shared_hp_tasks  = TRUE,
    shared_hp_clusts = TRUE,
    prior_mean_k     = prior_means_mt,
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

  # --- 5. HYPERPOSTÉRIEUR GLOBAL (une seule fois pour toutes les prédictions) ---
  # Collecter tous les inputs nécessaires pour l'hyperpost global
  all_pred_inputs_global <- list()
  all_test_inputs_global <- list()

  output_ids <- unique(data_raw$Output_ID)

  for (tid in test_task_ids) {
    for (oid in output_ids) {
      mt_task_id <- paste0("T", tid, "_O", oid)

      # Contexte de cette tâche MT
      pred_data_oid <- pred_tasks_data_raw[[tid]] %>%
        dplyr::filter(Output_ID == oid)
      if (nrow(pred_data_oid) > 0) {
        all_pred_inputs_global[[mt_task_id]] <- pred_data_oid %>%
          dplyr::mutate(Output_ID = "1") %>%
          dplyr::select(Output_ID, Input_ID, Input)
      }

      # Points test (seulement Output 1)
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

  # --- 6. PRÉDICTION PAR (TÂCHE, OUTPUT) ---
  # Pour chaque tâche test, on prédit chaque output séparément
  # Seul l'output 1 nous intéresse pour les métriques, mais on prédit tous les
  # outputs pour être complet
  all_predictions <- list()
  t_pred_total <- 0

  for (test_task_id in test_task_ids) {
    task_predictions <- list()

    for (oid in output_ids) {
      mt_task_id <- paste0("T", test_task_id, "_O", oid)

      # Données de contexte pour cette (tâche, output)
      pred_data_oid <- pred_tasks_data_raw[[test_task_id]] %>%
        dplyr::filter(Output_ID == oid) %>%
        dplyr::select(-any_of(c("Task_ID", "Cluster_ID"))) %>%
        dplyr::mutate(Output_ID = as.factor("1"))

      if (nrow(pred_data_oid) == 0) next

      # Grid d'inputs pour la prédiction
      if (oid == "1") {
        # Pour output 1 : prédire sur les points test
        test_data_oid <- test_tasks_data_raw[[test_task_id]]
        if (nrow(test_data_oid) == 0) next

        test_grid_inputs <- test_data_oid %>%
          dplyr::mutate(Output_ID = as.factor("1")) %>%
          dplyr::select(Input, Input_ID, Output_ID) %>%
          dplyr::distinct()
      } else {
        # Pour les autres outputs : pas de points test, on skip la prédiction
        next
      }

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
      t_pred_task <- as.numeric(difftime(Sys.time(), t_pred_start, units = "secs"))
      t_pred_total <- t_pred_total + t_pred_task

      task_predictions[[oid]] <- list(
        prediction = pred_res,
        t_pred     = t_pred_task
      )
    }

    # Stocker la vérité (données brutes, non scalées) pour les métriques
    all_predictions[[test_task_id]] <- list(
      predictions_by_output = task_predictions,
      truth                 = test_tasks_data_raw[[test_task_id]],
      inputs_pred           = pred_tasks_data_raw[[test_task_id]],
      t_pred                = sum(sapply(task_predictions, function(x) x$t_pred))
    )
  }

  # --- 7. SAUVEGARDE ---
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
          file.path(dir_predictions_mt, paste0("predictions_seed_", SEED, ".rds")))

  cat(paste0("  Prédictions OK (", round(t_pred_total, 1), "s total)\n"))

  rm(trained_model_mt, hyperpost_global, all_predictions)
  gc()

}, error = function(e) {
  cat(paste0("\n!!! ERREUR [MT out=", N_OUT, " train=", N_TRAIN,
             " pred=", N_PRED, " ", PROBLEM, " seed=", SEED, "] : ",
             e$message, "\n"))
  if (!is.null(e$parent)) cat("  Cause: ", conditionMessage(e$parent), "\n")
})

# --- Fin ---
duration_total <- as.numeric(difftime(Sys.time(), t_global_start, units = "mins"))
cat(paste0("\n=== MT TERMINÉ [out=", N_OUT, " train=", N_TRAIN,
           " pred=", N_PRED, " ", PROBLEM, " seed=", SEED, "] ===\n"))
cat(paste0("Durée : ", round(duration_total, 2), " min\n"))

run_info$ended_at      <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_info$duration_mins <- round(duration_total, 2)
writeLines(jsonlite::toJSON(run_info, auto_unbox = TRUE, pretty = TRUE), run_info_file)
