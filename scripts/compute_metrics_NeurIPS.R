# ==========================================================================
# compute_metrics_NeurIPS.R
# NeurIPS : Calcul des métriques à partir des prédictions sauvegardées
#
# Usage :
#   Rscript compute_metrics_NeurIPS.R
#
# Lit toutes les prédictions (MOMT, MO, MT) pour chaque configuration,
# calcule RMSE, NLL, coverage 95%, temps, et produit les tableaux récapitulatifs.
#
# IMPORTANT : Les métriques sont calculées sur OUTPUT 1 uniquement,
# en espace original (dé-scalé pour MOMT et MO).
#
# Sortie :
#   - metrics_raw_NeurIPS.csv     : métriques brutes par config/model/seed
#   - metrics_summary_NeurIPS.csv : (mean±std) et (median±IQR) par config/model
# ==========================================================================

library(tidyverse)
library(matrixStats)

# ---- Couverture de l'IC ----
compute_coverage <- function(pred_mean, pred_var, truth, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  pred_sd <- sqrt(pmax(pred_var, 1e-10))
  lower <- pred_mean - z * pred_sd
  upper <- pred_mean + z * pred_sd
  mean(truth >= lower & truth <= upper)
}

# ---- NLL pour mixture de gaussiennes (MOMT et MT) ----
compute_nll_mixture <- function(pred_by_cluster, weights, truth) {
  if ("Input_1" %in% names(pred_by_cluster) && !"Input" %in% names(pred_by_cluster)) {
    pred_by_cluster <- pred_by_cluster %>% dplyr::rename(Input = Input_1)
  }

  required_cols <- c("Input", "Output_ID", "Mean", "Var", "Cluster")
  if (!all(required_cols %in% names(pred_by_cluster))) return(NA)

  pred_by_cluster <- pred_by_cluster %>%
    dplyr::mutate(Input = round(Input, 5), Output_ID = as.character(Output_ID)) %>%
    dplyr::distinct(Input, Output_ID, Cluster, .keep_all = TRUE)

  clusters <- unique(pred_by_cluster$Cluster)
  K <- length(clusters)
  if (is.null(weights)) {
    weights <- rep(1 / K, K)
    names(weights) <- clusters
  }
  if (is.null(names(weights))) names(weights) <- clusters
  weights <- weights / sum(weights)

  truth_clean <- truth %>%
    dplyr::mutate(Input = round(Input, 5), Output_ID = as.character(Output_ID)) %>%
    dplyr::distinct(Input, Output_ID, .keep_all = TRUE)

  pred_with_truth <- pred_by_cluster %>%
    dplyr::inner_join(truth_clean, by = c("Input", "Output_ID"))

  if (nrow(pred_with_truth) == 0) return(NA)

  nll_per_point <- pred_with_truth %>%
    dplyr::mutate(
      Var_safe = pmax(Var, 1e-10),
      log_density = dnorm(Output, mean = Mean, sd = sqrt(Var_safe), log = TRUE),
      w = weights[as.character(Cluster)]
    ) %>%
    dplyr::group_by(Input, Output_ID) %>%
    dplyr::summarise(
      log_mixture_density = matrixStats::logSumExp(log(w) + log_density),
      .groups = "drop"
    ) %>%
    dplyr::mutate(nll = -log_mixture_density)

  return(mean(nll_per_point$nll))
}

# ---- NLL simple gaussien ----
compute_nll_gaussian <- function(pred_mean, pred_var, truth) {
  pred_var_safe <- pmax(pred_var, 1e-12)
  mean(0.5 * log(2 * pi * pred_var_safe) + 0.5 * (truth - pred_mean)^2 / pred_var_safe)
}

# ---- Extraction des métriques pour MOMT ----
extract_metrics_momt <- function(pred_entry, scale_factors) {
  tryCatch({
    pred_res   <- pred_entry$prediction
    truth_data <- pred_entry$truth  # données scalées
    cluster_id <- pred_entry$cluster_id

    # Filtrer Output 1 seulement pour les métriques
    truth_o1 <- truth_data %>% dplyr::filter(Output_ID == "1")
    if (nrow(truth_o1) == 0) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))

    # Extraire la prédiction agrégée (mixture)
    if (is.data.frame(pred_res) || tibble::is_tibble(pred_res)) {
      pred_df <- pred_res
    } else if ("mixture_pred" %in% names(pred_res)) {
      pred_df <- pred_res$mixture_pred
    } else if ("mixture" %in% names(pred_res) && is.data.frame(pred_res$mixture)) {
      pred_df <- pred_res$mixture
    } else if ("pred" %in% names(pred_res)) {
      pred_df <- pred_res$pred
    } else {
      return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    }

    if ("Input_1" %in% names(pred_df) && !"Input" %in% names(pred_df)) {
      pred_df <- pred_df %>% dplyr::rename(Input = Input_1)
    }

    # Filtrer Output 1
    pred_df_o1 <- pred_df %>%
      dplyr::filter(Output_ID == "1") %>%
      dplyr::mutate(Input = round(Input, 5)) %>%
      dplyr::distinct(Input, .keep_all = TRUE)

    mean_col <- intersect(names(pred_df_o1), c("Mean", "mean", "Prediction", "Mu"))
    var_col  <- intersect(names(pred_df_o1), c("Var", "var", "Variance", "Sigma2"))
    if (length(mean_col) == 0 || length(var_col) == 0) {
      return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    }

    # Dé-scaler les prédictions
    sf <- scale_factors %>%
      dplyr::filter(Output_ID == "1", Cluster_ID == cluster_id)
    if (nrow(sf) == 0) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    s <- sf$scale_factor[1]

    pred_df_o1[[mean_col[1]]] <- pred_df_o1[[mean_col[1]]] * s
    pred_df_o1[[var_col[1]]]  <- pred_df_o1[[var_col[1]]] * s^2

    # Dé-scaler la vérité aussi (elle est scalée)
    truth_o1_unscaled <- truth_o1 %>%
      dplyr::mutate(Output = Output * s, Input = round(Input, 5)) %>%
      dplyr::distinct(Input, .keep_all = TRUE)

    # Jointure
    merged <- truth_o1_unscaled %>%
      dplyr::inner_join(pred_df_o1, by = "Input", suffix = c("", "_pred"))

    if (nrow(merged) == 0) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))

    pred_mean <- merged[[mean_col[1]]]
    pred_var  <- merged[[var_col[1]]]
    truth     <- merged$Output

    rmse <- sqrt(mean((pred_mean - truth)^2))
    coverage_95 <- compute_coverage(pred_mean, pred_var, truth)

    # NLL mixture
    nll <- NA
    if ("pred_by_cluster" %in% names(pred_res)) {
      pred_by_clus <- pred_res$pred_by_cluster %>%
        dplyr::filter(Output_ID == "1")
      # Dé-scaler les prédictions par cluster
      if ("Mean" %in% names(pred_by_clus)) {
        pred_by_clus <- pred_by_clus %>%
          dplyr::mutate(Mean = Mean * s, Var = Var * s^2)
      }
      weights <- if ("weights" %in% names(pred_res)) pred_res$weights else pred_res$mixture
      if (is.data.frame(weights) && nrow(weights) == 1) {
        weights <- as.numeric(weights)
        names(weights) <- unique(pred_by_clus$Cluster)
      }
      nll <- compute_nll_mixture(pred_by_clus, weights, truth_o1_unscaled)
    } else {
      nll <- compute_nll_gaussian(pred_mean, pred_var, truth)
    }

    tibble(rmse = rmse, nll = nll, coverage_95 = coverage_95)
  }, error = function(e) {
    tibble(rmse = NA, nll = NA, coverage_95 = NA)
  })
}

# ---- Extraction des métriques pour MT ----
extract_metrics_mt <- function(pred_entry) {
  tryCatch({
    preds_by_output <- pred_entry$predictions_by_output
    truth_data      <- pred_entry$truth  # données brutes (non scalées)

    # On ne s'intéresse qu'à l'output 1
    if (!("1" %in% names(preds_by_output))) {
      return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    }

    pred_res <- preds_by_output[["1"]]$prediction
    truth_o1 <- truth_data %>% dplyr::filter(Output_ID == "1")
    if (nrow(truth_o1) == 0) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))

    # Extraire la prédiction agrégée
    if (is.data.frame(pred_res) || tibble::is_tibble(pred_res)) {
      pred_df <- pred_res
    } else if ("mixture_pred" %in% names(pred_res)) {
      pred_df <- pred_res$mixture_pred
    } else if ("mixture" %in% names(pred_res) && is.data.frame(pred_res$mixture)) {
      pred_df <- pred_res$mixture
    } else if ("pred" %in% names(pred_res)) {
      pred_df <- pred_res$pred
    } else {
      return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    }

    if ("Input_1" %in% names(pred_df) && !"Input" %in% names(pred_df)) {
      pred_df <- pred_df %>% dplyr::rename(Input = Input_1)
    }

    # MT : Output_ID = "1" (toujours)
    pred_df <- pred_df %>%
      dplyr::mutate(Input = round(Input, 5)) %>%
      dplyr::distinct(Input, .keep_all = TRUE)

    mean_col <- intersect(names(pred_df), c("Mean", "mean", "Prediction", "Mu"))
    var_col  <- intersect(names(pred_df), c("Var", "var", "Variance", "Sigma2"))
    if (length(mean_col) == 0 || length(var_col) == 0) {
      return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    }

    truth_o1_clean <- truth_o1 %>%
      dplyr::mutate(Input = round(Input, 5)) %>%
      dplyr::distinct(Input, .keep_all = TRUE)

    merged <- truth_o1_clean %>%
      dplyr::inner_join(pred_df, by = "Input", suffix = c("", "_pred"))

    if (nrow(merged) == 0) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))

    pred_mean <- merged[[mean_col[1]]]
    pred_var  <- merged[[var_col[1]]]
    truth     <- merged$Output

    rmse <- sqrt(mean((pred_mean - truth)^2))
    coverage_95 <- compute_coverage(pred_mean, pred_var, truth)

    # NLL mixture pour MT
    nll <- NA
    if ("pred_by_cluster" %in% names(pred_res)) {
      pred_by_clus <- pred_res$pred_by_cluster
      weights <- if ("weights" %in% names(pred_res)) pred_res$weights else pred_res$mixture
      if (is.data.frame(weights) && nrow(weights) == 1) {
        weights <- as.numeric(weights)
        names(weights) <- unique(pred_by_clus$Cluster)
      }
      nll <- compute_nll_mixture(pred_by_clus, weights, truth_o1_clean)
    } else {
      nll <- compute_nll_gaussian(pred_mean, pred_var, truth)
    }

    tibble(rmse = rmse, nll = nll, coverage_95 = coverage_95)
  }, error = function(e) {
    tibble(rmse = NA, nll = NA, coverage_95 = NA)
  })
}

# ---- Extraction des métriques pour MO ----
extract_metrics_mo <- function(pred_entry, scale_factors) {
  tryCatch({
    pred_df    <- pred_entry$prediction  # déjà filtré sur l'Output_ID MO de la tâche
    truth_data <- pred_entry$truth  # données scalées, Output_ID = "1"
    cluster_id <- pred_entry$cluster_id

    if (nrow(pred_df) == 0 || nrow(truth_data) == 0) {
      return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    }

    # Dé-scaler
    sf <- scale_factors %>%
      dplyr::filter(Output_ID == "1", Cluster_ID == cluster_id)
    if (nrow(sf) == 0) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    s <- sf$scale_factor[1]

    if ("Input_1" %in% names(pred_df) && !"Input" %in% names(pred_df)) {
      pred_df <- pred_df %>% dplyr::rename(Input = Input_1)
    }

    mean_col <- intersect(names(pred_df), c("Mean", "mean", "Prediction", "Mu"))
    var_col  <- intersect(names(pred_df), c("Var", "var", "Variance", "Sigma2"))
    if (length(mean_col) == 0 || length(var_col) == 0) {
      return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    }

    pred_df <- pred_df %>%
      dplyr::mutate(
        Input = round(Input, 5),
        !!mean_col[1] := .data[[mean_col[1]]] * s,
        !!var_col[1]  := .data[[var_col[1]]] * s^2
      ) %>%
      dplyr::distinct(Input, .keep_all = TRUE)

    truth_clean <- truth_data %>%
      dplyr::filter(Output_ID == "1") %>%
      dplyr::mutate(Output = Output * s, Input = round(Input, 5)) %>%
      dplyr::distinct(Input, .keep_all = TRUE)

    merged <- truth_clean %>%
      dplyr::inner_join(pred_df, by = "Input", suffix = c("", "_pred"))

    if (nrow(merged) == 0) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))

    pred_mean <- merged[[mean_col[1]]]
    pred_var  <- merged[[var_col[1]]]
    truth     <- merged$Output

    rmse <- sqrt(mean((pred_mean - truth)^2))
    coverage_95 <- compute_coverage(pred_mean, pred_var, truth)
    nll <- compute_nll_gaussian(pred_mean, pred_var, truth)

    tibble(rmse = rmse, nll = nll, coverage_95 = coverage_95)
  }, error = function(e) {
    tibble(rmse = NA, nll = NA, coverage_95 = NA)
  })
}

# ==========================================================================
# BOUCLE PRINCIPALE
# ==========================================================================

username   <- Sys.getenv("USER")
base_dir   <- file.path("/scratch", username, "NeurIPS_experiments")
output_dir <- file.path(base_dir, "Metrics")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Toutes les configurations
configs <- expand.grid(
  n_out   = c(2, 4, 8),
  n_train = c(15, 30, 100),
  n_pred  = c(1, 10, 100),
  problem = c("interpolation", "forecasting"),
  stringsAsFactors = FALSE
) %>%
  # Filtrer les 7 configs uniques
  dplyr::filter(
    # Défaut : (2, 30, 1)
    (n_out == 2 & n_train == 30 & n_pred == 1) |
    # Variations nb outputs
    (n_out == 4 & n_train == 30 & n_pred == 1) |
    (n_out == 8 & n_train == 30 & n_pred == 1) |
    # Variations nb train
    (n_out == 2 & n_train == 15 & n_pred == 1) |
    (n_out == 2 & n_train == 100 & n_pred == 1) |
    # Variations nb pred
    (n_out == 2 & n_train == 30 & n_pred == 10) |
    (n_out == 2 & n_train == 30 & n_pred == 100)
  )

models <- c("MOMT", "MO", "MT")

# Seeds par modèle et configuration :
#   MOMT : 50 seeds pour toutes les configs
#   MT   : 50 seeds pour toutes les configs
#   MO   : variable selon la config (8/30/1, 2/100/1 et 2/30/100 retirées;
#          4/30/1 à 10 seeds)
get_seeds <- function(model, n_out, n_train, n_pred) {
  if (model %in% c("MOMT", "MT")) {
    return(1:50)
  }
  # MO : configs retirées
  if (n_out == 8 & n_train == 30 & n_pred == 1) return(integer(0))
  if (n_out == 2 & n_train == 100 & n_pred == 1) return(integer(0))
  if (n_out == 2 & n_train == 30 & n_pred == 100) return(integer(0))
  # MO : config à 10 seeds
  if (n_out == 4 & n_train == 30 & n_pred == 1) return(1:10)
  # MO : les autres configs → 50 seeds
  return(1:50)
}

cat("=== Calcul des métriques NeurIPS ===\n")
cat(paste0("  Date      : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste0("  Configs   : ", nrow(configs), "\n"))
cat(paste0("  Modèles   : ", paste(models, collapse = ", "), "\n"))
cat(paste0("  Seeds     : variables (MOMT/MT=50, MO=0/10/50 selon config)\n\n"))

all_metrics <- list()

for (i in 1:nrow(configs)) {
  cfg <- configs[i, ]
  config_label <- paste0("out", cfg$n_out, "_train", cfg$n_train, "_pred", cfg$n_pred)
  config_path  <- file.path(base_dir, config_label, cfg$problem)

  for (model in models) {
    seeds <- get_seeds(model, cfg$n_out, cfg$n_train, cfg$n_pred)
    if (length(seeds) == 0) next
    for (seed in seeds) {
      pred_dir  <- file.path(config_path, paste0("Predictions_", model))
      pred_file <- file.path(pred_dir, paste0("predictions_seed_", seed, ".rds"))

      if (!file.exists(pred_file)) {
        cat(paste0("  [SKIP] ", config_label, "/", cfg$problem, "/", model,
                   " seed=", seed, "\n"))
        next
      }

      pred_data <- readRDS(pred_file)

      # Charger les scale_factors si nécessaire (MOMT et MO)
      scale_factors <- NULL
      if (model %in% c("MOMT", "MO")) {
        ds_file <- file.path(config_path, "Datasets",
                             paste0("datasets_seed_", seed, ".rds"))
        if (file.exists(ds_file)) {
          ds <- readRDS(ds_file)
          scale_factors <- ds$scale_factors
        }
      }

      # Temps
      t_training  <- if (!is.null(pred_data$t_training)) pred_data$t_training else 0
      t_pred      <- if (!is.null(pred_data$t_pred_total)) pred_data$t_pred_total else 0
      t_hyperpost <- if (!is.null(pred_data$t_hyperpost_global)) {
        pred_data$t_hyperpost_global
      } else 0

      # Ajustement selon le modèle
      if (model == "MO") {
        t_train_mo <- if (!is.null(pred_data$t_train_total)) pred_data$t_train_total else t_training
        t_pred <- t_pred + t_train_mo
        t_training <- 0
      } else if (model %in% c("MOMT", "MT")) {
        t_pred <- t_pred + t_hyperpost
      }

      # Extraire métriques par tâche
      predictions <- pred_data$predictions
      task_metrics_list <- list()

      for (tid in names(predictions)) {
        if (model == "MOMT") {
          metrics_task <- extract_metrics_momt(predictions[[tid]], scale_factors)
        } else if (model == "MT") {
          metrics_task <- extract_metrics_mt(predictions[[tid]])
        } else if (model == "MO") {
          metrics_task <- extract_metrics_mo(predictions[[tid]], scale_factors)
        }
        task_metrics_list[[tid]] <- metrics_task
      }

      task_metrics_df <- bind_rows(task_metrics_list, .id = "task_id")

      rmse_all     <- mean(task_metrics_df$rmse, na.rm = TRUE)
      nll_all      <- mean(task_metrics_df$nll, na.rm = TRUE)
      coverage_all <- mean(task_metrics_df$coverage_95, na.rm = TRUE)

      all_metrics[[length(all_metrics) + 1]] <- tibble(
        n_out       = cfg$n_out,
        n_train     = cfg$n_train,
        n_pred      = cfg$n_pred,
        problem     = cfg$problem,
        model       = model,
        seed        = seed,
        rmse        = rmse_all,
        nll         = nll_all,
        coverage_95 = coverage_all,
        t_training  = t_training,
        t_pred      = t_pred,
        t_hyperpost = t_hyperpost
      )

      cat(paste0("  [OK] ", config_label, "/", cfg$problem, "/", model,
                 " seed=", seed, " RMSE=", round(rmse_all, 4), "\n"))
    }
  }
}

metrics_raw <- bind_rows(all_metrics)

# --- Sauvegarde des métriques brutes ---
write_csv(metrics_raw, file.path(output_dir, "metrics_raw_NeurIPS.csv"))
saveRDS(metrics_raw, file.path(output_dir, "metrics_raw_NeurIPS.rds"))
cat(paste0("\nMétriques brutes : ", nrow(metrics_raw), " lignes\n"))

# --- Résumé statistique ---
metrics_summary <- metrics_raw %>%
  dplyr::group_by(n_out, n_train, n_pred, problem, model) %>%
  dplyr::summarise(
    rmse_mean   = mean(rmse, na.rm = TRUE),
    rmse_std    = sd(rmse, na.rm = TRUE),
    rmse_median = median(rmse, na.rm = TRUE),
    rmse_q1     = quantile(rmse, 0.25, na.rm = TRUE),
    rmse_q3     = quantile(rmse, 0.75, na.rm = TRUE),
    rmse_iqr    = IQR(rmse, na.rm = TRUE),

    nll_mean   = mean(nll, na.rm = TRUE),
    nll_std    = sd(nll, na.rm = TRUE),
    nll_median = median(nll, na.rm = TRUE),
    nll_q1     = quantile(nll, 0.25, na.rm = TRUE),
    nll_q3     = quantile(nll, 0.75, na.rm = TRUE),
    nll_iqr    = IQR(nll, na.rm = TRUE),

    cov95_mean   = mean(coverage_95, na.rm = TRUE),
    cov95_std    = sd(coverage_95, na.rm = TRUE),
    cov95_median = median(coverage_95, na.rm = TRUE),
    cov95_q1     = quantile(coverage_95, 0.25, na.rm = TRUE),
    cov95_q3     = quantile(coverage_95, 0.75, na.rm = TRUE),
    cov95_iqr    = IQR(coverage_95, na.rm = TRUE),

    t_train_mean   = mean(t_training, na.rm = TRUE),
    t_train_std    = sd(t_training, na.rm = TRUE),
    t_train_median = median(t_training, na.rm = TRUE),
    t_train_iqr    = IQR(t_training, na.rm = TRUE),

    t_pred_mean   = mean(t_pred, na.rm = TRUE),
    t_pred_std    = sd(t_pred, na.rm = TRUE),
    t_pred_median = median(t_pred, na.rm = TRUE),
    t_pred_iqr    = IQR(t_pred, na.rm = TRUE),

    n_seeds = n(),
    .groups = "drop"
  )

# Colonnes formatées pour affichage
metrics_display <- metrics_summary %>%
  dplyr::mutate(
    `RMSE (mean±std)`     = paste0(round(rmse_mean, 4), " ± ", round(rmse_std, 4)),
    `RMSE (median[IQR])`  = paste0(round(rmse_median, 4), " [", round(rmse_q1, 4), ", ", round(rmse_q3, 4), "]"),
    `NLL (mean±std)`      = paste0(round(nll_mean, 4), " ± ", round(nll_std, 4)),
    `NLL (median[IQR])`   = paste0(round(nll_median, 4), " [", round(nll_q1, 4), ", ", round(nll_q3, 4), "]"),
    `Cov95 (mean±std)`    = paste0(round(cov95_mean, 4), " ± ", round(cov95_std, 4)),
    `Cov95 (median[IQR])` = paste0(round(cov95_median, 4), " [", round(cov95_q1, 4), ", ", round(cov95_q3, 4), "]"),
    `T_train (mean±std)`  = paste0(round(t_train_mean, 2), " ± ", round(t_train_std, 2)),
    `T_pred (mean±std)`   = paste0(round(t_pred_mean, 2), " ± ", round(t_pred_std, 2))
  ) %>%
  dplyr::select(
    n_out, n_train, n_pred, problem, model,
    `RMSE (mean±std)`, `RMSE (median[IQR])`,
    `NLL (mean±std)`, `NLL (median[IQR])`,
    `Cov95 (mean±std)`, `Cov95 (median[IQR])`,
    `T_train (mean±std)`, `T_pred (mean±std)`
  )

metrics_summary <- metrics_summary %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))

write_excel_csv2(metrics_summary, file.path(output_dir, "metrics_summary_NeurIPS.csv"))
write_excel_csv2(metrics_display, file.path(output_dir, "metrics_display_NeurIPS.csv"))
saveRDS(metrics_summary, file.path(output_dir, "metrics_summary_NeurIPS.rds"))

# --- Affichage par paramètre d'intérêt ---
cat("\n--- Résumé par paramètre d'intérêt ---\n\n")

# Nb outputs (n_train=30, n_pred=1)
cat("=== Paramètre : nb_outputs ===\n")
sub <- metrics_display %>%
  dplyr::filter(n_train == 30, n_pred == 1) %>%
  dplyr::arrange(problem, n_out, model)
print(as.data.frame(sub), row.names = FALSE)
cat("\n")

# Nb train (n_out=2, n_pred=1)
cat("=== Paramètre : nb_tasks_train ===\n")
sub <- metrics_display %>%
  dplyr::filter(n_out == 2, n_pred == 1) %>%
  dplyr::arrange(problem, n_train, model)
print(as.data.frame(sub), row.names = FALSE)
cat("\n")

# Nb pred (n_out=2, n_train=30)
cat("=== Paramètre : nb_tasks_pred ===\n")
sub <- metrics_display %>%
  dplyr::filter(n_out == 2, n_train == 30) %>%
  dplyr::arrange(problem, n_pred, model)
print(as.data.frame(sub), row.names = FALSE)
cat("\n")

cat(paste0("\nFichiers sauvegardés dans : ", output_dir, "\n"))
cat(paste0("Terminé : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
