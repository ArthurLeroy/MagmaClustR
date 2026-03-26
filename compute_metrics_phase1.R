# ==========================================================================
# compute_metrics_phase1.R
# Phase 1 : Calcul des métriques à partir des prédictions sauvegardées
#
# Usage :
#   Rscript compute_metrics_phase1.R
#
# Lit toutes les prédictions (MOMT, MO, MT) pour chaque (param, config, seed),
# calcule RMSE, NLL, coverage 95%, temps, et produit les tableaux récapitulatifs.
#
# Sortie :
#   - metrics_raw_phase1.csv     : métriques brutes par (param, config, model, seed)
#   - metrics_summary_phase1.csv : (mean±std) et (median±IQR) par (param, config, model)
# ==========================================================================

library(tidyverse)

username  <- Sys.getenv("USER")
base_dir  <- file.path("/scratch", username, "Phase1_experiments")
output_dir <- file.path(base_dir, "Metrics")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

params  <- c("probleme", "hp_clusters", "hp_tasks", "inputs", "n_out", "n_train")
configs <- c("default", "variation")
models  <- c("MOMT", "MO", "MT")
seeds   <- 1:5

cat("=== Calcul des métriques Phase 1 ===\n")
cat(paste0("  Date : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n"))

# --- Fonction d'extraction des métriques par tâche ---
extract_task_metrics <- function(pred_entry) {
  tryCatch({
    pred_res   <- pred_entry$prediction
    truth_data <- pred_entry$truth

    # Extraire le data frame de prédiction
    if (is.data.frame(pred_res) || tibble::is_tibble(pred_res)) {
      pred_df <- pred_res
    } else if ("pred_mixture" %in% names(pred_res)) {
      pred_df <- pred_res$pred_mixture
    } else if ("mixture" %in% names(pred_res)) {
      pred_df <- pred_res$mixture
    } else if ("pred" %in% names(pred_res)) {
      pred_df <- pred_res$pred
    } else {
      # Dernier recours : prendre le premier élément data.frame
      for (nm in names(pred_res)) {
        if (is.data.frame(pred_res[[nm]])) { pred_df <- pred_res[[nm]]; break }
      }
    }

    if (!exists("pred_df")) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))

    # Identifier les colonnes de mean et variance
    mean_col <- intersect(names(pred_df), c("Mean", "mean", "Prediction", "prediction", "Mu"))
    var_col  <- intersect(names(pred_df), c("Var", "var", "Variance", "variance", "Sigma2"))

    if (length(mean_col) == 0 || length(var_col) == 0) {
      return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    }

    # Jointure prédiction <-> vérité
    truth_df <- truth_data %>%
      dplyr::mutate(Input = round(Input, 8), Output_ID = as.character(Output_ID))

    pred_df_clean <- pred_df %>%
      dplyr::mutate(Input = round(Input, 8))

    if ("Output_ID" %in% names(pred_df_clean)) {
      pred_df_clean <- pred_df_clean %>% dplyr::mutate(Output_ID = as.character(Output_ID))
      merged <- truth_df %>%
        dplyr::inner_join(pred_df_clean, by = c("Input", "Output_ID"), suffix = c("", "_pred"))
    } else {
      merged <- truth_df %>%
        dplyr::inner_join(pred_df_clean, by = "Input", suffix = c("", "_pred"))
    }

    if (nrow(merged) == 0) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))

    pred_mean <- merged[[mean_col[1]]]
    pred_var  <- merged[[var_col[1]]]
    truth     <- merged$Output

    # Protéger contre les variances nulles ou négatives
    pred_var <- pmax(pred_var, 1e-12)

    # RMSE
    rmse <- sqrt(mean((pred_mean - truth)^2))

    # NLL (par point, distribution gaussienne)
    nll <- mean(0.5 * log(2 * pi * pred_var) + 0.5 * (truth - pred_mean)^2 / pred_var)

    # Coverage 95%
    lower <- pred_mean - 1.96 * sqrt(pred_var)
    upper <- pred_mean + 1.96 * sqrt(pred_var)
    coverage_95 <- mean(truth >= lower & truth <= upper)

    tibble(rmse = rmse, nll = nll, coverage_95 = coverage_95)
  }, error = function(e) {
    tibble(rmse = NA, nll = NA, coverage_95 = NA)
  })
}

# --- Boucle principale ---
all_metrics <- list()

for (param in params) {
  for (config in configs) {
    for (model in models) {
      for (seed in seeds) {
        pred_dir  <- file.path(base_dir, param, config, paste0("Predictions_", model))
        pred_file <- file.path(pred_dir, paste0("predictions_seed_", seed, ".rds"))

        if (!file.exists(pred_file)) {
          cat(paste0("  [SKIP] ", param, "/", config, "/", model, " seed=", seed, " (fichier manquant)\n"))
          next
        }

        pred_data <- readRDS(pred_file)

        # Extraire temps
        t_training <- if (!is.null(pred_data$t_training)) pred_data$t_training else NA
        t_pred     <- if (!is.null(pred_data$t_pred_total)) pred_data$t_pred_total else NA
        t_hyperpost <- if (!is.null(pred_data$t_hyperpost_global)) pred_data$t_hyperpost_global else NA

        # Pour MO, le temps d'entraînement est t_train_total
        if (model == "MO" && is.na(t_training) && !is.null(pred_data$t_train_total)) {
          t_training <- pred_data$t_train_total
        }

        # Extraire métriques par tâche
        predictions <- pred_data$predictions
        task_metrics_list <- list()

        for (tid in names(predictions)) {
          metrics_task <- extract_task_metrics(predictions[[tid]])
          task_metrics_list[[tid]] <- metrics_task
        }

        task_metrics_df <- bind_rows(task_metrics_list, .id = "task_id")

        # Moyennes sur les tâches pour cette seed
        rmse_mean_tasks     <- mean(task_metrics_df$rmse, na.rm = TRUE)
        nll_mean_tasks      <- mean(task_metrics_df$nll, na.rm = TRUE)
        coverage_mean_tasks <- mean(task_metrics_df$coverage_95, na.rm = TRUE)

        all_metrics[[length(all_metrics) + 1]] <- tibble(
          param       = param,
          config      = config,
          model       = model,
          seed        = seed,
          rmse        = rmse_mean_tasks,
          nll         = nll_mean_tasks,
          coverage_95 = coverage_mean_tasks,
          t_training  = t_training,
          t_pred      = t_pred,
          t_hyperpost = t_hyperpost
        )
      }
    }
  }
}

metrics_raw <- bind_rows(all_metrics)

# --- Sauvegarde des métriques brutes ---
write_csv(metrics_raw, file.path(output_dir, "metrics_raw_phase1.csv"))
saveRDS(metrics_raw, file.path(output_dir, "metrics_raw_phase1.rds"))
cat(paste0("\nMétriques brutes sauvegardées : ", nrow(metrics_raw), " lignes\n"))

# --- Calcul des statistiques récapitulatives ---
metrics_summary <- metrics_raw %>%
  dplyr::group_by(param, config, model) %>%
  dplyr::summarise(
    # RMSE
    rmse_mean   = mean(rmse, na.rm = TRUE),
    rmse_std    = sd(rmse, na.rm = TRUE),
    rmse_median = median(rmse, na.rm = TRUE),
    rmse_q1     = quantile(rmse, 0.25, na.rm = TRUE),
    rmse_q3     = quantile(rmse, 0.75, na.rm = TRUE),
    rmse_iqr    = IQR(rmse, na.rm = TRUE),

    # NLL
    nll_mean   = mean(nll, na.rm = TRUE),
    nll_std    = sd(nll, na.rm = TRUE),
    nll_median = median(nll, na.rm = TRUE),
    nll_q1     = quantile(nll, 0.25, na.rm = TRUE),
    nll_q3     = quantile(nll, 0.75, na.rm = TRUE),
    nll_iqr    = IQR(nll, na.rm = TRUE),

    # Coverage 95%
    cov95_mean   = mean(coverage_95, na.rm = TRUE),
    cov95_std    = sd(coverage_95, na.rm = TRUE),
    cov95_median = median(coverage_95, na.rm = TRUE),
    cov95_q1     = quantile(coverage_95, 0.25, na.rm = TRUE),
    cov95_q3     = quantile(coverage_95, 0.75, na.rm = TRUE),
    cov95_iqr    = IQR(coverage_95, na.rm = TRUE),

    # Temps d'entraînement (secondes)
    t_train_mean   = mean(t_training, na.rm = TRUE),
    t_train_std    = sd(t_training, na.rm = TRUE),
    t_train_median = median(t_training, na.rm = TRUE),
    t_train_iqr    = IQR(t_training, na.rm = TRUE),

    # Temps de prédiction (secondes)
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
    `T_train (mean±std)`  = paste0(round(t_train_mean, 1), " ± ", round(t_train_std, 1)),
    `T_train (med[IQR])`  = paste0(round(t_train_median, 1), " [", round(t_train_iqr, 1), "]"),
    `T_pred (mean±std)`   = paste0(round(t_pred_mean, 1), " ± ", round(t_pred_std, 1)),
    `T_pred (med[IQR])`   = paste0(round(t_pred_median, 1), " [", round(t_pred_iqr, 1), "]")
  )

# Sauvegarde
write_csv(metrics_summary, file.path(output_dir, "metrics_summary_phase1.csv"))
write_csv(metrics_display, file.path(output_dir, "metrics_display_phase1.csv"))
saveRDS(metrics_summary, file.path(output_dir, "metrics_summary_phase1.rds"))

cat("\n--- Résumé par paramètre ---\n\n")

for (p in params) {
  cat(paste0("=== Paramètre : ", p, " ===\n"))
  sub <- metrics_display %>%
    dplyr::filter(param == p) %>%
    dplyr::select(config, model,
                  `RMSE (mean±std)`, `RMSE (median[IQR])`,
                  `NLL (mean±std)`, `NLL (median[IQR])`,
                  `Cov95 (mean±std)`, `Cov95 (median[IQR])`,
                  `T_train (mean±std)`, `T_pred (mean±std)`)
  print(as.data.frame(sub), row.names = FALSE)
  cat("\n")
}

cat(paste0("\nFichiers sauvegardés dans : ", output_dir, "\n"))
cat("  - metrics_raw_phase1.csv / .rds\n")
cat("  - metrics_summary_phase1.csv / .rds\n")
cat("  - metrics_display_phase1.csv\n")
cat(paste0("\nTerminé : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
