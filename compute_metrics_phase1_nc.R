# ==========================================================================
# compute_metrics_phase1_nc.R
# Phase 1 NC : Calcul des métriques à partir des prédictions sauvegardées
#
# Usage :
#   Rscript compute_metrics_phase1_nc.R
#
# Lit toutes les prédictions (MOMT, MO, MT) NC pour chaque (param, config, seed),
# calcule RMSE, NLL, coverage 95%, temps, et produit les tableaux récapitulatifs.
#
# Sortie :
#   - metrics_raw_phase1_nc.csv     : métriques brutes par (param, config, model, seed)
#   - metrics_summary_phase1_nc.csv : (mean±std) et (median±IQR) par (param, config, model)
# ==========================================================================

library(tidyverse)

# ---- Couverture de l'intervalle de confiance (IC) ----
compute_coverage <- function(pred_mean, pred_var, truth, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  pred_sd <- sqrt(pmax(pred_var, 1e-10))
  lower <- pred_mean - z * pred_sd
  upper <- pred_mean + z * pred_sd
  mean(truth >= lower & truth <= upper)
}

username  <- Sys.getenv("USER")
base_dir  <- file.path("/scratch", username, "Phase1_NC_experiments")
output_dir <- file.path(base_dir, "Metrics")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

params  <- c("probleme", "hp_clusters", "hp_tasks", "inputs", "n_train")
configs <- c("default", "variation")
models  <- c("MOMT", "MO", "MT")
seeds   <- 1:5

cat("=== Calcul des métriques Phase 1 NC ===\n")
cat(paste0("  Date : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n"))

# --- Fonction d'extraction des métriques par tâche ---
extract_task_metrics <- function(pred_entry) {
  tryCatch({
    pred_res   <- pred_entry$prediction
    truth_data <- pred_entry$truth

    # 1. Extraire le data frame de prédiction
    pred_df <- NULL
    if (is.data.frame(pred_res) || tibble::is_tibble(pred_res)) {
      pred_df <- pred_res
    } else if ("pred" %in% names(pred_res)) {
      pred_df <- pred_res$pred
    } else if ("pred_gp" %in% names(pred_res)) {
      pred_df <- pred_res$pred_gp$pred
    } else {
      for (nm in names(pred_res)) {
        if (is.data.frame(pred_res[[nm]])) { pred_df <- pred_res[[nm]]; break }
      }
    }

    if (is.null(pred_df)) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))

    if ("Input_1" %in% names(pred_df) && !"Input" %in% names(pred_df)) {
      pred_df <- pred_df %>% dplyr::rename(Input = Input_1)
    }

    mean_col <- intersect(names(pred_df), c("Mean", "mean", "Prediction", "prediction", "Mu"))
    var_col  <- intersect(names(pred_df), c("Var", "var", "Variance", "variance", "Sigma2"))

    if (length(mean_col) == 0 || length(var_col) == 0) {
      return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    }

    # 2. Jointure prédiction <-> vérité
    truth_df <- truth_data %>%
      dplyr::mutate(Input = round(Input, 5), Output_ID = as.character(Output_ID)) %>%
      dplyr::distinct(Input, Output_ID, .keep_all = TRUE)

    pred_df_clean <- pred_df %>%
      dplyr::mutate(Input = round(Input, 5)) %>%
      dplyr::distinct(Input, Output_ID, .keep_all = TRUE)

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

    # 3. Calcul RMSE et Coverage 95%
    rmse <- sqrt(mean((pred_mean - truth)^2))
    coverage_95 <- compute_coverage(pred_mean, pred_var, truth)

    # 4. NLL : toujours Gaussien simple (pas de mixture sans clustering)
    pred_var_safe <- pmax(pred_var, 1e-12)
    nll <- mean(0.5 * log(2 * pi * pred_var_safe) + 0.5 * (truth - pred_mean)^2 / pred_var_safe)

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
        t_training  <- if (!is.null(pred_data$t_training)) pred_data$t_training else NA
        t_pred      <- if (!is.null(pred_data$t_pred_total)) pred_data$t_pred_total else NA
        t_hyperpost <- if (!is.null(pred_data$t_hyperpost_global)) pred_data$t_hyperpost_global else NA

        # Ajustement logique selon le modèle
        if (model == "MO") {
          t_train_mo <- if (!is.null(pred_data$t_train_total)) pred_data$t_train_total else 0
          t_pred <- t_pred + t_train_mo
          t_training <- 0
        } else if (model %in% c("MOMT", "MT")) {
          t_pred <- t_pred + t_hyperpost
        }

        # Extraire métriques par tâche
        predictions <- pred_data$predictions
        task_metrics_list <- list()

        for (tid in names(predictions)) {
          metrics_task <- extract_task_metrics(predictions[[tid]])
          task_metrics_list[[tid]] <- metrics_task
        }

        task_metrics_df <- bind_rows(task_metrics_list, .id = "task_id")

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
write_csv(metrics_raw, file.path(output_dir, "metrics_raw_phase1_nc.csv"))
saveRDS(metrics_raw, file.path(output_dir, "metrics_raw_phase1_nc.rds"))
cat(paste0("\nMétriques brutes sauvegardées : ", nrow(metrics_raw), " lignes\n"))

# --- Calcul des statistiques récapitulatives ---
metrics_summary <- metrics_raw %>%
  dplyr::group_by(param, config, model) %>%
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

metrics_display <- metrics_summary %>%
  dplyr::mutate(
    `RMSE (mean±std)`     = paste0(round(rmse_mean, 2), " ± ", round(rmse_std, 2)),
    `RMSE (median[IQR])`  = paste0(round(rmse_median, 2), " [", round(rmse_q1, 2), ", ", round(rmse_q3, 2), "]"),
    `NLL (mean±std)`      = paste0(round(nll_mean, 2), " ± ", round(nll_std, 2)),
    `NLL (median[IQR])`   = paste0(round(nll_median, 2), " [", round(nll_q1, 2), ", ", round(nll_q3, 2), "]"),
    `Cov95 (mean±std)`    = paste0(round(cov95_mean, 2), " ± ", round(cov95_std, 2)),
    `Cov95 (median[IQR])` = paste0(round(cov95_median, 2), " [", round(cov95_q1, 2), ", ", round(cov95_q3, 2), "]"),
    `T_train (mean±std)`  = paste0(round(t_train_mean, 2), " ± ", round(t_train_std, 2)),
    `T_train (med[IQR])`  = paste0(round(t_train_median, 2), " [", round(t_train_iqr, 2), "]"),
    `T_pred (mean±std)`   = paste0(round(t_pred_mean, 2), " ± ", round(t_pred_std, 2)),
    `T_pred (med[IQR])`   = paste0(round(t_pred_median, 2), " [", round(t_pred_iqr, 2), "]")
  ) %>%
  dplyr::select(
    param, config, model,
    `RMSE (mean±std)`, `RMSE (median[IQR])`,
    `NLL (mean±std)`, `NLL (median[IQR])`,
    `Cov95 (mean±std)`, `Cov95 (median[IQR])`,
    `T_train (mean±std)`, `T_pred (mean±std)`
  )

metrics_summary <- metrics_summary %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2)))

write_excel_csv2(metrics_summary, file.path(output_dir, "metrics_summary_phase1_nc.csv"))
write_excel_csv2(metrics_display, file.path(output_dir, "metrics_display_phase1_nc.csv"))
saveRDS(metrics_summary, file.path(output_dir, "metrics_summary_phase1_nc.rds"))

cat("\n--- Résumé par paramètre ---\n\n")

for (p in params) {
  cat(paste0("=== ", p, " ===\n"))
  sub <- metrics_display %>% dplyr::filter(param == p)
  print(sub, n = Inf, width = Inf)
  cat("\n")
}

cat(paste0("\n=== Métriques NC terminées : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ===\n"))
