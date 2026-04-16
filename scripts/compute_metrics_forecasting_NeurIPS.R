# ==========================================================================
# compute_metrics_forecasting_NeurIPS.R
# NeurIPS : Calcul des metriques pour le forecasting uniquement
#
# Usage :
#   Rscript compute_metrics_forecasting_NeurIPS.R
#
# Lit les predictions MOMT et MT dans l'arborescence :
#   /scratch/$USER/NeurIPS_experiments/{config_label}/forecasting/
#
# Sorties :
#   - metrics_raw_forecasting_NeurIPS.csv
#   - metrics_summary_forecasting_NeurIPS.csv
#   - metrics_display_forecasting_NeurIPS.csv
# ==========================================================================

library(tidyverse)
library(matrixStats)

compute_coverage <- function(pred_mean, pred_var, truth, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  pred_sd <- sqrt(pmax(pred_var, 1e-10))
  lower <- pred_mean - z * pred_sd
  upper <- pred_mean + z * pred_sd
  mean(truth >= lower & truth <= upper)
}

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
  k <- length(clusters)
  if (is.null(weights)) {
    weights <- rep(1 / k, k)
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

  mean(nll_per_point$nll)
}

compute_nll_gaussian <- function(pred_mean, pred_var, truth) {
  pred_var_safe <- pmax(pred_var, 1e-12)
  mean(0.5 * log(2 * pi * pred_var_safe) + 0.5 * (truth - pred_mean)^2 / pred_var_safe)
}

extract_metrics_momt <- function(pred_entry, scale_factors) {
  tryCatch({
    pred_res <- pred_entry$prediction
    truth_data <- pred_entry$truth
    cluster_id <- pred_entry$cluster_id

    truth_o1 <- truth_data %>% dplyr::filter(Output_ID == "1")
    if (nrow(truth_o1) == 0) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))

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

    pred_df_o1 <- pred_df %>%
      dplyr::filter(Output_ID == "1") %>%
      dplyr::mutate(Input = round(Input, 5)) %>%
      dplyr::distinct(Input, .keep_all = TRUE)

    mean_col <- intersect(names(pred_df_o1), c("Mean", "mean", "Prediction", "Mu"))
    var_col <- intersect(names(pred_df_o1), c("Var", "var", "Variance", "Sigma2"))
    if (length(mean_col) == 0 || length(var_col) == 0) {
      return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    }

    sf <- scale_factors %>%
      dplyr::filter(Output_ID == "1", Cluster_ID == cluster_id)
    if (nrow(sf) == 0) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    s <- sf$scale_factor[1]

    pred_df_o1[[mean_col[1]]] <- pred_df_o1[[mean_col[1]]] * s
    pred_df_o1[[var_col[1]]] <- pred_df_o1[[var_col[1]]] * s^2

    truth_o1_unscaled <- truth_o1 %>%
      dplyr::mutate(Output = Output * s, Input = round(Input, 5)) %>%
      dplyr::distinct(Input, .keep_all = TRUE)

    merged <- truth_o1_unscaled %>%
      dplyr::inner_join(pred_df_o1, by = "Input", suffix = c("", "_pred"))

    if (nrow(merged) == 0) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))

    pred_mean <- merged[[mean_col[1]]]
    pred_var <- merged[[var_col[1]]]
    truth <- merged$Output

    rmse <- sqrt(mean((pred_mean - truth)^2))
    coverage_95 <- compute_coverage(pred_mean, pred_var, truth)

    if ("pred_by_cluster" %in% names(pred_res)) {
      pred_by_clus <- pred_res$pred_by_cluster %>% dplyr::filter(Output_ID == "1")
      if ("Mean" %in% names(pred_by_clus)) {
        pred_by_clus <- pred_by_clus %>% dplyr::mutate(Mean = Mean * s, Var = Var * s^2)
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

extract_metrics_mt <- function(pred_entry) {
  tryCatch({
    preds_by_output <- pred_entry$predictions_by_output
    truth_data <- pred_entry$truth

    if (!("1" %in% names(preds_by_output))) {
      return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    }

    pred_res <- preds_by_output[["1"]]$prediction
    truth_o1 <- truth_data %>% dplyr::filter(Output_ID == "1")
    if (nrow(truth_o1) == 0) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))

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

    pred_df <- pred_df %>%
      dplyr::mutate(Input = round(Input, 5)) %>%
      dplyr::distinct(Input, .keep_all = TRUE)

    mean_col <- intersect(names(pred_df), c("Mean", "mean", "Prediction", "Mu"))
    var_col <- intersect(names(pred_df), c("Var", "var", "Variance", "Sigma2"))
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
    pred_var <- merged[[var_col[1]]]
    truth <- merged$Output

    rmse <- sqrt(mean((pred_mean - truth)^2))
    coverage_95 <- compute_coverage(pred_mean, pred_var, truth)

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

username <- Sys.getenv("USER")
base_dir <- file.path("/scratch", username, "NeurIPS_experiments")
output_dir <- file.path(base_dir, "Metrics")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

base_configs <- tibble::tribble(
  ~n_out, ~n_train, ~n_pred, ~n_clust,
  4, 30, 1, 1,
  8, 30, 1, 1,
  2, 30, 1, 1,
  2, 15, 1, 1,
  2, 100, 1, 1,
  2, 30, 10, 1,
  2, 30, 100, 1,
  2, 30, 1, 2,
  2, 30, 1, 3,
  2, 30, 1, 4
)

models <- c("MOMT", "MT")
problem <- "forecasting"

get_config_label <- function(n_out, n_train, n_pred, n_clust) {
  paste0("out", n_out, "_train", n_train, "_pred", n_pred, "_clust", n_clust)
}

get_prediction_dir <- function(config_label, model) {
  file.path(base_dir, config_label, problem, paste0("Predictions_", model))
}

get_dataset_file <- function(config_label, seed) {
  file.path(base_dir, config_label, problem, "Datasets", paste0("datasets_seed_", seed, ".rds"))
}

get_available_seeds <- function(pred_dir) {
  if (!dir.exists(pred_dir)) return(integer(0))

  pred_files <- list.files(
    pred_dir,
    pattern = "^predictions_seed_[0-9]+\\.rds$",
    full.names = FALSE
  )
  if (length(pred_files) == 0) return(integer(0))

  pred_files %>%
    stringr::str_match("^predictions_seed_([0-9]+)\\.rds$") %>%
    .[, 2] %>%
    as.integer() %>%
    sort()
}

cat("=== Calcul des metriques NeurIPS (forecasting uniquement) ===\n")
cat(paste0("  Date      : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste0("  Configs   : ", nrow(base_configs), "\n"))
cat(paste0("  Modeles   : ", paste(models, collapse = ", "), "\n"))
cat("  Seeds     : detectees automatiquement depuis les fichiers de prediction\n")
cat(paste0("  Resultats : ", base_dir, "/{config}/forecasting\n\n"))

all_metrics <- list()

for (i in 1:nrow(base_configs)) {
  cfg <- base_configs[i, ]
  config_label <- get_config_label(cfg$n_out, cfg$n_train, cfg$n_pred, cfg$n_clust)

  for (model in models) {
    pred_dir <- get_prediction_dir(config_label, model)
    seeds <- get_available_seeds(pred_dir)
    if (length(seeds) == 0) next

    for (seed in seeds) {
      pred_file <- file.path(pred_dir, paste0("predictions_seed_", seed, ".rds"))
      if (!file.exists(pred_file)) next

      pred_data <- readRDS(pred_file)

      scale_factors <- NULL
      if (model == "MOMT") {
        ds_file <- get_dataset_file(config_label, seed)
        if (file.exists(ds_file)) {
          ds <- readRDS(ds_file)
          scale_factors <- ds$scale_factors
        }
      }

      t_training <- if (!is.null(pred_data$t_training)) pred_data$t_training else 0
      t_pred <- if (!is.null(pred_data$t_pred_total)) pred_data$t_pred_total else 0
      t_hyperpost <- if (!is.null(pred_data$t_hyperpost_global)) pred_data$t_hyperpost_global else 0
      t_pred <- t_pred + t_hyperpost

      predictions <- pred_data$predictions
      task_metrics_list <- list()

      for (tid in names(predictions)) {
        if (model == "MOMT") {
          metrics_task <- extract_metrics_momt(predictions[[tid]], scale_factors)
        } else {
          metrics_task <- extract_metrics_mt(predictions[[tid]])
        }
        task_metrics_list[[tid]] <- metrics_task
      }

      task_metrics_df <- bind_rows(task_metrics_list, .id = "task_id")

      rmse_all <- mean(task_metrics_df$rmse, na.rm = TRUE)
      nll_all <- mean(task_metrics_df$nll, na.rm = TRUE)
      coverage_all <- mean(task_metrics_df$coverage_95, na.rm = TRUE)

      all_metrics[[length(all_metrics) + 1]] <- tibble(
        n_out = cfg$n_out,
        n_train = cfg$n_train,
        n_pred = cfg$n_pred,
        n_clust = cfg$n_clust,
        problem = problem,
        model = model,
        seed = seed,
        rmse = rmse_all,
        nll = nll_all,
        coverage_95 = coverage_all,
        t_training = t_training,
        t_pred = t_pred,
        t_hyperpost = t_hyperpost
      )

      cat(paste0("  [OK] ", config_label, "/", problem, "/", model,
                 " seed=", seed, " RMSE=", round(rmse_all, 4), "\n"))
    }
  }
}

metrics_raw <- bind_rows(all_metrics)

if (nrow(metrics_raw) == 0) {
  stop("Aucune prediction forecasting exploitable trouvee dans les repertoires configures.")
}

write_csv(metrics_raw, file.path(output_dir, "metrics_raw_forecasting_NeurIPS.csv"))
saveRDS(metrics_raw, file.path(output_dir, "metrics_raw_forecasting_NeurIPS.rds"))

metrics_summary <- metrics_raw %>%
  dplyr::group_by(n_out, n_train, n_pred, n_clust, problem, model) %>%
  dplyr::summarise(
    rmse_mean = mean(rmse, na.rm = TRUE),
    rmse_std = sd(rmse, na.rm = TRUE),
    rmse_median = median(rmse, na.rm = TRUE),
    rmse_q1 = quantile(rmse, 0.25, na.rm = TRUE),
    rmse_q3 = quantile(rmse, 0.75, na.rm = TRUE),
    rmse_iqr = IQR(rmse, na.rm = TRUE),
    nll_mean = mean(nll, na.rm = TRUE),
    nll_std = sd(nll, na.rm = TRUE),
    nll_median = median(nll, na.rm = TRUE),
    nll_q1 = quantile(nll, 0.25, na.rm = TRUE),
    nll_q3 = quantile(nll, 0.75, na.rm = TRUE),
    nll_iqr = IQR(nll, na.rm = TRUE),
    cov95_mean = mean(coverage_95, na.rm = TRUE),
    cov95_std = sd(coverage_95, na.rm = TRUE),
    cov95_median = median(coverage_95, na.rm = TRUE),
    cov95_q1 = quantile(coverage_95, 0.25, na.rm = TRUE),
    cov95_q3 = quantile(coverage_95, 0.75, na.rm = TRUE),
    cov95_iqr = IQR(coverage_95, na.rm = TRUE),
    t_train_mean = mean(t_training, na.rm = TRUE),
    t_train_std = sd(t_training, na.rm = TRUE),
    t_train_median = median(t_training, na.rm = TRUE),
    t_train_iqr = IQR(t_training, na.rm = TRUE),
    t_pred_mean = mean(t_pred, na.rm = TRUE),
    t_pred_std = sd(t_pred, na.rm = TRUE),
    t_pred_median = median(t_pred, na.rm = TRUE),
    t_pred_iqr = IQR(t_pred, na.rm = TRUE),
    n_seeds = n(),
    .groups = "drop"
  )

metrics_display <- metrics_summary %>%
  dplyr::mutate(
    `RMSE (mean+-std)` = paste0(round(rmse_mean, 4), " +- ", round(rmse_std, 4)),
    `RMSE (median[IQR])` = paste0(round(rmse_median, 4), " [", round(rmse_q1, 4), ", ", round(rmse_q3, 4), "]"),
    `NLL (mean+-std)` = paste0(round(nll_mean, 4), " +- ", round(nll_std, 4)),
    `NLL (median[IQR])` = paste0(round(nll_median, 4), " [", round(nll_q1, 4), ", ", round(nll_q3, 4), "]"),
    `Cov95 (mean+-std)` = paste0(round(cov95_mean, 4), " +- ", round(cov95_std, 4)),
    `Cov95 (median[IQR])` = paste0(round(cov95_median, 4), " [", round(cov95_q1, 4), ", ", round(cov95_q3, 4), "]"),
    `T_train (mean+-std)` = paste0(round(t_train_mean, 2), " +- ", round(t_train_std, 2)),
    `T_pred (mean+-std)` = paste0(round(t_pred_mean, 2), " +- ", round(t_pred_std, 2))
  ) %>%
  dplyr::select(
    n_out, n_train, n_pred, n_clust, problem, model,
    `RMSE (mean+-std)`, `RMSE (median[IQR])`,
    `NLL (mean+-std)`, `NLL (median[IQR])`,
    `Cov95 (mean+-std)`, `Cov95 (median[IQR])`,
    `T_train (mean+-std)`, `T_pred (mean+-std)`
  )

metrics_summary <- metrics_summary %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))

write_excel_csv2(metrics_summary, file.path(output_dir, "metrics_summary_forecasting_NeurIPS.csv"))
write_excel_csv2(metrics_display, file.path(output_dir, "metrics_display_forecasting_NeurIPS.csv"))
saveRDS(metrics_summary, file.path(output_dir, "metrics_summary_forecasting_NeurIPS.rds"))

cat("\n--- Resume forecasting par parametre d'interet ---\n\n")

cat("=== Parametre : nb_outputs ===\n")
sub <- metrics_display %>%
  dplyr::filter(n_train == 30, n_pred == 1, n_clust == 1) %>%
  dplyr::arrange(n_out, model)
print(as.data.frame(sub), row.names = FALSE)
cat("\n")

cat("=== Parametre : nb_tasks_train ===\n")
sub <- metrics_display %>%
  dplyr::filter(n_out == 2, n_pred == 1, n_clust == 1) %>%
  dplyr::arrange(n_train, model)
print(as.data.frame(sub), row.names = FALSE)
cat("\n")

cat("=== Parametre : nb_tasks_pred ===\n")
sub <- metrics_display %>%
  dplyr::filter(n_out == 2, n_train == 30, n_clust == 1) %>%
  dplyr::arrange(n_pred, model)
print(as.data.frame(sub), row.names = FALSE)
cat("\n")

cat("=== Parametre : nb_clusters ===\n")
sub <- metrics_display %>%
  dplyr::filter(n_out == 2, n_train == 30, n_pred == 1) %>%
  dplyr::arrange(n_clust, model)
print(as.data.frame(sub), row.names = FALSE)
cat("\n")

cat(paste0("\nFichiers sauvegardes dans : ", output_dir, "\n"))
cat(paste0("Termine : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))