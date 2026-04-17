# ==========================================================================
# compute_metrics_forecasting_NeurIPS.R
# NeurIPS : Calcul des metriques pour le forecasting uniquement (MOMT + MT)
#
# Usage :
#   Rscript compute_metrics_forecasting_NeurIPS.R
#
# Basé sur compute_orginal.R (qui extrait correctement MOMT et MT).
# Ajout du paramètre n_clust et gestion du format pred_magma() (n_clust=1).
#
# IMPORTANT : Les métriques sont calculées sur OUTPUT 1 uniquement,
# en espace original (dé-scalé pour MOMT).
#
# Sorties :
#   - metrics_raw_forecasting_NeurIPS.csv
#   - metrics_summary_forecasting_NeurIPS.csv
#   - metrics_display_forecasting_NeurIPS.csv
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

# ---- NLL multivariée gaussienne (matrice de covariance complète) ----
compute_nll_multivariate_gaussian <- function(pred_mean, cov_matrix, truth) {
  n <- length(truth)
  if (n == 0 || is.null(cov_matrix)) return(NA)
  if (nrow(cov_matrix) != n || ncol(cov_matrix) != n) return(NA)
  diff <- truth - pred_mean
  cov_safe <- cov_matrix + diag(1e-8, n)
  L <- tryCatch(chol(cov_safe), error = function(e) NULL)
  if (is.null(L)) return(NA)
  log_det <- 2 * sum(log(diag(L)))
  alpha <- backsolve(L, diff, transpose = TRUE)
  quad_form <- sum(alpha^2)
  0.5 * (n * log(2 * pi) + log_det + quad_form)
}

# ---- NLL multivariée pour mixture de gaussiennes ----
compute_nll_multivariate_mixture <- function(cluster_preds, weights, truth) {
  K <- length(cluster_preds)
  n <- length(truth)
  if (K == 0 || n == 0) return(NA)
  clusters <- names(cluster_preds)
  if (is.null(names(weights))) names(weights) <- clusters
  w <- weights[clusters]
  w <- w / sum(w)
  log_densities <- numeric(K)
  for (i in seq_along(clusters)) {
    k <- clusters[i]
    mu <- cluster_preds[[k]]$mean
    Sigma <- cluster_preds[[k]]$cov
    if (is.null(Sigma) || length(mu) != n) return(NA)
    if (nrow(Sigma) != n || ncol(Sigma) != n) return(NA)
    Sigma_safe <- Sigma + diag(1e-8, n)
    L <- tryCatch(chol(Sigma_safe), error = function(e) NULL)
    if (is.null(L)) return(NA)
    log_det <- 2 * sum(log(diag(L)))
    diff <- truth - mu
    alpha <- backsolve(L, diff, transpose = TRUE)
    quad_form <- sum(alpha^2)
    log_densities[i] <- -0.5 * (n * log(2 * pi) + log_det + quad_form)
  }
  log_mixture <- matrixStats::logSumExp(log(w) + log_densities)
  return(-log_mixture)
}

# ---- Extraction de la matrice de covariance ----
extract_pred_cov <- function(pred_res) {
  if ("pred_gp" %in% names(pred_res)) {
    pg <- pred_res$pred_gp
    if (is.list(pg) && !is.data.frame(pg) && "cov" %in% names(pg)) {
      return(pg$cov)
    }
  }
  return(NULL)
}

# ---- Extraction du dataframe de prédiction ----
# Gère les formats de sortie de pred_magmaclust() ET pred_magma()
#   pred_magmaclust() → list(mixture_pred = tibble, pred = ..., mixture = ...)
#   pred_magma()      → list(pred_gp = list(pred = tibble, cov = ...), ponderation_matrix = ...)
#   pred_magma() sans get_full_cov → list(pred_gp = tibble, ponderation_matrix = ...)
extract_pred_df <- function(pred_res) {
  # Cas 1 : dataframe direct
  if (is.data.frame(pred_res) || tibble::is_tibble(pred_res)) {
    return(pred_res)
  }
  # Cas 2 : pred_magmaclust() → $mixture_pred
  if ("mixture_pred" %in% names(pred_res)) {
    return(pred_res$mixture_pred)
  }
  # Cas 3 : pred_magma() → $pred_gp (list ou tibble)
  if ("pred_gp" %in% names(pred_res)) {
    pg <- pred_res$pred_gp
    if (is.data.frame(pg) || tibble::is_tibble(pg)) {
      return(pg)
    }
    if (is.list(pg) && "pred" %in% names(pg)) {
      return(pg$pred)
    }
    return(NULL)
  }
  # Cas 4 : $mixture comme dataframe
  if ("mixture" %in% names(pred_res) && is.data.frame(pred_res$mixture)) {
    return(pred_res$mixture)
  }
  # Cas 5 : $pred
  if ("pred" %in% names(pred_res)) {
    return(pred_res$pred)
  }
  return(NULL)
}

# ---- Extraction des métriques pour MOMT ----
extract_metrics_momt <- function(pred_entry, scale_factors) {
  tryCatch({
    pred_res   <- pred_entry$prediction
    truth_data <- pred_entry$truth
    cluster_id <- pred_entry$cluster_id

    # Filtrer Output 1 seulement pour les métriques
    truth_o1 <- truth_data %>% dplyr::filter(Output_ID == "1")
    if (nrow(truth_o1) == 0) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))

    # Extraire la prédiction agrégée (gère pred_magmaclust ET pred_magma)
    pred_df <- extract_pred_df(pred_res)
    if (is.null(pred_df)) return(tibble(rmse = NA, nll = NA, nll_multi = NA, coverage_95 = NA))
    pred_df$row_idx <- seq_len(nrow(pred_df))

    if ("Input_1" %in% names(pred_df) && !"Input" %in% names(pred_df)) {
      pred_df <- pred_df %>% dplyr::rename(Input = Input_1)
    }

    # Filtrer Output 1
    if ("Output_ID" %in% names(pred_df)) {
      pred_df_o1 <- pred_df %>%
        dplyr::mutate(Output_ID = as.character(Output_ID)) %>%
        dplyr::filter(Output_ID == "1") %>%
        dplyr::mutate(Input = round(Input, 5)) %>%
        dplyr::distinct(Input, .keep_all = TRUE)
    } else {
      pred_df_o1 <- pred_df %>%
        dplyr::mutate(Input = round(Input, 5)) %>%
        dplyr::distinct(Input, .keep_all = TRUE)
    }

    mean_col <- intersect(names(pred_df_o1), c("Mean", "mean", "Prediction", "Mu"))
    var_col  <- intersect(names(pred_df_o1), c("Var", "var", "Variance", "Sigma2"))
    if (length(mean_col) == 0 || length(var_col) == 0) {
      return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    }

    # Dé-scaler les prédictions
    sf <- scale_factors %>%
      dplyr::filter(Output_ID == "1", Cluster_ID == cluster_id)
    if (nrow(sf) == 0) {
      # Fallback : essayer sans filtre Cluster_ID (n_clust=1 → un seul cluster)
      sf <- scale_factors %>% dplyr::filter(Output_ID == "1")
    }
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

    # NLL mixture (pred_magmaclust) ou NLL gaussian (pred_magma)
    nll <- NA
    if ("pred_by_cluster" %in% names(pred_res)) {
      pred_by_clus <- pred_res$pred_by_cluster %>%
        dplyr::filter(Output_ID == "1")
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

    # NLL multivariée (matrice de covariance complète)
    nll_multi <- NA
    if ("pred_by_cluster" %in% names(pred_res) && "pred" %in% names(pred_res) && is.list(pred_res$pred)) {
      # Cas pred_magmaclust : NLL mixture multivariée
      cluster_preds_multi <- list()
      truth_vec_multi <- NULL
      for (cl_name in names(pred_res$pred)) {
        cl_obj <- pred_res$pred[[cl_name]]
        cl_cov <- extract_pred_cov(cl_obj)
        cl_df  <- extract_pred_df(cl_obj)
        if (is.null(cl_cov) || is.null(cl_df)) next
        cl_df$row_idx <- seq_len(nrow(cl_df))
        if ("Input_1" %in% names(cl_df) && !"Input" %in% names(cl_df))
          cl_df <- cl_df %>% dplyr::rename(Input = Input_1)
        if ("Output_ID" %in% names(cl_df)) {
          cl_df_o1 <- cl_df %>%
            dplyr::mutate(Output_ID = as.character(Output_ID)) %>%
            dplyr::filter(Output_ID == "1") %>%
            dplyr::mutate(Input = round(Input, 5)) %>%
            dplyr::distinct(Input, .keep_all = TRUE)
        } else {
          cl_df_o1 <- cl_df %>%
            dplyr::mutate(Input = round(Input, 5)) %>%
            dplyr::distinct(Input, .keep_all = TRUE)
        }
        cl_mc <- intersect(names(cl_df_o1), c("Mean", "mean", "Prediction", "Mu"))
        if (length(cl_mc) == 0) next
        cl_aligned <- cl_df_o1 %>%
          dplyr::inner_join(truth_o1_unscaled, by = "Input", suffix = c("_pred", ""))
        if (nrow(cl_aligned) == 0) next
        ci <- cl_aligned$row_idx
        cluster_preds_multi[[cl_name]] <- list(
          mean = cl_aligned[[cl_mc[1]]] * s,
          cov  = cl_cov[ci, ci, drop = FALSE] * s^2
        )
        if (is.null(truth_vec_multi)) truth_vec_multi <- cl_aligned$Output
      }
      if (length(cluster_preds_multi) > 0 && !is.null(truth_vec_multi)) {
        w_m <- if ("weights" %in% names(pred_res)) pred_res$weights else pred_res$mixture
        if (is.data.frame(w_m) && nrow(w_m) == 1) {
          w_m <- as.numeric(w_m); names(w_m) <- names(pred_res$pred)
        }
        nll_multi <- compute_nll_multivariate_mixture(cluster_preds_multi, w_m, truth_vec_multi)
      }
    } else {
      # Cas pred_magma : NLL gaussienne multivariée
      cov_matrix <- extract_pred_cov(pred_res)
      if (!is.null(cov_matrix) && nrow(merged) > 0 && "row_idx" %in% names(merged)) {
        ci <- merged$row_idx
        cov_sub <- cov_matrix[ci, ci, drop = FALSE] * s^2
        nll_multi <- compute_nll_multivariate_gaussian(pred_mean, cov_sub, truth)
      }
    }

    tibble(rmse = rmse, nll = nll, nll_multi = nll_multi, coverage_95 = coverage_95)
  }, error = function(e) {
    cat(paste0("    [ERREUR MOMT] ", conditionMessage(e), "\n"))
    tibble(rmse = NA, nll = NA, nll_multi = NA, coverage_95 = NA)
  })
}

# ---- Extraction des métriques pour MT ----
extract_metrics_mt <- function(pred_entry) {
  tryCatch({
    preds_by_output <- pred_entry$predictions_by_output
    truth_data      <- pred_entry$truth

    # On ne s'intéresse qu'à l'output 1
    if (!("1" %in% names(preds_by_output))) {
      return(tibble(rmse = NA, nll = NA, coverage_95 = NA))
    }

    pred_res <- preds_by_output[["1"]]$prediction
    truth_o1 <- truth_data %>% dplyr::filter(Output_ID == "1")
    if (nrow(truth_o1) == 0) return(tibble(rmse = NA, nll = NA, coverage_95 = NA))

    # Extraire la prédiction (gère pred_magmaclust ET pred_magma)
    pred_df <- extract_pred_df(pred_res)
    if (is.null(pred_df)) return(tibble(rmse = NA, nll = NA, nll_multi = NA, coverage_95 = NA))
    pred_df$row_idx <- seq_len(nrow(pred_df))

    if ("Input_1" %in% names(pred_df) && !"Input" %in% names(pred_df)) {
      pred_df <- pred_df %>% dplyr::rename(Input = Input_1)
    }

    # MT : Output_ID = "1" toujours, pas besoin de filtrer
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

    # NLL mixture (pred_magmaclust) ou NLL gaussian (pred_magma)
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

    # NLL multivariée (matrice de covariance complète)
    nll_multi <- NA
    if ("pred_by_cluster" %in% names(pred_res) && "pred" %in% names(pred_res) && is.list(pred_res$pred)) {
      # Cas pred_magmaclust : NLL mixture multivariée
      cluster_preds_multi <- list()
      truth_vec_multi <- NULL
      for (cl_name in names(pred_res$pred)) {
        cl_obj <- pred_res$pred[[cl_name]]
        cl_cov <- extract_pred_cov(cl_obj)
        cl_df  <- extract_pred_df(cl_obj)
        if (is.null(cl_cov) || is.null(cl_df)) next
        cl_df$row_idx <- seq_len(nrow(cl_df))
        if ("Input_1" %in% names(cl_df) && !"Input" %in% names(cl_df))
          cl_df <- cl_df %>% dplyr::rename(Input = Input_1)
        cl_df <- cl_df %>%
          dplyr::mutate(Input = round(Input, 5)) %>%
          dplyr::distinct(Input, .keep_all = TRUE)
        cl_mc <- intersect(names(cl_df), c("Mean", "mean", "Prediction", "Mu"))
        if (length(cl_mc) == 0) next
        cl_aligned <- cl_df %>%
          dplyr::inner_join(truth_o1_clean, by = "Input", suffix = c("_pred", ""))
        if (nrow(cl_aligned) == 0) next
        ci <- cl_aligned$row_idx
        cluster_preds_multi[[cl_name]] <- list(
          mean = cl_aligned[[cl_mc[1]]],
          cov  = cl_cov[ci, ci, drop = FALSE]
        )
        if (is.null(truth_vec_multi)) truth_vec_multi <- cl_aligned$Output
      }
      if (length(cluster_preds_multi) > 0 && !is.null(truth_vec_multi)) {
        w_m <- if ("weights" %in% names(pred_res)) pred_res$weights else pred_res$mixture
        if (is.data.frame(w_m) && nrow(w_m) == 1) {
          w_m <- as.numeric(w_m); names(w_m) <- names(pred_res$pred)
        }
        nll_multi <- compute_nll_multivariate_mixture(cluster_preds_multi, w_m, truth_vec_multi)
      }
    } else {
      # Cas pred_magma : NLL gaussienne multivariée
      cov_matrix <- extract_pred_cov(pred_res)
      if (!is.null(cov_matrix) && nrow(merged) > 0 && "row_idx" %in% names(merged)) {
        ci <- merged$row_idx
        cov_sub <- cov_matrix[ci, ci, drop = FALSE]
        nll_multi <- compute_nll_multivariate_gaussian(pred_mean, cov_sub, truth)
      }
    }

    tibble(rmse = rmse, nll = nll, nll_multi = nll_multi, coverage_95 = coverage_95)
  }, error = function(e) {
    cat(paste0("    [ERREUR MT] ", conditionMessage(e), "\n"))
    tibble(rmse = NA, nll = NA, nll_multi = NA, coverage_95 = NA)
  })
}

# ==========================================================================
# BOUCLE PRINCIPALE
# ==========================================================================

username   <- Sys.getenv("USER")
base_dir   <- file.path("/scratch", username, "NeurIPS_experiments_forecasting")
output_dir <- file.path(base_dir, "Metrics")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Configurations : 7 configs de base (n_clust=1) + 3 configs clustering
base_configs <- tibble::tribble(
  ~n_out, ~n_train, ~n_pred, ~n_clust,
  2, 30, 1, 1,
  4, 30, 1, 1,
  8, 30, 1, 1,
  2, 15, 1, 1,
  2, 100, 1, 1,
  2, 30, 10, 1,
  2, 30, 100, 1,
  2, 30, 1, 2,
  2, 30, 1, 3,
  2, 30, 1, 4
)

models <- c("MOMT", "MT")
seeds  <- 1:5

cat("=== Calcul des metriques NeurIPS (forecasting uniquement) ===\n")
cat(paste0("  Date      : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste0("  Configs   : ", nrow(base_configs), "\n"))
cat(paste0("  Modeles   : ", paste(models, collapse = ", "), "\n"))
cat(paste0("  Seeds     : ", paste(seeds, collapse = ", "), "\n\n"))

all_metrics <- list()

for (i in 1:nrow(base_configs)) {
  cfg <- base_configs[i, ]
  config_label <- paste0("out", cfg$n_out, "_train", cfg$n_train,
                          "_pred", cfg$n_pred, "_clust", cfg$n_clust)
  dataset_dir  <- file.path(base_dir, "Datasets", config_label)

  for (model in models) {

    for (seed in seeds) {
      pred_dir  <- file.path(base_dir, paste0("Predictions_", model), config_label)
      pred_file <- file.path(pred_dir, paste0("predictions_seed_", seed, ".rds"))

      if (!file.exists(pred_file)) {
        cat(paste0("  [SKIP] ", config_label, "/", model,
                   " seed=", seed, "\n"))
        next
      }

      pred_data <- readRDS(pred_file)

      # Charger les scale_factors si nécessaire (MOMT)
      scale_factors <- NULL
      if (model == "MOMT") {
        ds_file <- file.path(dataset_dir, paste0("datasets_seed_", seed, ".rds"))
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

      # MOMT et MT : temps pred inclut hyperpost
      t_pred <- t_pred + t_hyperpost

      # Extraire métriques par tâche
      predictions <- pred_data$predictions
      task_metrics_list <- list()

      for (tid in names(predictions)) {
        if (model == "MOMT") {
          metrics_task <- extract_metrics_momt(predictions[[tid]], scale_factors)
        } else if (model == "MT") {
          metrics_task <- extract_metrics_mt(predictions[[tid]])
        }
        task_metrics_list[[tid]] <- metrics_task
      }

      task_metrics_df <- bind_rows(task_metrics_list, .id = "task_id")

      rmse_all     <- mean(task_metrics_df$rmse, na.rm = TRUE)
      nll_all      <- mean(task_metrics_df$nll, na.rm = TRUE)
      nll_multi_all <- mean(task_metrics_df$nll_multi, na.rm = TRUE)
      coverage_all <- mean(task_metrics_df$coverage_95, na.rm = TRUE)

      all_metrics[[length(all_metrics) + 1]] <- tibble(
        n_out       = cfg$n_out,
        n_train     = cfg$n_train,
        n_pred      = cfg$n_pred,
        n_clust     = cfg$n_clust,
        problem     = "forecasting",
        model       = model,
        seed        = seed,
        rmse        = rmse_all,
        nll         = nll_all,
        nll_multi   = nll_multi_all,
        coverage_95 = coverage_all,
        t_training  = t_training,
        t_pred      = t_pred,
        t_hyperpost = t_hyperpost
      )

      cat(paste0("  [OK] ", config_label, "/", model,
                 " seed=", seed, " RMSE=", round(rmse_all, 4), "\n"))
    }
  }
}

metrics_raw <- bind_rows(all_metrics)

if (nrow(metrics_raw) == 0) {
  stop("Aucune prediction forecasting exploitable trouvee dans les repertoires configures.")
}

# --- Sauvegarde des métriques brutes ---
write_csv(metrics_raw, file.path(output_dir, "metrics_raw_forecasting_NeurIPS.csv"))
saveRDS(metrics_raw, file.path(output_dir, "metrics_raw_forecasting_NeurIPS.rds"))
cat(paste0("\nMetriques brutes : ", nrow(metrics_raw), " lignes\n"))

# --- Résumé statistique ---
metrics_summary <- metrics_raw %>%
  dplyr::group_by(n_out, n_train, n_pred, n_clust, problem, model) %>%
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

    nll_multi_mean   = mean(nll_multi, na.rm = TRUE),
    nll_multi_std    = sd(nll_multi, na.rm = TRUE),
    nll_multi_median = median(nll_multi, na.rm = TRUE),
    nll_multi_q1     = quantile(nll_multi, 0.25, na.rm = TRUE),
    nll_multi_q3     = quantile(nll_multi, 0.75, na.rm = TRUE),
    nll_multi_iqr    = IQR(nll_multi, na.rm = TRUE),

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
    `RMSE (mean+-std)`     = paste0(round(rmse_mean, 4), " +- ", round(rmse_std, 4)),
    `RMSE (median[IQR])`  = paste0(round(rmse_median, 4), " [", round(rmse_q1, 4), ", ", round(rmse_q3, 4), "]"),
    `NLL (mean+-std)`      = paste0(round(nll_mean, 4), " +- ", round(nll_std, 4)),
    `NLL (median[IQR])`   = paste0(round(nll_median, 4), " [", round(nll_q1, 4), ", ", round(nll_q3, 4), "]"),
    `NLL_Multi (mean+-std)`      = paste0(round(nll_multi_mean, 4), " +- ", round(nll_multi_std, 4)),
    `NLL_Multi (median[IQR])`   = paste0(round(nll_multi_median, 4), " [", round(nll_multi_q1, 4), ", ", round(nll_multi_q3, 4), "]"),
    `Cov95 (mean+-std)`    = paste0(round(cov95_mean, 4), " +- ", round(cov95_std, 4)),
    `Cov95 (median[IQR])` = paste0(round(cov95_median, 4), " [", round(cov95_q1, 4), ", ", round(cov95_q3, 4), "]"),
    `T_train (mean+-std)`  = paste0(round(t_train_mean, 2), " +- ", round(t_train_std, 2)),
    `T_pred (mean+-std)`   = paste0(round(t_pred_mean, 2), " +- ", round(t_pred_std, 2))
  ) %>%
  dplyr::select(
    n_out, n_train, n_pred, n_clust, problem, model,
    `RMSE (mean+-std)`, `RMSE (median[IQR])`,
    `NLL (mean+-std)`, `NLL (median[IQR])`,
    `NLL_Multi (mean+-std)`, `NLL_Multi (median[IQR])`,
    `Cov95 (mean+-std)`, `Cov95 (median[IQR])`,
    `T_train (mean+-std)`, `T_pred (mean+-std)`
  )

metrics_summary <- metrics_summary %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))

write_excel_csv2(metrics_summary, file.path(output_dir, "metrics_summary_forecasting_NeurIPS.csv"))
write_excel_csv2(metrics_display, file.path(output_dir, "metrics_display_forecasting_NeurIPS.csv"))
saveRDS(metrics_summary, file.path(output_dir, "metrics_summary_forecasting_NeurIPS.rds"))

# --- Affichage par paramètre d'intérêt ---
cat("\n--- Resume forecasting par parametre d'interet ---\n\n")

# Nb outputs (n_train=30, n_pred=1, n_clust=1)
cat("=== Parametre : nb_outputs ===\n")
sub <- metrics_display %>%
  dplyr::filter(n_train == 30, n_pred == 1, n_clust == 1) %>%
  dplyr::arrange(n_out, model)
print(as.data.frame(sub), row.names = FALSE)
cat("\n")

# Nb train (n_out=2, n_pred=1, n_clust=1)
cat("=== Parametre : nb_tasks_train ===\n")
sub <- metrics_display %>%
  dplyr::filter(n_out == 2, n_pred == 1, n_clust == 1) %>%
  dplyr::arrange(n_train, model)
print(as.data.frame(sub), row.names = FALSE)
cat("\n")

# Nb pred (n_out=2, n_train=30, n_clust=1)
cat("=== Parametre : nb_tasks_pred ===\n")
sub <- metrics_display %>%
  dplyr::filter(n_out == 2, n_train == 30, n_clust == 1) %>%
  dplyr::arrange(n_pred, model)
print(as.data.frame(sub), row.names = FALSE)
cat("\n")

# Nb clusters (n_out=2, n_train=30, n_pred=1)
cat("=== Parametre : nb_clusters ===\n")
sub <- metrics_display %>%
  dplyr::filter(n_out == 2, n_train == 30, n_pred == 1) %>%
  dplyr::arrange(n_clust, model)
print(as.data.frame(sub), row.names = FALSE)
cat("\n")

cat(paste0("\nFichiers sauvegardes dans : ", output_dir, "\n"))
cat(paste0("Termine : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))