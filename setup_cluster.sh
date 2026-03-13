#!/bin/bash
#===============================================================================
# setup_cluster.sh — Configuration initiale du cluster IMT
#
# Ce script prépare l'environnement :
#   1. Crée l'arborescence sur /scratch
#   2. Installe les packages R nécessaires
#   3. Vérifie que MagmaClustR se charge correctement
#
# Usage : bash setup_cluster.sh
#===============================================================================

set -e  # Arrêter en cas d'erreur

echo "=============================================="
echo " Configuration du cluster IMT pour XP1"
echo " Utilisateur : ${USER}"
echo " Date        : $(date)"
echo "=============================================="

# --- 1. Charger l'environnement ---
echo ""
echo ">>> 1. Chargement de Spack et R..."
load_spack
spack load r@4.4.0 cmake
echo "R version :"
R --version | head -1

# --- 2. Créer l'arborescence ---
echo ""
echo ">>> 2. Création de l'arborescence sur /scratch..."

# Répertoire des packages R
mkdir -p /scratch/${USER}/R
export R_LIBS=/scratch/${USER}/R
echo "R_LIBS=${R_LIBS}"

# Répertoire du code source
mkdir -p /scratch/${USER}/MagmaClustR_dev

# Répertoire des résultats
RESULTS="/scratch/${USER}/NeurIPS_experiments/Experience_1"
for SUBDIR in Datasets Models_MOMT Predictions_MOMT Models_MT Predictions_MT Models_MO Predictions_MO; do
  for N_OUT in 2 3 4 6 8; do
    for N_TRAIN in 15 30 300; do
      for N_PRED in 1 10 100 1000; do
        mkdir -p "${RESULTS}/${SUBDIR}/n_out_${N_OUT}/train_${N_TRAIN}_pred_${N_PRED}"
      done
    done
  done
done
echo "Arborescence des résultats créée."

# Répertoire des logs
mkdir -p /scratch/${USER}/logs/phase1_details
mkdir -p /scratch/${USER}/logs/phase2_details

# Répertoire run_info (métadonnées JSON)
mkdir -p "${RESULTS}/run_info"
echo "Répertoires de logs et run_info créés."

# --- 3. Installer les packages R ---
echo ""
echo ">>> 3. Installation des packages R..."

R --vanilla -e '
  # Fixer le miroir CRAN
  options(repos = c(CRAN = "https://cloud.r-project.org"))

  # Liste des packages nécessaires
  pkgs <- c(
    "devtools", "Rcpp", "tidyverse", "mvtnorm", "Metrics",
    "matrixStats", "broom", "plyr", "purrr", "rlang",
    "tidyselect", "readr", "jsonlite", "MASS", "stringr"
  )

  # Installer uniquement ceux qui manquent
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    cat("Installation de :", paste(missing, collapse = ", "), "\n")
    install.packages(missing, lib = Sys.getenv("R_LIBS"))
  } else {
    cat("Tous les packages sont déjà installés.\n")
  }

  # Vérifier que tout est chargeable
  cat("\nVérification du chargement des packages :\n")
  for (pkg in pkgs) {
    ok <- requireNamespace(pkg, quietly = TRUE)
    cat(sprintf("  %-15s : %s\n", pkg, ifelse(ok, "OK", "ÉCHEC")))
  }
'

# --- 4. Tester le chargement de MagmaClustR ---
echo ""
echo ">>> 4. Test de chargement de MagmaClustR (devtools::load_all)..."

cd /scratch/${USER}/MagmaClustR_dev

R --vanilla -e '
  library(devtools)
  tryCatch({
    load_all()
    cat("MagmaClustR chargé avec succès !\n")
    cat("Fonctions disponibles : train_magmaclust, pred_magmaclust, etc.\n")
  }, error = function(e) {
    cat("ERREUR lors du chargement de MagmaClustR :\n")
    cat(e$message, "\n")
    cat("\nVérifiez que le code source est bien copié dans /scratch/", Sys.getenv("USER"), "/MagmaClustR_dev/\n")
    quit(status = 1)
  })
'

echo ""
echo "=============================================="
echo " Configuration terminée avec succès !"
echo "=============================================="
echo ""
echo "Prochaines étapes :"
echo "  1. Copier les scripts cluster dans /scratch/${USER}/MagmaClustR_dev/"
echo "  2. Soumettre Phase 1 : sbatch job_phase1_momt.sh"
echo "  3. Vérifier les résultats de Phase 1"
echo "  4. Soumettre Phase 2 : sbatch job_phase2_mt_mo.sh"
