#!/bin/bash -l
#===============================================================================
# job_phase2_mt_mo.sh — PHASE 2 : Entraînement MT + MO
#
# Lance 30 processus R en parallèle :
#   - 15 processus MT (5 n_out × 3 n_train)
#   - 15 processus MO (5 n_out × 3 n_train)
#
# Queue "huge" : 32 threads, 7 jours max
# 30 processus R ≈ 30 cœurs utilisés
#
# IMPORTANT : lancer UNIQUEMENT après avoir vérifié les résultats de la Phase 1.
#
# Soumission : sbatch job_phase2_mt_mo.sh
#===============================================================================

#SBATCH --job-name=xp1_phase2
#SBATCH --qos=huge
#SBATCH -c 32
#SBATCH --output=/scratch/%u/logs/phase2_mt_mo_%j.out
#SBATCH --error=/scratch/%u/logs/phase2_mt_mo_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alexia.grenouillat@math.univ-toulouse.fr

echo "=============================================="
echo " PHASE 2 : MT + MO"
echo " Date    : $(date)"
echo " Noeud   : $(hostname)"
echo " Job ID  : ${SLURM_JOB_ID}"
echo " CPUs    : ${SLURM_CPUS_ON_NODE}"
echo " 30 processus R en parallèle (15 MT + 15 MO)"
echo "=============================================="

# --- Environnement ---
source ~/.bashrc             # Force le chargement de ton environnement personnel
shopt -s expand_aliases      # Active les alias au cas où load_spack en soit un

load_spack
spack load r@4.4.0 cmake
export R_LIBS=/scratch/${USER}/R
# # --- Environnement ---
# load_spack
# spack load r@4.4.0 cmake
# export R_LIBS=/scratch/${USER}/R

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd /scratch/${USER}/MagmaClustR

# --- Vérification préalable : les données MOMT existent-elles ? ---
RESULTS_DIR="/scratch/${USER}/NeurIPS_experiments/Experience_1"
MISSING=0

for N_OUT in 3 4 6 8; do
  for N_TRAIN in 15 30 300; do
    for N_PRED in 1 10 100 1000; do
      DATASET_DIR="${RESULTS_DIR}/Datasets/n_out_${N_OUT}/train_${N_TRAIN}_pred_${N_PRED}"
      if [ ! -d "${DATASET_DIR}" ]; then
        echo "[ATTENTION] Répertoire manquant : ${DATASET_DIR}"
        MISSING=$((MISSING + 1))
      fi
    done
  done
done

if [ ${MISSING} -gt 0 ]; then
  echo ""
  echo "ERREUR : ${MISSING} répertoires de données MOMT manquants."
  echo "Vérifiez que la Phase 1 s'est terminée correctement."
  echo "Abandon."
  exit 1
fi

echo "Vérification OK : toutes les données MOMT sont présentes."
echo ""

# --- Logs ---
LOGDIR="/scratch/${USER}/logs/phase2_details"
mkdir -p "${LOGDIR}"

PIDS=()

# --- Lancer les processus MT (5 n_out × 3 n_train = 15) ---
# Note : n_out=2 est un cas "single output" qui ne dépend pas de MOMT.
# Le script MT gère ce cas en créant ses propres données si nécessaire.
for N_OUT in 2 3 4 6 8; do
  for N_TRAIN in 15 30 300; do
    LOGFILE="${LOGDIR}/mt_nout${N_OUT}_ntrain${N_TRAIN}.log"
    echo "[$(date +%H:%M:%S)] Lancement MT  n_out=${N_OUT} n_train=${N_TRAIN} → ${LOGFILE}"

    Rscript --vanilla Benchmark_XP_1_MT_cluster.R \
      --n_out=${N_OUT} --n_train=${N_TRAIN} \
      > "${LOGFILE}" 2>&1 &

    PIDS+=($!)
  done
done

# --- Lancer les processus MO (5 n_out × 3 n_train = 15) ---
for N_OUT in 2 3 4 6 8; do
  for N_TRAIN in 15 30 300; do
    LOGFILE="${LOGDIR}/mo_nout${N_OUT}_ntrain${N_TRAIN}.log"
    echo "[$(date +%H:%M:%S)] Lancement MO  n_out=${N_OUT} n_train=${N_TRAIN} → ${LOGFILE}"

    Rscript --vanilla Benchmark_XP_1_MO_cluster.R \
      --n_out=${N_OUT} --n_train=${N_TRAIN} \
      > "${LOGFILE}" 2>&1 &

    PIDS+=($!)
  done
done

echo ""
echo "30 processus lancés. PIDs : ${PIDS[*]}"
echo "En attente de la fin de tous les processus..."
echo ""

# --- Attendre et collecter les résultats ---
FAILED=0
for i in "${!PIDS[@]}"; do
  wait ${PIDS[$i]}
  EXIT_CODE=$?
  if [ ${EXIT_CODE} -ne 0 ]; then
    echo "[ERREUR] Processus PID=${PIDS[$i]} terminé avec code ${EXIT_CODE}"
    FAILED=$((FAILED + 1))
  fi
done

echo ""
echo "=============================================="
echo " PHASE 2 TERMINÉE"
echo " Date  : $(date)"
echo " Échecs: ${FAILED} / 30"
echo "=============================================="

if [ ${FAILED} -gt 0 ]; then
  echo "ATTENTION : ${FAILED} processus ont échoué. Vérifiez les logs dans ${LOGDIR}/"
  exit 1
fi

exit 0
