#!/bin/bash -l
#===============================================================================
# job_NeurIPS_phase1_MOMT.sh
# NeurIPS : Phase 1 — Entraînement MOMT pour les 7 configs × 2 problèmes
#
# 14 expériences MOMT en parallèle (14 cœurs).
# Chaque expérience enchaîne les 5 seeds séquentiellement.
#
# Soumission : sbatch job_NeurIPS_phase1_MOMT.sh
#===============================================================================

#SBATCH --job-name=neurips_momt
#SBATCH --qos=large
#SBATCH -c 14
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/%u/logs/neurips_momt_%j.out
#SBATCH --error=/scratch/%u/logs/neurips_momt_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alexia.grenouillat@math.univ-toulouse.fr

echo "=============================================="
echo " NeurIPS Phase 1 : MOMT (14 configs × 5 seeds)"
echo " Date    : $(date)"
echo " Noeud   : $(hostname)"
echo " Job ID  : ${SLURM_JOB_ID}"
echo " CPUs    : ${SLURM_CPUS_ON_NODE}"
echo "=============================================="

# --- Environnement ---
source ~/.bashrc
load_spack
spack load r@4.4.0
export R_LIBS=/scratch/${USER}/R

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd /scratch/${USER}/MagmaClustR

# --- Logs ---
LOGDIR="/scratch/${USER}/logs/neurips_momt"
mkdir -p "${LOGDIR}"

# --- Installer le package ---
LIB_TEMP="/scratch/${USER}/R_temp_$(date +%Y%m%d)"
mkdir -p "${LIB_TEMP}"
Rscript -e "install.packages('/scratch/${USER}/MagmaClustR', repos = NULL, type = 'source', lib = '${LIB_TEMP}')"
export R_LIBS="${LIB_TEMP}:${R_LIBS}"

# --- Copier les scripts ---
SCRIPT_DIR="/scratch/${USER}/NeurIPS_experiments/scripts"
UTILS_DIR="/scratch/${USER}/NeurIPS_experiments/utils"
mkdir -p "${SCRIPT_DIR}" "${UTILS_DIR}"

# S'assurer que les scripts sont sur le serveur
# (à copier manuellement ou via rsync avant soumission)

RSCRIPT="/opt/spack/opt/spack/linux-debian11-zen2/gcc-13.2.0/r-4.4.0-tohpugilej6myswwe73dlbkypu7qqn4p/bin/Rscript"

# Les 7 configurations uniques (n_out, n_train, n_pred)
# Config 1 (défaut) : 2 30 1
# Config 2 : 4 30 1
# Config 3 : 8 30 1
# Config 4 : 2 15 1
# Config 5 : 2 100 1
# Config 6 : 2 30 10
# Config 7 : 2 30 100

declare -a CONFIGS=(
  "2 30 1"
  "4 30 1"
  "8 30 1"
  "2 15 1"
  "2 100 1"
  "2 30 10"
  "2 30 100"
)

echo ""
echo "--- Phase 1 : MOMT (14 processus = 7 configs × 2 problèmes) ---"
echo ""

PIDS=()

for PROBLEM in interpolation forecasting; do
  for CONFIG_STR in "${CONFIGS[@]}"; do
    read -r N_OUT N_TRAIN N_PRED <<< "${CONFIG_STR}"
    LABEL="out${N_OUT}_train${N_TRAIN}_pred${N_PRED}_${PROBLEM}"
    LOGFILE="${LOGDIR}/momt_${LABEL}.log"

    echo "[$(date +%H:%M:%S)] MOMT ${LABEL} → ${LOGFILE}"

    # Lancer les 5 seeds séquentiellement dans un sous-shell
    (
      for SEED in $(seq 1 5); do
        echo "[$(date +%H:%M:%S)] MOMT ${LABEL} seed=${SEED}" >> "${LOGFILE}"
        ${RSCRIPT} --vanilla "${SCRIPT_DIR}/Benchmark_NeurIPS_MOMT.R" \
          --n_out=${N_OUT} --n_train=${N_TRAIN} --n_pred=${N_PRED} \
          --problem=${PROBLEM} --seed=${SEED} \
          >> "${LOGFILE}" 2>&1

        EXIT_CODE=$?
        if [ ${EXIT_CODE} -ne 0 ]; then
          echo "[ERREUR] MOMT ${LABEL} seed=${SEED} code=${EXIT_CODE}" >> "${LOGFILE}"
        fi
      done
    ) &

    PIDS+=($!)
  done
done

echo ""
echo "${#PIDS[@]} processus MOMT lancés. PIDs : ${PIDS[*]}"
echo "En attente..."

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
echo " Phase 1 MOMT TERMINÉE"
echo " Date    : $(date)"
echo " Échecs  : ${FAILED} / ${#PIDS[@]}"
echo "=============================================="

# --- Vérification des fichiers ---
RESULTS_DIR="/scratch/${USER}/NeurIPS_experiments"
MISSING=0

for PROBLEM in interpolation forecasting; do
  for CONFIG_STR in "${CONFIGS[@]}"; do
    read -r N_OUT N_TRAIN N_PRED <<< "${CONFIG_STR}"
    for SEED in $(seq 1 5); do
      DATASET_FILE="${RESULTS_DIR}/out${N_OUT}_train${N_TRAIN}_pred${N_PRED}/${PROBLEM}/Datasets/datasets_seed_${SEED}.rds"
      if [ ! -f "${DATASET_FILE}" ]; then
        echo "[MANQUANT] ${DATASET_FILE}"
        MISSING=$((MISSING + 1))
      fi
    done
  done
done

echo ""
echo "Fichiers manquants : ${MISSING}"
if [ ${MISSING} -gt 0 ]; then
  echo "ATTENTION : des fichiers MOMT sont manquants. Vérifiez les logs."
  exit 1
fi

echo "Vérification OK : toutes les données MOMT sont présentes."
exit 0
