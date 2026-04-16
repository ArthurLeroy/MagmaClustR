#!/bin/bash -l
#===============================================================================
# job_NeurIPS_nclust_phase2_MT_interpolation.sh
# NeurIPS : Phase 2 — MT INTERPOLATION ONLY
#
# 10 configs × 1 problème (interpolation) × 5 seeds = 50 jobs
#
# Soumission : sbatch job_NeurIPS_nclust_phase2_MT_interpolation.sh
#===============================================================================

#SBATCH --job-name=neurips_nclust_mt_itp
#SBATCH --qos=huge
#SBATCH -c 16
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/%u/logs/neurips_nclust_mt_itp_%j.out
#SBATCH --error=/scratch/%u/logs/neurips_nclust_mt_itp_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alexia.grenouillat@math.univ-toulouse.fr

N_WORKERS=16
N_SEEDS=5

echo "=============================================="
echo " NeurIPS Phase 2 : MT INTERPOLATION (10 configs)"
echo " Date    : $(date)"
echo " Noeud   : $(hostname)"
echo " Job ID  : ${SLURM_JOB_ID}"
echo " CPUs    : ${SLURM_CPUS_ON_NODE}"
echo " Workers : ${N_WORKERS}"
echo "=============================================="

# --- Environnement ---
source ~/.bashrc
load_spack
spack load r@4.4.0
export R_LIBS=/scratch/${USER}/R

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd /scratch/${USER}/NeurIPS_experiments_interpolation

# --- Logs ---
LOGDIR="/scratch/${USER}/logs/neurips_nclust_mt_itp"
mkdir -p "${LOGDIR}"

# --- Installer le package ---
LIB_TEMP_V2="/scratch/${USER}/R_temp_v2_$(date +%Y%m%d)"
mkdir -p "${LIB_TEMP_V2}"
Rscript -e "install.packages('/scratch/${USER}/NeurIPS_experiments_interpolation', repos = NULL, type = 'source', lib = '${LIB_TEMP_V2}')"
export R_LIBS="${LIB_TEMP_V2}:${R_LIBS}"

RSCRIPT="/opt/spack/opt/spack/linux-debian11-zen2/gcc-13.2.0/r-4.4.0-tohpugilej6myswwe73dlbkypu7qqn4p/bin/Rscript"
SCRIPT_DIR="/scratch/${USER}/NeurIPS_experiments_interpolation/scripts"

declare -a MT_CONFIGS=(
  "2 30 1 2"
  "2 30 1 3"
  "2 30 1 4"
  "4 30 1 1"
  "8 30 1 1"
  "2 30 1 1"
  "2 15 1 1"
  "2 100 1 1"
  "2 30 10 1"
  "2 30 100 1"
)

# --- Génération de la file de jobs ---
JOBFILE="/scratch/${USER}/logs/neurips_nclust_mt_itp/jobqueue_${SLURM_JOB_ID}.txt"
COUNTERFILE="/scratch/${USER}/logs/neurips_nclust_mt_itp/counter_${SLURM_JOB_ID}"

> "${JOBFILE}"
for CONFIG_STR in "${MT_CONFIGS[@]}"; do
  read -r N_OUT N_TRAIN N_PRED N_CLUST <<< "${CONFIG_STR}"
  for SEED in $(seq 1 ${N_SEEDS}); do
    echo "${N_OUT} ${N_TRAIN} ${N_PRED} ${N_CLUST} interpolation ${SEED}" >> "${JOBFILE}"
  done
done

TOTAL_JOBS=$(wc -l < "${JOBFILE}")
echo "0" > "${COUNTERFILE}"

echo ""
echo "--- File de jobs : ${TOTAL_JOBS} jobs MT pour ${N_WORKERS} workers ---"
echo ""

# --- Fonction worker ---
worker() {
  local WORKER_ID=$1
  local JOB_NUM
  local LINE

  exec 200>"${COUNTERFILE}.lock"
  while true; do
    flock -x 200
    JOB_NUM=$(cat "${COUNTERFILE}")
    echo $((JOB_NUM + 1)) > "${COUNTERFILE}"
    flock -u 200

    if [ "${JOB_NUM}" -ge "${TOTAL_JOBS}" ]; then
      break
    fi

    LINE=$(sed -n "$((JOB_NUM + 1))p" "${JOBFILE}")
    read -r N_OUT N_TRAIN N_PRED N_CLUST PROBLEM SEED <<< "${LINE}"
    LABEL="MT_out${N_OUT}_train${N_TRAIN}_pred${N_PRED}_clust${N_CLUST}_${PROBLEM}_seed${SEED}"
    LOGFILE="${LOGDIR}/${LABEL}.log"

    echo "[$(date +%H:%M:%S)] Worker ${WORKER_ID} → ${LABEL}"

    ${RSCRIPT} --vanilla "${SCRIPT_DIR}/Benchmark_NeurIPS_MT_nclust_interpolation.R" \
      --n_out=${N_OUT} --n_train=${N_TRAIN} --n_pred=${N_PRED} --n_clust=${N_CLUST} \
      --problem=${PROBLEM} --seed=${SEED} \
      > "${LOGFILE}" 2>&1

    EXIT_CODE=$?
    if [ ${EXIT_CODE} -ne 0 ]; then
      echo "[ERREUR] Worker ${WORKER_ID} : ${LABEL} code=${EXIT_CODE}"
    fi
  done
  exec 200>&-
}

# --- Lancement des workers ---
PIDS=()
for W in $(seq 1 ${N_WORKERS}); do
  worker ${W} &
  PIDS+=($!)
done

echo "${#PIDS[@]} workers lancés. PIDs : ${PIDS[*]}"
echo "En attente de la fin de tous les jobs..."

FAILED=0
for PID in "${PIDS[@]}"; do
  wait ${PID}
  EXIT_CODE=$?
  if [ ${EXIT_CODE} -ne 0 ]; then
    echo "[ERREUR] Worker PID=${PID} terminé avec code ${EXIT_CODE}"
    FAILED=$((FAILED + 1))
  fi
done

echo ""
echo "=============================================="
echo " NeurIPS Phase 2 MT INTERPOLATION TERMINÉE"
echo " Date       : $(date)"
echo " Échecs workers : ${FAILED} / ${N_WORKERS}"
echo "=============================================="

RESULTS_DIR="/scratch/${USER}/NeurIPS_experiments_interpolation"
MISSING=0

for CONFIG_STR in "${MT_CONFIGS[@]}"; do
  read -r N_OUT N_TRAIN N_PRED N_CLUST <<< "${CONFIG_STR}"
  for SEED in $(seq 1 ${N_SEEDS}); do
    PRED_FILE="${RESULTS_DIR}/Predictions_MT/out${N_OUT}_train${N_TRAIN}_pred${N_PRED}_clust${N_CLUST}/predictions_seed_${SEED}.rds"
    if [ ! -f "${PRED_FILE}" ]; then
      echo "[MANQUANT MT] ${PRED_FILE}"
      MISSING=$((MISSING + 1))
    fi
  done
done

echo ""
echo "Fichiers manquants : ${MISSING} / ${TOTAL_JOBS}"
if [ ${MISSING} -gt 0 ]; then
  echo "ATTENTION : des fichiers MT sont manquants. Vérifiez ${LOGDIR}/"
  exit 1
fi

echo "Vérification OK : toutes les prédictions MT interpolation sont présentes."
exit 0