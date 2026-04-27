#!/bin/bash -l
#===============================================================================
# job_NeurIPS_nclust_rerun_predictions_MOMT_interpolation.sh
# NeurIPS : relance uniquement les prédictions MOMT interpolation.
#
# 10 configs x 50 seeds = 500 jobs
#
# Soumission : sbatch job_NeurIPS_nclust_rerun_predictions_MOMT_interpolation.sh
#===============================================================================

#SBATCH --job-name=rerun_momt_itp
#SBATCH --qos=huge
#SBATCH -c 32
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/%u/logs/rerun_momt_itp_%j.out
#SBATCH --error=/scratch/%u/logs/rerun_momt_itp_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alexia.grenouillat@math.univ-toulouse.fr

N_WORKERS=32
N_SEEDS=50
PROJECT_ROOT="/scratch/${USER}/NeurIPS_experiments_interpolation"
SCRIPT_DIR="${PROJECT_ROOT}/scripts"
SCRIPT_NAME="Rerun_Predictions_NeurIPS_MOMT_nclust_interpolation.R"
RESULTS_DIR="${PROJECT_ROOT}/Predictions_MOMT"
LOGDIR="/scratch/${USER}/logs/rerun_momt_itp"

echo "=============================================="
echo " RERUN MOMT INTERPOLATION"
echo " Date    : $(date)"
echo " Noeud   : $(hostname)"
echo " Job ID  : ${SLURM_JOB_ID}"
echo " CPUs    : ${SLURM_CPUS_ON_NODE}"
echo " Workers : ${N_WORKERS}"
echo "=============================================="

source ~/.bashrc
load_spack
spack load r@4.4.0
export R_LIBS=/scratch/${USER}/R

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

mkdir -p "${LOGDIR}"

if [ ! -d "${SCRIPT_DIR}" ]; then
  echo "[ERREUR] Répertoire scripts introuvable : ${SCRIPT_DIR}"
  exit 1
fi

RSCRIPT="$(command -v Rscript)"
if [ -z "${RSCRIPT}" ]; then
  echo "[ERREUR] Rscript introuvable dans l'environnement courant."
  exit 1
fi

declare -a CONFIGS=(
  "4 30 1 1"
  "8 30 1 1"
  "2 30 1 1"
  "2 15 1 1"
  "2 100 1 1"
  "2 30 10 1"
  "2 30 100 1"
  "2 30 1 2"
  "2 30 1 3"
  "2 30 1 4"
)

JOBFILE="${LOGDIR}/jobqueue_${SLURM_JOB_ID}.txt"
COUNTERFILE="${LOGDIR}/counter_${SLURM_JOB_ID}"

> "${JOBFILE}"
for CONFIG_STR in "${CONFIGS[@]}"; do
  read -r N_OUT N_TRAIN N_PRED N_CLUST <<< "${CONFIG_STR}"
  for SEED in $(seq 1 ${N_SEEDS}); do
    echo "${N_OUT} ${N_TRAIN} ${N_PRED} ${N_CLUST} interpolation ${SEED}" >> "${JOBFILE}"
  done
done

TOTAL_JOBS=$(wc -l < "${JOBFILE}")
echo "0" > "${COUNTERFILE}"

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
    LABEL="momt_itp_out${N_OUT}_train${N_TRAIN}_pred${N_PRED}_clust${N_CLUST}_seed${SEED}"
    LOGFILE="${LOGDIR}/${LABEL}.log"

    echo "[$(date +%H:%M:%S)] Worker ${WORKER_ID} -> ${LABEL}"

    "${RSCRIPT}" --vanilla "${SCRIPT_DIR}/${SCRIPT_NAME}" \
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

PIDS=()
for W in $(seq 1 ${N_WORKERS}); do
  worker ${W} &
  PIDS+=($!)
done

FAILED=0
for PID in "${PIDS[@]}"; do
  wait ${PID}
  EXIT_CODE=$?
  if [ ${EXIT_CODE} -ne 0 ]; then
    FAILED=$((FAILED + 1))
  fi
done

MISSING=0
while read -r N_OUT N_TRAIN N_PRED N_CLUST PROBLEM SEED; do
  PRED_FILE="${RESULTS_DIR}/out${N_OUT}_train${N_TRAIN}_pred${N_PRED}_clust${N_CLUST}/predictions_seed_${SEED}.rds"
  if [ ! -f "${PRED_FILE}" ]; then
    echo "[MANQUANT] ${PRED_FILE}"
    MISSING=$((MISSING + 1))
  fi
done < "${JOBFILE}"

echo "Échecs workers : ${FAILED} / ${N_WORKERS}"
echo "Fichiers manquants : ${MISSING} / ${TOTAL_JOBS}"
if [ ${MISSING} -gt 0 ]; then
  exit 1
fi
exit 0