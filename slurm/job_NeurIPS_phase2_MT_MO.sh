#!/bin/bash -l
#===============================================================================
# job_NeurIPS_phase2_MT_MO.sh
# NeurIPS : Phase 2 — MT et MO pour les 7 configs × 2 problèmes
#
# Pré-requis : les données MOMT doivent être présentes (Phase 1 terminée).
#
# MT : 7 configs × 2 problèmes × 50 seeds = 700 jobs
# MO : 4 configs × 2 problèmes × seeds variables = 320 jobs
#      (n_out=4,n_train=30,n_pred=1 → 10 seeds ; les autres → 50 seeds)
#      (configs retirées : 8/30/1, 2/100/1, 2/30/100)
#
# Total : 1020 jobs répartis sur 32 cœurs via file d'attente.
#
# Soumission : sbatch job_NeurIPS_phase2_MT_MO.sh
#===============================================================================

#SBATCH --job-name=neurips_mt_mo
#SBATCH --qos=huge
#SBATCH -c 32
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/%u/logs/neurips_mt_mo_%j.out
#SBATCH --error=/scratch/%u/logs/neurips_mt_mo_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alexia.grenouillat@math.univ-toulouse.fr

N_WORKERS=32

echo "=============================================="
echo " NeurIPS Phase 2 : MT + MO"
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

cd /scratch/${USER}/MagmaClustR

# --- Logs ---
LOGDIR="/scratch/${USER}/logs/neurips_mt_mo"
mkdir -p "${LOGDIR}"

# --- Installer le package ---
LIB_TEMP="/scratch/${USER}/R_temp_$(date +%Y%m%d)"
mkdir -p "${LIB_TEMP}"
Rscript -e "install.packages('/scratch/${USER}/MagmaClustR', repos = NULL, type = 'source', lib = '${LIB_TEMP}')"
export R_LIBS="${LIB_TEMP}:${R_LIBS}"

RSCRIPT="/opt/spack/opt/spack/linux-debian11-zen2/gcc-13.2.0/r-4.4.0-tohpugilej6myswwe73dlbkypu7qqn4p/bin/Rscript"
SCRIPT_DIR="/scratch/${USER}/NeurIPS_experiments/scripts"

# --- Configurations MT (7 configs × 50 seeds) ---
# Ordre choisi : les configurations MT les plus lourdes en premier.
declare -a MT_CONFIGS=(
  "4 30 1"
  "2 100 1"
  "8 30 1"
  "2 30 1"
  "2 15 1"
  "2 30 10"
  "2 30 100"
)

# --- Configurations MO (4 configs, seeds variables) ---
# Ordre choisi : les configurations MO les plus lourdes en premier.
#   Format : "n_out n_train n_pred n_seeds"
declare -a MO_CONFIGS=(
  "4 30 1 10"
  "2 30 10 50"
  "2 30 1 50"
  "2 15 1 50"
)

# --- Génération de la file de jobs ---
JOBFILE="/scratch/${USER}/logs/neurips_mt_mo/jobqueue_${SLURM_JOB_ID}.txt"
COUNTERFILE="/scratch/${USER}/logs/neurips_mt_mo/counter_${SLURM_JOB_ID}"

> "${JOBFILE}"

# Jobs MT : 7 configs × 2 problèmes × 50 seeds = 700
for PROBLEM in interpolation forecasting; do
  for CONFIG_STR in "${MT_CONFIGS[@]}"; do
    read -r N_OUT N_TRAIN N_PRED <<< "${CONFIG_STR}"
    for SEED in $(seq 1 50); do
      echo "MT ${N_OUT} ${N_TRAIN} ${N_PRED} ${PROBLEM} ${SEED}" >> "${JOBFILE}"
    done
  done
done

# Jobs MO : 4 configs × 2 problèmes × seeds variables = 320
for PROBLEM in interpolation forecasting; do
  for CONFIG_STR in "${MO_CONFIGS[@]}"; do
    read -r N_OUT N_TRAIN N_PRED N_SEEDS <<< "${CONFIG_STR}"
    for SEED in $(seq 1 ${N_SEEDS}); do
      echo "MO ${N_OUT} ${N_TRAIN} ${N_PRED} ${PROBLEM} ${SEED}" >> "${JOBFILE}"
    done
  done
done

TOTAL_JOBS=$(wc -l < "${JOBFILE}")
echo "0" > "${COUNTERFILE}"

echo ""
echo "--- File de jobs : ${TOTAL_JOBS} jobs (MT + MO) pour ${N_WORKERS} workers ---"
echo ""

# --- Fonction worker ---
worker() {
  local WORKER_ID=$1
  while true; do
    # Lecture atomique du prochain numéro de job
    JOB_NUM=$(
      flock 200
      N=$(cat "${COUNTERFILE}")
      echo $((N + 1)) > "${COUNTERFILE}"
      echo "${N}"
    ) 200>"${COUNTERFILE}.lock"

    if [ "${JOB_NUM}" -ge "${TOTAL_JOBS}" ]; then
      break
    fi

    LINE=$(sed -n "$((JOB_NUM + 1))p" "${JOBFILE}")
    read -r MODEL N_OUT N_TRAIN N_PRED PROBLEM SEED <<< "${LINE}"
    LABEL="${MODEL}_out${N_OUT}_train${N_TRAIN}_pred${N_PRED}_${PROBLEM}_seed${SEED}"
    LOGFILE="${LOGDIR}/${LABEL}.log"

    echo "[$(date +%H:%M:%S)] Worker ${WORKER_ID} → ${LABEL}"

    if [ "${MODEL}" = "MT" ]; then
      ${RSCRIPT} --vanilla "${SCRIPT_DIR}/Benchmark_NeurIPS_MT.R" \
        --n_out=${N_OUT} --n_train=${N_TRAIN} --n_pred=${N_PRED} \
        --problem=${PROBLEM} --seed=${SEED} \
        > "${LOGFILE}" 2>&1
    else
      ${RSCRIPT} --vanilla "${SCRIPT_DIR}/Benchmark_NeurIPS_MO.R" \
        --n_out=${N_OUT} --n_train=${N_TRAIN} --n_pred=${N_PRED} \
        --problem=${PROBLEM} --seed=${SEED} \
        > "${LOGFILE}" 2>&1
    fi

    EXIT_CODE=$?
    if [ ${EXIT_CODE} -ne 0 ]; then
      echo "[ERREUR] Worker ${WORKER_ID} : ${LABEL} code=${EXIT_CODE}"
    fi
  done
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

# ==========================================
# RÉSUMÉ ET VÉRIFICATION
# ==========================================
echo ""
echo "=============================================="
echo " NeurIPS Phase 2 TERMINÉE"
echo " Date       : $(date)"
echo " Échecs workers : ${FAILED} / ${N_WORKERS}"
echo "=============================================="

RESULTS_DIR="/scratch/${USER}/NeurIPS_experiments"
MISSING=0

# Vérification MT
for PROBLEM in interpolation forecasting; do
  for CONFIG_STR in "${MT_CONFIGS[@]}"; do
    read -r N_OUT N_TRAIN N_PRED <<< "${CONFIG_STR}"
    for SEED in $(seq 1 50); do
      PRED_FILE="${RESULTS_DIR}/out${N_OUT}_train${N_TRAIN}_pred${N_PRED}/${PROBLEM}/Predictions_MT/predictions_seed_${SEED}.rds"
      if [ ! -f "${PRED_FILE}" ]; then
        echo "[MANQUANT MT] ${PRED_FILE}"
        MISSING=$((MISSING + 1))
      fi
    done
  done
done

# Vérification MO
for PROBLEM in interpolation forecasting; do
  for CONFIG_STR in "${MO_CONFIGS[@]}"; do
    read -r N_OUT N_TRAIN N_PRED N_SEEDS <<< "${CONFIG_STR}"
    for SEED in $(seq 1 ${N_SEEDS}); do
      PRED_FILE="${RESULTS_DIR}/out${N_OUT}_train${N_TRAIN}_pred${N_PRED}/${PROBLEM}/Predictions_MO/predictions_seed_${SEED}.rds"
      if [ ! -f "${PRED_FILE}" ]; then
        echo "[MANQUANT MO] ${PRED_FILE}"
        MISSING=$((MISSING + 1))
      fi
    done
  done
done

echo ""
echo "Fichiers manquants : ${MISSING} / ${TOTAL_JOBS}"
if [ ${MISSING} -gt 0 ]; then
  echo "ATTENTION : des fichiers MT/MO sont manquants. Vérifiez ${LOGDIR}/"
  exit 1
fi

echo "Vérification OK : toutes les prédictions MT/MO sont présentes."
exit 0
