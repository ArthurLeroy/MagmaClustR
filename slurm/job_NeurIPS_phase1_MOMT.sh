#!/bin/bash -l
#===============================================================================
# job_NeurIPS_phase1_MOMT.sh
# NeurIPS : Phase 1 — Entraînement MOMT pour les 7 configs × 2 problèmes
#
# 700 jobs (7 configs × 2 problèmes × 50 seeds) répartis sur 32 cœurs.
# Chaque cœur pioche le prochain job disponible dans une file d'attente.
#
# Soumission : sbatch job_NeurIPS_phase1_MOMT.sh
#===============================================================================

#SBATCH --job-name=neurips_momt
#SBATCH --qos=huge
#SBATCH -c 32
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/%u/logs/neurips_momt_%j.out
#SBATCH --error=/scratch/%u/logs/neurips_momt_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alexia.grenouillat@math.univ-toulouse.fr

N_WORKERS=32
N_SEEDS=50

echo "=============================================="
echo " NeurIPS Phase 1 : MOMT (7 configs × 2 problèmes × ${N_SEEDS} seeds)"
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
LOGDIR="/scratch/${USER}/logs/neurips_momt"
mkdir -p "${LOGDIR}"

# --- Installer le package ---
LIB_TEMP="/scratch/${USER}/R_temp_$(date +%Y%m%d)"
mkdir -p "${LIB_TEMP}"
Rscript -e "install.packages('/scratch/${USER}/MagmaClustR', repos = NULL, type = 'source', lib = '${LIB_TEMP}')"
export R_LIBS="${LIB_TEMP}:${R_LIBS}"

SCRIPT_DIR="/scratch/${USER}/NeurIPS_experiments/scripts"
UTILS_DIR="/scratch/${USER}/NeurIPS_experiments/utils"
mkdir -p "${SCRIPT_DIR}" "${UTILS_DIR}"

RSCRIPT="/opt/spack/opt/spack/linux-debian11-zen2/gcc-13.2.0/r-4.4.0-tohpugilej6myswwe73dlbkypu7qqn4p/bin/Rscript"

# Les 7 configurations uniques (n_out, n_train, n_pred)
# Ordre choisi : les configurations MOMT les plus lourdes en premier.
declare -a CONFIGS=(
  "4 30 1"
  "8 30 1"
  "2 30 1"
  "2 15 1"
  "2 100 1"
  "2 30 10"
  "2 30 100"
)

# --- Génération de la file de jobs ---
JOBFILE="/scratch/${USER}/logs/neurips_momt/jobqueue_${SLURM_JOB_ID}.txt"
COUNTERFILE="/scratch/${USER}/logs/neurips_momt/counter_${SLURM_JOB_ID}"

> "${JOBFILE}"
for PROBLEM in interpolation forecasting; do
  for CONFIG_STR in "${CONFIGS[@]}"; do
    read -r N_OUT N_TRAIN N_PRED <<< "${CONFIG_STR}"
    for SEED in $(seq 1 ${N_SEEDS}); do
      echo "${N_OUT} ${N_TRAIN} ${N_PRED} ${PROBLEM} ${SEED}" >> "${JOBFILE}"
    done
  done
done

TOTAL_JOBS=$(wc -l < "${JOBFILE}")
echo "0" > "${COUNTERFILE}"

echo ""
echo "--- File de jobs : ${TOTAL_JOBS} jobs pour ${N_WORKERS} workers ---"
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
    read -r N_OUT N_TRAIN N_PRED PROBLEM SEED <<< "${LINE}"
    LABEL="out${N_OUT}_train${N_TRAIN}_pred${N_PRED}_${PROBLEM}_seed${SEED}"
    LOGFILE="${LOGDIR}/momt_${LABEL}.log"

    echo "[$(date +%H:%M:%S)] Worker ${WORKER_ID} → MOMT ${LABEL}"

    ${RSCRIPT} --vanilla "${SCRIPT_DIR}/Benchmark_NeurIPS_MOMT.R" \
      --n_out=${N_OUT} --n_train=${N_TRAIN} --n_pred=${N_PRED} \
      --problem=${PROBLEM} --seed=${SEED} \
      > "${LOGFILE}" 2>&1

    EXIT_CODE=$?
    if [ ${EXIT_CODE} -ne 0 ]; then
      echo "[ERREUR] Worker ${WORKER_ID} : MOMT ${LABEL} code=${EXIT_CODE}"
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

echo ""
echo "=============================================="
echo " Phase 1 MOMT TERMINÉE"
echo " Date    : $(date)"
echo " Échecs workers : ${FAILED} / ${N_WORKERS}"
echo "=============================================="

# --- Vérification des fichiers ---
RESULTS_DIR="/scratch/${USER}/NeurIPS_experiments"
MISSING=0

for PROBLEM in interpolation forecasting; do
  for CONFIG_STR in "${CONFIGS[@]}"; do
    read -r N_OUT N_TRAIN N_PRED <<< "${CONFIG_STR}"
    for SEED in $(seq 1 ${N_SEEDS}); do
      DATASET_FILE="${RESULTS_DIR}/out${N_OUT}_train${N_TRAIN}_pred${N_PRED}/${PROBLEM}/Datasets/datasets_seed_${SEED}.rds"
      if [ ! -f "${DATASET_FILE}" ]; then
        echo "[MANQUANT] ${DATASET_FILE}"
        MISSING=$((MISSING + 1))
      fi
    done
  done
done

echo ""
echo "Fichiers manquants : ${MISSING} / ${TOTAL_JOBS}"
if [ ${MISSING} -gt 0 ]; then
  echo "ATTENTION : des fichiers MOMT sont manquants. Vérifiez les logs dans ${LOGDIR}/"
  exit 1
fi

echo "Vérification OK : toutes les données MOMT sont présentes."
exit 0
