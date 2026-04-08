#!/bin/bash -l
#===============================================================================
# job_NeurIPS_nclust_phase1_MOMT.sh
# NeurIPS : Phase 1 n_clust — Entraînement MOMT pour n_clust=1,2,3,4
#
# Configuration fixe : n_out=2, n_train=30, n_pred=1
# 4 valeurs de n_clust × 2 problèmes × 50 seeds = 400 jobs
#
# Soumission : sbatch job_NeurIPS_nclust_phase1_MOMT.sh
#===============================================================================

#SBATCH --job-name=neurips_nclust_momt
#SBATCH --qos=huge
#SBATCH -c 16
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/%u/logs/neurips_nclust_momt_%j.out
#SBATCH --error=/scratch/%u/logs/neurips_nclust_momt_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alexia.grenouillat@math.univ-toulouse.fr

N_WORKERS=16
N_SEEDS=50

# Configuration fixe pour l'étude n_clust
N_OUT=2
N_TRAIN=30
N_PRED=1

echo "=============================================="
echo " NeurIPS Phase 1 n_clust : MOMT"
echo " n_clust = 1, 2, 3, 4"
echo " Config fixe : out=${N_OUT}, train=${N_TRAIN}, pred=${N_PRED}"
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

cd /scratch/${USER}/NeurIPS_experiments_v2

# --- Logs ---
LOGDIR="/scratch/${USER}/logs/neurips_nclust_momt"
mkdir -p "${LOGDIR}"

# --- Installer le package (Version 2) ---
LIB_TEMP_V2="/scratch/${USER}/R_temp_v2_$(date +%Y%m%d)"
mkdir -p "${LIB_TEMP_V2}"
# Pointez vers le nouveau dossier source que vous avez cloné
Rscript -e "install.packages('/scratch/${USER}/NeurIPS_experiments_v2', repos = NULL, type = 'source', lib = '${LIB_TEMP_V2}')"
# On force R à chercher d'abord dans cette librairie
export R_LIBS="${LIB_TEMP_V2}:${R_LIBS}"

SCRIPT_DIR="/scratch/${USER}/NeurIPS_experiments/scripts"
UTILS_DIR="/scratch/${USER}/NeurIPS_experiments/utils"
mkdir -p "${SCRIPT_DIR}" "${UTILS_DIR}"

RSCRIPT="/opt/spack/opt/spack/linux-debian11-zen2/gcc-13.2.0/r-4.4.0-tohpugilej6myswwe73dlbkypu7qqn4p/bin/Rscript"

# Valeurs de n_clust à tester (du plus lourd au plus léger)
declare -a NCLUST_VALUES=(4 3 2 1)

# --- Génération de la file de jobs ---
JOBFILE="/scratch/${USER}/logs/neurips_nclust_momt/jobqueue_${SLURM_JOB_ID}.txt"
COUNTERFILE="/scratch/${USER}/logs/neurips_nclust_momt/counter_${SLURM_JOB_ID}"

> "${JOBFILE}"
for PROBLEM in interpolation forecasting; do
  for N_CLUST in "${NCLUST_VALUES[@]}"; do
    for SEED in $(seq 1 ${N_SEEDS}); do
      echo "${N_CLUST} ${PROBLEM} ${SEED}" >> "${JOBFILE}"
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
    read -r N_CLUST PROBLEM SEED <<< "${LINE}"
    LABEL="clust${N_CLUST}_${PROBLEM}_seed${SEED}"
    LOGFILE="${LOGDIR}/momt_${LABEL}.log"

    echo "[$(date +%H:%M:%S)] Worker ${WORKER_ID} → MOMT ${LABEL}"

    ${RSCRIPT} --vanilla "${SCRIPT_DIR}/Benchmark_NeurIPS_MOMT_nclust.R" \
      --n_out=${N_OUT} --n_train=${N_TRAIN} --n_pred=${N_PRED} \
      --n_clust=${N_CLUST} --problem=${PROBLEM} --seed=${SEED} \
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
echo " Phase 1 n_clust MOMT TERMINÉE"
echo " Date    : $(date)"
echo " Échecs workers : ${FAILED} / ${N_WORKERS}"
echo "=============================================="

# --- Vérification des fichiers ---
RESULTS_DIR="/scratch/${USER}/NeurIPS_experiments_v2"
MISSING=0

for PROBLEM in interpolation forecasting; do
  for N_CLUST in "${NCLUST_VALUES[@]}"; do
    for SEED in $(seq 1 ${N_SEEDS}); do
      DATASET_FILE="${RESULTS_DIR}/out${N_OUT}_train${N_TRAIN}_pred${N_PRED}_clust${N_CLUST}/${PROBLEM}/Datasets/datasets_seed_${SEED}.rds"
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

echo "Vérification OK : toutes les données MOMT n_clust sont présentes."
exit 0
