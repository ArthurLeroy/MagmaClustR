#!/bin/bash -l
#===============================================================================
# job_NeurIPS_nclust_phase1_MOMT.sh
#===============================================================================

#SBATCH --job-name=neurips_nclust_momt_v2
#SBATCH --qos=huge
#SBATCH -c 16
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/%u/logs/neurips_nclust_momt_v2_%j.out
#SBATCH --error=/scratch/%u/logs/neurips_nclust_momt_v2_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alexia.grenouillat@math.univ-toulouse.fr

N_WORKERS=16
N_SEEDS=50

N_OUT=2
N_TRAIN=30
N_PRED=1

echo "=============================================="
echo " NeurIPS Phase 1 n_clust : MOMT"
echo " Date    : $(date)"
echo " Noeud   : $(hostname)"
echo " Job ID  : ${SLURM_JOB_ID}"
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
LOGDIR="/scratch/${USER}/logs/neurips_nclust_momt_v2"
mkdir -p "${LOGDIR}"

# --- Installer le package ---
LIB_TEMP_V2="/scratch/${USER}/R_temp_v2_$(date +%Y%m%d)"
mkdir -p "${LIB_TEMP_V2}"
Rscript -e "install.packages('/scratch/${USER}/NeurIPS_experiments_v2', repos = NULL, type = 'source', lib = '${LIB_TEMP_V2}')"
export R_LIBS="${LIB_TEMP_V2}:${R_LIBS}"

SCRIPT_DIR="/scratch/${USER}/NeurIPS_experiments_v2/scripts"
UTILS_DIR="/scratch/${USER}/NeurIPS_experiments_v2/utils"
mkdir -p "${SCRIPT_DIR}" "${UTILS_DIR}"

RSCRIPT="/opt/spack/opt/spack/linux-debian11-zen2/gcc-13.2.0/r-4.4.0-tohpugilej6myswwe73dlbkypu7qqn4p/bin/Rscript"

# =========================================================
# CONFIGURATION DES EXPERIENCES
# Modifie ces tableaux pour faire varier tes valeurs !
# =========================================================
declare -a NOUT_VALUES=(${N_OUT})
declare -a NTRAIN_VALUES=(${N_TRAIN})
declare -a NCLUST_VALUES=(4 3 2 1)
declare -a NPRED_VALUES=(${N_PRED})

# --- Génération de la file de jobs ---
JOBFILE="/scratch/${USER}/logs/neurips_nclust_momt_v2/jobqueue_${SLURM_JOB_ID}.txt"

> "${JOBFILE}"
# ORDRE DES BOUCLES : N_OUT est à l'intérieur, donc N_OUT varie EN PREMIER !
for PROBLEM in interpolation forecasting; do
  for N_CLUST in "${NCLUST_VALUES[@]}"; do
    for N_TRAIN_VAL in "${NTRAIN_VALUES[@]}"; do
      for N_PRED_VAL in "${NPRED_VALUES[@]}"; do
        for N_OUT_VAL in "${NOUT_VALUES[@]}"; do
          for SEED in $(seq 1 ${N_SEEDS}); do
            echo "${N_OUT_VAL} ${N_TRAIN_VAL} ${N_CLUST} ${N_PRED_VAL} ${PROBLEM} ${SEED}" >> "${JOBFILE}"
          done
        done
      done
    done
  done
done

TOTAL_JOBS=$(wc -l < "${JOBFILE}")
echo ""
echo "--- File de jobs : ${TOTAL_JOBS} jobs pour ${N_WORKERS} workers ---"
echo ""

export RSCRIPT SCRIPT_DIR LOGDIR

# --- Lancement via xargs ---
cat "${JOBFILE}" | xargs -n 6 -P ${N_WORKERS} bash -c '
  N_OUT_VAL=$1
  N_TRAIN_VAL=$2
  N_CLUST=$3
  N_PRED_VAL=$4
  PROBLEM=$5
  SEED=$6
  
  LABEL="out${N_OUT_VAL}_train${N_TRAIN_VAL}_clust${N_CLUST}_pred${N_PRED_VAL}_${PROBLEM}_seed${SEED}"
  LOGFILE="${LOGDIR}/momt_${LABEL}.log"
  
  echo "[$(date +%H:%M:%S)] Démarrage → MOMT ${LABEL}"
  
  ${RSCRIPT} --vanilla "${SCRIPT_DIR}/Benchmark_NeurIPS_MOMT_nclust.R" \
    --n_out=${N_OUT_VAL} --n_train=${N_TRAIN_VAL} --n_pred=${N_PRED_VAL} \
    --n_clust=${N_CLUST} --problem=${PROBLEM} --seed=${SEED} \
    > "${LOGFILE}" 2>&1
    
  EXIT_CODE=$?
  if [ ${EXIT_CODE} -ne 0 ]; then
    echo "[ERREUR] MOMT ${LABEL} a échoué (code=${EXIT_CODE})"
  fi
' _

echo ""
echo "=============================================="
echo " Phase 1 n_clust MOMT TERMINÉE"
echo " Date    : $(date)"
echo "=============================================="

# --- Vérification des fichiers ---
RESULTS_DIR="/scratch/${USER}/NeurIPS_experiments_v2"
MISSING=0

for PROBLEM in interpolation forecasting; do
  for N_CLUST in "${NCLUST_VALUES[@]}"; do
    for N_TRAIN_VAL in "${NTRAIN_VALUES[@]}"; do
      for N_PRED_VAL in "${NPRED_VALUES[@]}"; do
        for N_OUT_VAL in "${NOUT_VALUES[@]}"; do
          for SEED in $(seq 1 ${N_SEEDS}); do
            DATASET_FILE="${RESULTS_DIR}/out${N_OUT_VAL}_train${N_TRAIN_VAL}_pred${N_PRED_VAL}_clust${N_CLUST}/${PROBLEM}/Datasets/datasets_seed_${SEED}.rds"
            if [ ! -f "${DATASET_FILE}" ]; then
              echo "[MANQUANT] ${DATASET_FILE}"
              MISSING=$((MISSING + 1))
            fi
          done
        done
      done
    done
  done
done

echo ""
echo "Fichiers manquants : ${MISSING} / ${TOTAL_JOBS}"
if [ ${MISSING} -gt 0 ]; then
  exit 1
fi
exit 0