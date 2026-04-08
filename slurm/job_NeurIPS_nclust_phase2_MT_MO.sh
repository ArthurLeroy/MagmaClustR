#!/bin/bash -l
#===============================================================================
# job_NeurIPS_nclust_phase2_MT_MO.sh
# NeurIPS : Phase 2 n_clust — MT et MO pour n_clust=1,2,3,4
#===============================================================================

#SBATCH --job-name=neurips_nclust_mt_mo_v2
#SBATCH --qos=huge
#SBATCH -c 16
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/%u/logs/neurips_nclust_mt_mo_v2_%j.out
#SBATCH --error=/scratch/%u/logs/neurips_nclust_mt_mo_v2_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alexia.grenouillat@math.univ-toulouse.fr

N_WORKERS=16
N_SEEDS=50

# Configuration fixe pour l'étude n_clust
N_OUT=2
N_TRAIN=30
N_PRED=1

echo "=============================================="
echo " NeurIPS Phase 2 n_clust : MT + MO"
echo " n_clust = 1, 2, 3, 4"
echo " Config fixe : out=${N_OUT}, train=${N_TRAIN}, pred=${N_PRED}"
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
LOGDIR="/scratch/${USER}/logs/neurips_nclust_mt_mo_v2"
mkdir -p "${LOGDIR}"

# --- Installer le package (Version 2) ---
LIB_TEMP_V2="/scratch/${USER}/R_temp_v2_$(date +%Y%m%d)"
mkdir -p "${LIB_TEMP_V2}"
Rscript -e "install.packages('/scratch/${USER}/NeurIPS_experiments_v2', repos = NULL, type = 'source', lib = '${LIB_TEMP_V2}')"
export R_LIBS="${LIB_TEMP_V2}:${R_LIBS}"

RSCRIPT="/opt/spack/opt/spack/linux-debian11-zen2/gcc-13.2.0/r-4.4.0-tohpugilej6myswwe73dlbkypu7qqn4p/bin/Rscript"
SCRIPT_DIR="/scratch/${USER}/NeurIPS_experiments_v2/scripts"

# Valeurs de n_clust (du plus lourd au plus léger)
declare -a NCLUST_VALUES=(4 3 2 1)

# --- Génération de la file de jobs ---
JOBFILE="/scratch/${USER}/logs/neurips_nclust_mt_mo_v2/jobqueue_${SLURM_JOB_ID}.txt"

> "${JOBFILE}"

# Jobs MT
for PROBLEM in interpolation forecasting; do
  for N_CLUST in "${NCLUST_VALUES[@]}"; do
    for SEED in $(seq 1 ${N_SEEDS}); do
      echo "MT ${N_CLUST} ${PROBLEM} ${SEED}" >> "${JOBFILE}"
    done
  done
done

# Jobs MO
for PROBLEM in interpolation forecasting; do
  for N_CLUST in "${NCLUST_VALUES[@]}"; do
    for SEED in $(seq 1 ${N_SEEDS}); do
      echo "MO ${N_CLUST} ${PROBLEM} ${SEED}" >> "${JOBFILE}"
    done
  done
done

TOTAL_JOBS=$(wc -l < "${JOBFILE}")

echo ""
echo "--- File de jobs : ${TOTAL_JOBS} jobs (MT + MO) pour ${N_WORKERS} workers ---"
echo ""

# On exporte les variables fixes pour que le sous-processus xargs puisse les lire
export RSCRIPT SCRIPT_DIR LOGDIR N_OUT N_TRAIN N_PRED

# --- Lancement via xargs (Robuste, sans collisions) ---
cat "${JOBFILE}" | xargs -n 4 -P ${N_WORKERS} bash -c '
  MODEL=$1
  N_CLUST=$2
  PROBLEM=$3
  SEED=$4
  
  LABEL="${MODEL}_clust${N_CLUST}_${PROBLEM}_seed${SEED}"
  LOGFILE="${LOGDIR}/${LABEL}.log"
  
  echo "[$(date +%H:%M:%S)] Démarrage → ${LABEL}"
  
  if [ "${MODEL}" = "MT" ]; then
    ${RSCRIPT} --vanilla "${SCRIPT_DIR}/Benchmark_NeurIPS_MT_nclust.R" \
      --n_out=${N_OUT} --n_train=${N_TRAIN} --n_pred=${N_PRED} \
      --n_clust=${N_CLUST} --problem=${PROBLEM} --seed=${SEED} \
      > "${LOGFILE}" 2>&1
  else
    ${RSCRIPT} --vanilla "${SCRIPT_DIR}/Benchmark_NeurIPS_MO_nclust.R" \
      --n_out=${N_OUT} --n_train=${N_TRAIN} --n_pred=${N_PRED} \
      --n_clust=${N_CLUST} --problem=${PROBLEM} --seed=${SEED} \
      > "${LOGFILE}" 2>&1
  fi
    
  EXIT_CODE=$?
  if [ ${EXIT_CODE} -ne 0 ]; then
    echo "[ERREUR] ${LABEL} a échoué (code=${EXIT_CODE})"
  fi
' _

echo ""
echo "=============================================="
echo " NeurIPS Phase 2 n_clust TERMINÉE"
echo " Date       : $(date)"
echo "=============================================="

# --- Vérification des fichiers ---
RESULTS_DIR="/scratch/${USER}/NeurIPS_experiments_v2"
MISSING=0

# Vérification MT
for PROBLEM in interpolation forecasting; do
  for N_CLUST in "${NCLUST_VALUES[@]}"; do
    for SEED in $(seq 1 ${N_SEEDS}); do
      PRED_FILE="${RESULTS_DIR}/out${N_OUT}_train${N_TRAIN}_pred${N_PRED}_clust${N_CLUST}/${PROBLEM}/Predictions_MT/predictions_seed_${SEED}.rds"
      if [ ! -f "${PRED_FILE}" ]; then
        echo "[MANQUANT MT] ${PRED_FILE}"
        MISSING=$((MISSING + 1))
      fi
    done
  done
done

# Vérification MO
for PROBLEM in interpolation forecasting; do
  for N_CLUST in "${NCLUST_VALUES[@]}"; do
    for SEED in $(seq 1 ${N_SEEDS}); do
      PRED_FILE="${RESULTS_DIR}/out${N_OUT}_train${N_TRAIN}_pred${N_PRED}_clust${N_CLUST}/${PROBLEM}/Predictions_MO/predictions_seed_${SEED}.rds"
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

echo "Vérification OK : toutes les prédictions MT/MO n_clust sont présentes."
exit 0
