#!/bin/bash -l
#===============================================================================
# job_NeurIPS_phase2_MT_MO.sh
# NeurIPS : Phase 2 — MT et MO pour les 7 configs × 2 problèmes
#
# Pré-requis : les données MOMT doivent être présentes (Phase 1 terminée).
#
# 14 MT + 14 MO = 28 processus lancés en parallèle (28 cœurs).
#
# Soumission : sbatch job_NeurIPS_phase2_MT_MO.sh
#===============================================================================

#SBATCH --job-name=neurips_mt_mo
#SBATCH --qos=huge
#SBATCH -c 28
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/%u/logs/neurips_mt_mo_%j.out
#SBATCH --error=/scratch/%u/logs/neurips_mt_mo_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alexia.grenouillat@math.univ-toulouse.fr

echo "=============================================="
echo " NeurIPS Phase 2 : MT + MO"
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
LOGDIR="/scratch/${USER}/logs/neurips_mt_mo"
mkdir -p "${LOGDIR}"

# --- Installer le package ---
LIB_TEMP="/scratch/${USER}/R_temp_$(date +%Y%m%d)"
mkdir -p "${LIB_TEMP}"
Rscript -e "install.packages('/scratch/${USER}/MagmaClustR', repos = NULL, type = 'source', lib = '${LIB_TEMP}')"
export R_LIBS="${LIB_TEMP}:${R_LIBS}"

RSCRIPT="/opt/spack/opt/spack/linux-debian11-zen2/gcc-13.2.0/r-4.4.0-tohpugilej6myswwe73dlbkypu7qqn4p/bin/Rscript"
SCRIPT_DIR="/scratch/${USER}/NeurIPS_experiments/scripts"

declare -a CONFIGS=(
  "2 30 1"
  "4 30 1"
  "8 30 1"
  "2 15 1"
  "2 100 1"
  "2 30 10"
  "2 30 100"
)

# ==========================================
# LANCEMENT MT + MO EN PARALLÈLE (28 processus)
# ==========================================
echo ""
echo "--- Lancement simultané : 14 MT + 14 MO = 28 processus ---"
echo ""

PIDS_MT=()
PIDS_MO=()

for PROBLEM in interpolation forecasting; do
  for CONFIG_STR in "${CONFIGS[@]}"; do
    read -r N_OUT N_TRAIN N_PRED <<< "${CONFIG_STR}"
    LABEL="out${N_OUT}_train${N_TRAIN}_pred${N_PRED}_${PROBLEM}"

    # --- MT ---
    LOGFILE_MT="${LOGDIR}/mt_${LABEL}.log"
    echo "[$(date +%H:%M:%S)] MT ${LABEL} → ${LOGFILE_MT}"

    (
      for SEED in $(seq 1 5); do
        echo "[$(date +%H:%M:%S)] MT ${LABEL} seed=${SEED}" >> "${LOGFILE_MT}"
        ${RSCRIPT} --vanilla "${SCRIPT_DIR}/Benchmark_NeurIPS_MT.R" \
          --n_out=${N_OUT} --n_train=${N_TRAIN} --n_pred=${N_PRED} \
          --problem=${PROBLEM} --seed=${SEED} \
          >> "${LOGFILE_MT}" 2>&1

        EXIT_CODE=$?
        if [ ${EXIT_CODE} -ne 0 ]; then
          echo "[ERREUR] MT ${LABEL} seed=${SEED} code=${EXIT_CODE}" >> "${LOGFILE_MT}"
        fi
      done
    ) &
    PIDS_MT+=($!)

    # --- MO ---
    LOGFILE_MO="${LOGDIR}/mo_${LABEL}.log"
    echo "[$(date +%H:%M:%S)] MO ${LABEL} → ${LOGFILE_MO}"

    (
      for SEED in $(seq 1 5); do
        echo "[$(date +%H:%M:%S)] MO ${LABEL} seed=${SEED}" >> "${LOGFILE_MO}"
        ${RSCRIPT} --vanilla "${SCRIPT_DIR}/Benchmark_NeurIPS_MO.R" \
          --n_out=${N_OUT} --n_train=${N_TRAIN} --n_pred=${N_PRED} \
          --problem=${PROBLEM} --seed=${SEED} \
          >> "${LOGFILE_MO}" 2>&1

        EXIT_CODE=$?
        if [ ${EXIT_CODE} -ne 0 ]; then
          echo "[ERREUR] MO ${LABEL} seed=${SEED} code=${EXIT_CODE}" >> "${LOGFILE_MO}"
        fi
      done
    ) &
    PIDS_MO+=($!)

  done
done

echo ""
echo "${#PIDS_MT[@]} processus MT + ${#PIDS_MO[@]} processus MO lancés. En attente..."

# --- Attente MT ---
FAILED_MT=0
for i in "${!PIDS_MT[@]}"; do
  wait ${PIDS_MT[$i]}
  EXIT_CODE=$?
  if [ ${EXIT_CODE} -ne 0 ]; then
    echo "[ERREUR] MT PID=${PIDS_MT[$i]} terminé avec code ${EXIT_CODE}"
    FAILED_MT=$((FAILED_MT + 1))
  fi
done

# --- Attente MO ---
FAILED_MO=0
for i in "${!PIDS_MO[@]}"; do
  wait ${PIDS_MO[$i]}
  EXIT_CODE=$?
  if [ ${EXIT_CODE} -ne 0 ]; then
    echo "[ERREUR] MO PID=${PIDS_MO[$i]} terminé avec code ${EXIT_CODE}"
    FAILED_MO=$((FAILED_MO + 1))
  fi
done

# ==========================================
# RÉSUMÉ
# ==========================================
echo ""
echo "=============================================="
echo " NeurIPS Phase 2 TERMINÉE"
echo " Date       : $(date)"
echo " MT         : ${FAILED_MT} / 14 échecs"
echo " MO         : ${FAILED_MO} / 14 échecs"
echo "=============================================="

TOTAL_FAILED=$((FAILED_MT + FAILED_MO))
if [ ${TOTAL_FAILED} -gt 0 ]; then
  echo "ATTENTION : ${TOTAL_FAILED} processus ont échoué. Vérifiez ${LOGDIR}/"
  exit 1
fi

exit 0
