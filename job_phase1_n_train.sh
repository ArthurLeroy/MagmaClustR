#!/bin/bash -l
#===============================================================================
# job_phase1_n_train.sh — Phase 1 : Paramètre "n_train" (20 vs 100 tâches)
#
# Phase 1 : 20 MOMT (5 seeds × 2 configs)
# Phase 2a: 20 MT   (5 seeds × 2 configs)
# Phase 2b: 20 MO   (5 seeds × 2 configs)
#
# Soumission : sbatch job_phase1_n_train.sh
#===============================================================================

#SBATCH --job-name=p1_n_train
#SBATCH --qos=huge
#SBATCH -c 20
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/%u/logs/phase1_n_train_%j.out
#SBATCH --error=/scratch/%u/logs/phase1_n_train_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alexia.grenouillat@math.univ-toulouse.fr

PARAM="n_train"

echo "=============================================="
echo " Phase 1 : Paramètre '${PARAM}'"
echo " Date    : $(date)"
echo " Noeud   : $(hostname)"
echo " Job ID  : ${SLURM_JOB_ID}"
echo " CPUs    : ${SLURM_CPUS_ON_NODE}"
echo " 20 cœurs : 5 seeds × 2 configs"
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
LOGDIR="/scratch/${USER}/logs/phase1_${PARAM}"
mkdir -p "${LOGDIR}"

# --- Installer le package ---
LIB_TEMP="/scratch/${USER}/R_temp_$(date +%Y%m%d)"
mkdir -p "${LIB_TEMP}"
Rscript -e "install.packages('/scratch/${USER}/MagmaClustR', repos = NULL, type = 'source', lib = '${LIB_TEMP}')"
export R_LIBS="${LIB_TEMP}:${R_LIBS}"

RSCRIPT="/opt/spack/opt/spack/linux-debian11-zen2/gcc-13.2.0/r-4.4.0-tohpugilej6myswwe73dlbkypu7qqn4p/bin/Rscript"

# ==========================================
# PHASE 1 : MOMT (20 processus)
# ==========================================
echo ""
echo "--- PHASE 1 : MOMT (20 processus) ---"
echo ""

PIDS=()

for CONFIG in default variation; do
  for SEED in $(seq 1 5); do
    LOGFILE="${LOGDIR}/momt_${CONFIG}_seed${SEED}.log"
    echo "[$(date +%H:%M:%S)] MOMT config=${CONFIG} seed=${SEED} → ${LOGFILE}"

    ${RSCRIPT} --vanilla Benchmark_Phase1_MOMT.R \
      --param=${PARAM} --config=${CONFIG} --seed=${SEED} \
      > "${LOGFILE}" 2>&1 &

    PIDS+=($!)
  done
done

echo ""
echo "${#PIDS[@]} processus MOMT lancés. PIDs : ${PIDS[*]}"
echo "En attente..."

FAILED_P1=0
for i in "${!PIDS[@]}"; do
  wait ${PIDS[$i]}
  EXIT_CODE=$?
  if [ ${EXIT_CODE} -ne 0 ]; then
    echo "[ERREUR] MOMT PID=${PIDS[$i]} terminé avec code ${EXIT_CODE}"
    FAILED_P1=$((FAILED_P1 + 1))
  fi
done

echo ""
echo "--- PHASE 1 TERMINÉE (Échecs: ${FAILED_P1} / 10) ---"

if [ ${FAILED_P1} -gt 0 ]; then
  echo "ATTENTION : ${FAILED_P1} processus MOMT ont échoué. Vérifiez ${LOGDIR}/"
  echo "Abandon des phases suivantes."
  exit 1
fi

# --- Vérification ---
RESULTS_DIR="/scratch/${USER}/Phase1_experiments/${PARAM}"
MISSING=0

for CONFIG in default variation; do
  for SEED in $(seq 1 5); do
    DATASET_FILE="${RESULTS_DIR}/${CONFIG}/Datasets/datasets_seed_${SEED}.rds"
    if [ ! -f "${DATASET_FILE}" ]; then
      echo "[ATTENTION] Fichier manquant : ${DATASET_FILE}"
      MISSING=$((MISSING + 1))
    fi
  done
done

if [ ${MISSING} -gt 0 ]; then
  echo "ERREUR : ${MISSING} fichiers de données manquants. Abandon."
  exit 1
fi

echo "Vérification OK : toutes les données MOMT sont présentes."

# ==========================================
# PHASE 2a : MT (20 processus)
# ==========================================
echo ""
echo "--- PHASE 2a : MT (20 processus) ---"
echo ""

PIDS=()

for CONFIG in default variation; do
  for SEED in $(seq 1 5); do
    LOGFILE="${LOGDIR}/mt_${CONFIG}_seed${SEED}.log"
    echo "[$(date +%H:%M:%S)] MT config=${CONFIG} seed=${SEED} → ${LOGFILE}"

    ${RSCRIPT} --vanilla Benchmark_Phase1_MT.R \
      --param=${PARAM} --config=${CONFIG} --seed=${SEED} \
      > "${LOGFILE}" 2>&1 &

    PIDS+=($!)
  done
done

echo ""
echo "${#PIDS[@]} processus MT lancés. En attente..."

FAILED_P2A=0
for i in "${!PIDS[@]}"; do
  wait ${PIDS[$i]}
  EXIT_CODE=$?
  if [ ${EXIT_CODE} -ne 0 ]; then
    echo "[ERREUR] MT PID=${PIDS[$i]} terminé avec code ${EXIT_CODE}"
    FAILED_P2A=$((FAILED_P2A + 1))
  fi
done

echo "--- PHASE 2a TERMINÉE (Échecs: ${FAILED_P2A} / 20) ---"

# ==========================================
# PHASE 2b : MO (20 processus)
# ==========================================
echo ""
echo "--- PHASE 2b : MO (20 processus) ---"
echo ""

PIDS=()

for CONFIG in default variation; do
  for SEED in $(seq 1 5); do
    LOGFILE="${LOGDIR}/mo_${CONFIG}_seed${SEED}.log"
    echo "[$(date +%H:%M:%S)] MO config=${CONFIG} seed=${SEED} → ${LOGFILE}"

    ${RSCRIPT} --vanilla Benchmark_Phase1_MO.R \
      --param=${PARAM} --config=${CONFIG} --seed=${SEED} \
      > "${LOGFILE}" 2>&1 &

    PIDS+=($!)
  done
done

echo ""
echo "${#PIDS[@]} processus MO lancés. En attente..."

FAILED_P2B=0
for i in "${!PIDS[@]}"; do
  wait ${PIDS[$i]}
  EXIT_CODE=$?
  if [ ${EXIT_CODE} -ne 0 ]; then
    echo "[ERREUR] MO PID=${PIDS[$i]} terminé avec code ${EXIT_CODE}"
    FAILED_P2B=$((FAILED_P2B + 1))
  fi
done

echo "--- PHASE 2b TERMINÉE (Échecs: ${FAILED_P2B} / 20) ---"

# ==========================================
# RÉSUMÉ
# ==========================================
echo ""
echo "=============================================="
echo " Phase 1 '${PARAM}' TERMINÉE"
echo " Date       : $(date)"
echo " MOMT       : ${FAILED_P1} / 20 échecs"
echo " MT         : ${FAILED_P2A} / 20 échecs"
echo " MO         : ${FAILED_P2B} / 20 échecs"
echo "=============================================="

TOTAL_FAILED=$((FAILED_P1 + FAILED_P2A + FAILED_P2B))
if [ ${TOTAL_FAILED} -gt 0 ]; then
  echo "ATTENTION : ${TOTAL_FAILED} processus ont échoué. Vérifiez ${LOGDIR}/"
  exit 1
fi

exit 0
