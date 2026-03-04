#!/bin/bash -l
#===============================================================================
# job_phase1_momt.sh — PHASE 1 : Entraînement MOMT (several_outputs)
#
# Lance 15 processus R en parallèle (5 n_out × 3 n_train).
# Chaque processus exécute 100 itérations pour une combinaison (n_out, n_train).
#
# Queue "huge" : 32 threads, 7 jours max
# 15 processus R ≈ 15 cœurs utilisés (R est single-threaded par processus)
#
# Soumission : sbatch job_phase1_momt.sh
#===============================================================================

#SBATCH --job-name=xp1_phase1
#SBATCH --qos=huge
#SBATCH -c 32
#SBATCH --output=/scratch/%u/logs/phase1_momt_%j.out
#SBATCH --error=/scratch/%u/logs/phase1_momt_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alexia.grenouillat@math.univ-toulouse.fr

echo "=============================================="
echo " PHASE 1 : MOMT (several_outputs)"
echo " Date    : $(date)"
echo " Noeud   : $(hostname)"
echo " Job ID  : ${SLURM_JOB_ID}"
echo " CPUs    : ${SLURM_CPUS_ON_NODE}"
echo " 15 processus R en parallèle"
echo "=============================================="

# --- Environnement ---
# En mode batch, le shell n'est ni login ni interactif :
# les fichiers ~/.bashrc et ~/.bash_profile ne sont pas lus.
# On les source explicitement pour rendre load_spack disponible.
source ~/.bashrc

load_spack
spack load r@4.4.0
export R_LIBS=/scratch/${USER}/R

# Limiter BLAS à 1 thread par processus (éviter la surcharge)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd /scratch/${USER}/MagmaClustR

# --- Lancer les 15 combinaisons en parallèle ---
LOGDIR="/scratch/${USER}/logs/phase1_details"
mkdir -p "${LOGDIR}"

PIDS=()

for N_OUT in 2 3 4 6 8; do
  for N_TRAIN in 15 30 300; do
    LOGFILE="${LOGDIR}/momt_nout${N_OUT}_ntrain${N_TRAIN}.log"
    echo "[$(date +%H:%M:%S)] Lancement MOMT n_out=${N_OUT} n_train=${N_TRAIN} → ${LOGFILE}"

    stdbuf -oL Rscript --vanilla Benchmark_XP_1_several_outputs_cluster.R \
      --n_out=${N_OUT} --n_train=${N_TRAIN} \
      > "${LOGFILE}" 2>&1 &

    PIDS+=($!)
  done
done

echo ""
echo "15 processus lancés. PIDs : ${PIDS[*]}"
echo "En attente de la fin de tous les processus..."
echo ""

# --- Attendre tous les processus et collecter les codes de retour ---
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
echo " PHASE 1 TERMINÉE"
echo " Date  : $(date)"
echo " Échecs: ${FAILED} / 15"
echo "=============================================="

if [ ${FAILED} -gt 0 ]; then
  echo "ATTENTION : ${FAILED} processus ont échoué. Vérifiez les logs dans ${LOGDIR}/"
  exit 1
fi

exit 0
