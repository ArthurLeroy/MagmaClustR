#!/bin/bash -l
#===============================================================================
# job_phase1_momt.sh — PHASE 1 : Entraînement MOMT (several_outputs)
#
# Lance 12 processus R en parallèle (4 n_out × 3 n_train).
# Chaque processus exécute 100 itérations pour une combinaison (n_out, n_train).
#
# Queue "huge" : 16 threads, 7 jours max
# 12 processus R ≈ 12 cœurs utilisés (R est single-threaded par processus)
#
# Soumission : sbatch job_phase1_momt.sh
#===============================================================================

#SBATCH --job-name=xp1_phase1
#SBATCH --qos=large
#SBATCH -c 12
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
echo " 12 processus R en parallèle"
echo "=============================================="

# --- Environnement ---
source ~/.bashrc
load_spack
spack load r@4.4.0 
export R_LIBS=/scratch/${USER}/R

# Limiter BLAS à 1 thread par processus (éviter la surcharge)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd /scratch/${USER}/MagmaClustR

# --- Lancer les 12 combinaisons en parallèle ---
LOGDIR="/scratch/${USER}/logs/phase1_details"
mkdir -p "${LOGDIR}"

PIDS=()

# Créer un dossier de bibliothèque propre pour ce run
LIB_TEMP="/scratch/${USER}/R_temp_$(date +%Y%m%d)"
mkdir -p "${LIB_TEMP}"

# Installer le package cloné dans ce dossier spécifique
Rscript -e "install.packages('/scratch/${USER}/MagmaClustR', repos = NULL, type = 'source', lib = '${LIB_TEMP}')"

# Ajouter ce dossier AU DÉBUT du chemin de recherche de R
export R_LIBS="${LIB_TEMP}:${R_LIBS}"

for N_OUT in 2 3 4 6 8; do
  for N_TRAIN in 15 30 300; do
    LOGFILE="${LOGDIR}/momt_nout${N_OUT}_ntrain${N_TRAIN}.log"
    echo "[$(date +%H:%M:%S)] Lancement MOMT n_out=${N_OUT} n_train=${N_TRAIN} → ${LOGFILE}"

    /opt/spack/opt/spack/linux-debian11-zen2/gcc-13.2.0/r-4.4.0-tohpugilej6myswwe73dlbkypu7qqn4p/bin/Rscript --vanilla Benchmark_XP_1_MOMT_cluster.R \
      --n_out=${N_OUT} --n_train=${N_TRAIN} \
      > "${LOGFILE}" 2>&1 &

    PIDS+=($!)
  done
done

echo ""
echo "12 processus lancés. PIDs : ${PIDS[*]}"
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
echo " Échecs: ${FAILED} / 12"
echo "=============================================="

if [ ${FAILED} -gt 0 ]; then
  echo "ATTENTION : ${FAILED} processus ont échoué. Vérifiez les logs dans ${LOGDIR}/"
  exit 1
fi

exit 0
