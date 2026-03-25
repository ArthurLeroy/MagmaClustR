#!/bin/bash
# ==========================================================================
# submit_tactis2.sh — Script de soumission SLURM pour le benchmark TACTiS-2
#
# Usage :
#   sbatch submit_tactis2.sh
# ==========================================================================

#SBATCH --job-name=tactis2_bench
#SBATCH --partition=fourier
#SBATCH --qos=huge
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/agrenoui/NeurIPS_experiments/Experience_2/logs_TACTiS/slurm_%j.out
#SBATCH --error=/scratch/agrenoui/NeurIPS_experiments/Experience_2/logs_TACTiS/slurm_%j.err

# --- Créer les dossiers de sortie ---
mkdir -p /scratch/agrenoui/NeurIPS_experiments/Experience_2/Models_TACTiS
mkdir -p /scratch/agrenoui/NeurIPS_experiments/Experience_2/Predictions_TACTiS
mkdir -p /scratch/agrenoui/NeurIPS_experiments/Experience_2/logs_TACTiS

echo "=== Job SLURM: ${SLURM_JOB_ID} ==="
echo "  Node     : $(hostname)"
echo "  CPUs     : ${SLURM_CPUS_PER_TASK}"
echo "  Début    : $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# --- Chargement de l'environnement ---
load_spack
spack load miniconda3@24.7.1%gcc@13.2.0
# source activate tactis2
source activate /scratch/agrenoui/NeurIPS_experiments/Experience_2/TACTiS2/env_tactis2

# --- Vérification rapide ---
echo "Python : $(which python)"
python -c "import torch; print(f'PyTorch {torch.__version__}')"
python -c "from tactis.model.tactis import TACTiS; print('TACTiS OK')"
echo ""

# --- Lancement du benchmark ---
# SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# cd "${SCRIPT_DIR}"
cd $SLURM_SUBMIT_DIR

# Empêcher la surcharge des CPU en forçant 1 thread par worker
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

python -u run_benchmark.py

echo ""
echo "=== Fin : $(date '+%Y-%m-%d %H:%M:%S') ==="
