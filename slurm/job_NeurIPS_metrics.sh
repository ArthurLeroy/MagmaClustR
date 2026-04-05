#!/bin/bash -l
#===============================================================================
# job_NeurIPS_metrics.sh
# NeurIPS : Calcul des métriques après toutes les expériences
#
# Soumission : sbatch job_NeurIPS_metrics.sh
#===============================================================================

#SBATCH --job-name=neurips_metrics
#SBATCH --qos=huge
#SBATCH -c 1
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/%u/logs/neurips_metrics_%j.out
#SBATCH --error=/scratch/%u/logs/neurips_metrics_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alexia.grenouillat@math.univ-toulouse.fr

echo "=============================================="
echo " NeurIPS : Calcul des métriques"
echo " Date    : $(date)"
echo " Noeud   : $(hostname)"
echo "=============================================="

# --- Environnement ---
source ~/.bashrc
load_spack
spack load r@4.4.0
export R_LIBS=/scratch/${USER}/R

LIB_TEMP="/scratch/${USER}/R_temp_$(date +%Y%m%d)"
if [ -d "${LIB_TEMP}" ]; then
  export R_LIBS="${LIB_TEMP}:${R_LIBS}"
fi

cd /scratch/${USER}/MagmaClustR

RSCRIPT="/opt/spack/opt/spack/linux-debian11-zen2/gcc-13.2.0/r-4.4.0-tohpugilej6myswwe73dlbkypu7qqn4p/bin/Rscript"
SCRIPT_DIR="/scratch/${USER}/NeurIPS_experiments/scripts"

${RSCRIPT} --vanilla "${SCRIPT_DIR}/compute_metrics_NeurIPS.R"

echo ""
echo "=============================================="
echo " Métriques terminées : $(date)"
echo "=============================================="
