#!/bin/bash
# ==========================================================================
# setup_env.sh — Création de l'environnement conda pour TACTiS-2
#
# Usage (sur le cluster, après `load_spack` et `spack load miniconda3@24.7.1`) :
#   bash setup_env.sh
# ==========================================================================

set -e

# --- MODIFICATION ICI : On utilise un chemin (PATH) au lieu d'un nom (NAME) ---
ENV_PATH="/scratch/agrenoui/NeurIPS_experiments/Experience_2/TACTiS2/env_tactis2"

echo "=== Création de l'environnement conda dans '${ENV_PATH}' ==="

# --- MODIFICATION ICI : Suppression basée sur le dossier ---
conda deactivate 2>/dev/null || true
if [ -d "${ENV_PATH}" ]; then
    echo "Suppression de l'environnement existant..."
    conda env remove -p "${ENV_PATH}" -y
fi

# --- MODIFICATION ICI : Création avec -p au lieu de -n ---
conda create -p "${ENV_PATH}" python=3.10 -y

# --- MODIFICATION ICI : Activation avec le chemin complet ---
source activate "${ENV_PATH}"

echo "=== Installation des dépendances Python ==="

# PyTorch CPU (pas de GPU sur le cluster IMT)
pip install --no-cache-dir torch==2.5.1 --index-url https://download.pytorch.org/whl/cpu

# Dépendances scientifiques
pip install --no-cache-dir \
    numpy==2.2.5 \
    pandas==2.2.3 \
    scipy==1.15.2 \
    matplotlib==3.10.1 \
    scikit-learn==1.6.1

# pyreadr et rdata pour lire les fichiers .rds depuis Python
pip install --no-cache-dir pyreadr rdata

# TACTiS-2 : installation depuis le dépôt GitHub officiel
pip install --no-cache-dir git+https://github.com/ServiceNow/TACTiS.git

echo ""
echo "=== Vérification de l'installation ==="
python -c "
import torch; print(f'PyTorch {torch.__version__}')
import pandas; print(f'Pandas {pandas.__version__}')
import numpy; print(f'NumPy {numpy.__version__}')
import pyreadr; print(f'pyreadr {pyreadr.__version__}')
from tactis.model.tactis import TACTiS; print('TACTiS OK')
print('\nTout est installé correctement.')
"

echo ""
# --- MODIFICATION ICI : Message de fin mis à jour ---
echo "=== Environnement prêt ==="
echo "Pour l'utiliser : source activate ${ENV_PATH}"

# #!/bin/bash
# # ==========================================================================
# # setup_env.sh — Création de l'environnement conda pour TACTiS-2
# #
# # Usage (sur le cluster, après `load_spack` et `spack load miniconda3@24.7.1`) :
# #   bash setup_env.sh
# # ==========================================================================

# set -e

# ENV_NAME="tactis2"

# echo "=== Création de l'environnement conda '${ENV_NAME}' ==="

# # Supprimer l'ancien environnement s'il existe
# conda deactivate 2>/dev/null || true
# if conda env list | grep -q "^${ENV_NAME} "; then
#     echo "Suppression de l'environnement existant..."
#     conda env remove -n "${ENV_NAME}" -y
# fi

# # Créer l'environnement avec Python 3.10 (compatible avec le cluster)
# conda create -n "${ENV_NAME}" python=3.10 -y

# # Activer l'environnement
# source activate "${ENV_NAME}"

# echo "=== Installation des dépendances Python ==="

# # PyTorch CPU (pas de GPU sur le cluster IMT)
# pip install --no-cache-dir torch==2.5.1 --index-url https://download.pytorch.org/whl/cpu

# # Dépendances scientifiques
# pip install --no-cache-dir \
#     numpy==2.2.5 \
#     pandas==2.2.3 \
#     scipy==1.15.2 \
#     matplotlib==3.10.1 \
#     scikit-learn==1.6.1

# # pyreadr et rdata pour lire les fichiers .rds depuis Python
# pip install --no-cache-dir pyreadr rdata

# # TACTiS-2 : installation depuis le dépôt GitHub officiel
# pip install --no-cache-dir git+https://github.com/ServiceNow/TACTiS.git

# echo ""
# echo "=== Vérification de l'installation ==="
# python -c "
# import torch; print(f'PyTorch {torch.__version__}')
# import pandas; print(f'Pandas {pandas.__version__}')
# import numpy; print(f'NumPy {numpy.__version__}')
# import pyreadr; print(f'pyreadr {pyreadr.__version__}')
# from tactis.model.tactis import TACTiS; print('TACTiS OK')
# print('\\nTout est installé correctement.')
# "

# echo ""
# echo "=== Environnement '${ENV_NAME}' prêt ==="
# echo "Pour l'utiliser : source activate ${ENV_NAME}"
