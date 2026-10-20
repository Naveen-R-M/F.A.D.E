#!/bin/bash
#SBATCH --job-name=fade-egfr
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/rajagopalmohanraj.n/fade-egfr-validation/run-%j.log
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
unset PYTHONPATH PYTHONHOME FADE_PRIVATE_ENV || true
export FADE_LOCAL_HPC=1 HPC_USER=rajagopalmohanraj.n
export DIFFSBDD_PYTHON=/scratch/rajagopalmohanraj.n/conda-envs/fade-diffsbdd-venv/bin/python
export BOLTZ2_ENV_SETUP='export PATH=/scratch/rajagopalmohanraj.n/conda-envs/fade-boltz2/bin:$PATH'
export BOLTZ_CACHE_DIR=/scratch/rajagopalmohanraj.n/boltz_cache
export MAX_MOLECULES_TO_GENERATE=5 MAX_MOLECULES_TO_SCREEN=2 BOLTZ2_AFFINITY_MAX_LIGANDS=2
PY=/scratch/rajagopalmohanraj.n/conda-envs/fade-app310/bin/python
cd /scratch/rajagopalmohanraj.n/F.A.D.E
$PY -m pip check
$PY -m pytest backend/test_workflow_repairs.py -q
# Check private overlays after the interrupted shared repair, before launching jobs.
$DIFFSBDD_PYTHON -c 'import torch,torch_scatter,triton; from rdkit import Chem; print(torch.__version__,torch_scatter.__version__,triton.__version__)'
/scratch/rajagopalmohanraj.n/conda-envs/fade-boltz2/bin/python -c 'from boltz.main import cli; from rdkit import Chem; import torch; print(torch.__version__)'
cd backend
$PY -u /scratch/rajagopalmohanraj.n/fade-egfr-validation/run.py "$1"
