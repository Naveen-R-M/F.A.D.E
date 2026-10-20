#!/bin/bash
#SBATCH --job-name=egfr-atp-check
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/rajagopalmohanraj.n/fade-egfr-validation/atp-site/validation-%j.log
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
unset PYTHONPATH PYTHONHOME || true
PY=/scratch/rajagopalmohanraj.n/conda-envs/fade-app310/bin/python
$PY -c 'import Bio,numpy; print("Biopython",Bio.__version__,"numpy",numpy.__version__)'
$PY /scratch/rajagopalmohanraj.n/fade-egfr-validation/atp-site/validate.py
