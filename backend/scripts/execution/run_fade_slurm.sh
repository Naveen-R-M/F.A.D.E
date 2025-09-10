#!/bin/bash
#SBATCH --job-name=FADE_pipeline
#SBATCH --output=fade_%j.out
#SBATCH --error=fade_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=gpu-short

# Load required modules
module load miniconda3/24.11.1
module load schrodinger/2024-4

# Activate conda environment
source activate $SCRATCH/conda-envs/fade

# Set working directory
cd $SLURM_SUBMIT_DIR

# Run F.A.D.E pipeline
python main.py --query "$1" --output-dir "results_$SLURM_JOB_ID" --log-level INFO

# Optional: Send email notification when job completes
mail -s "F.A.D.E Job $SLURM_JOB_ID Completed" $USER@northeastern.edu << EOF
Your F.A.D.E job has completed.
Query: $1
Job ID: $SLURM_JOB_ID
Results directory: results_$SLURM_JOB_ID
EOF
