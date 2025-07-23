#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --constraint=gpu:v100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --job-name=alphafold3
#SBATCH --output=${LOG_DIR}/alphafold3_%j.out
#SBATCH --error=${LOG_DIR}/alphafold3_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER}@northeastern.edu

# Load modules
module load miniconda3/24.11.1
module load cuda/12.3.0

# Set environment variables
export PYTHONPATH=${PROJECT_ROOT}:${PYTHONPATH}

# Run AlphaFold3
singularity exec --nv ${ALPHAFOLD3_CONTAINER} \
  python /app/alphafold/run_alphafold.py \
  --fasta_paths=${FASTA_PATH} \
  --output_dir=${OUTPUT_DIR} \
  --model_preset=monomer \
  --db_preset=full_dbs \
  --max_template_date=2024-07-15 \
  --use_gpu_relax=True

# Exit with the command's exit code
exit $?
