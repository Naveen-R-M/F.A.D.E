#!/bin/bash
#SBATCH --job-name=af3_KRAS
#SBATCH --output=/scratch/rajagopalmohanraj.n/F.A.D.E/logs/%x_%j.out
#SBATCH --error=/scratch/rajagopalmohanraj.n/F.A.D.E/logs/%x_%j.err
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:h200:1
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# Load modules
module load miniconda3/24.11.1
module load schrodinger/2024-4

# Set up environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate fade

# Set AlphaFold3 paths - required for container execution
export DATA_PATH=/shared/container_repository/AlphaFold/database
export MODEL_PATH=/projects/SimBioSys/share/software/AF3/models

# Load configuration
CONFIG_FILE="/scratch/rajagopalmohanraj.n/F.A.D.E/data/inputs/configs/KRAS_alphafold.json"
PROTEIN_NAME="KRAS"

echo "Starting AlphaFold3 prediction for ${PROTEIN_NAME}"
echo "Using configuration file: ${CONFIG_FILE}"
echo "DATA_PATH: ${DATA_PATH}"
echo "MODEL_PATH: ${MODEL_PATH}"

# Parse configuration
SEQUENCE_FILE=$(cat ${CONFIG_FILE} | jq -r '.sequence_file')
OUTPUT_DIR=$(cat ${CONFIG_FILE} | jq -r '.output_dir')
MAX_TEMPLATE_DATE=$(cat ${CONFIG_FILE} | jq -r '.max_template_date')
USE_GPU=$(cat ${CONFIG_FILE} | jq -r '.use_gpu')
MODEL_PRESET=$(cat ${CONFIG_FILE} | jq -r '.model_preset')
DB_PRESET=$(cat ${CONFIG_FILE} | jq -r '.db_preset')
NUM_RECYCLE=$(cat ${CONFIG_FILE} | jq -r '.num_recycle')
ENABLE_AMBER_RELAX=$(cat ${CONFIG_FILE} | jq -r '.enable_amber_relax')

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Run AlphaFold3 using Singularity container with exported paths
echo "Running AlphaFold3..."

singularity exec --nv \
    --bind ${DATA_PATH}:/data \
    --bind ${MODEL_PATH}:/models \
    /shared/container_repository/AlphaFold/alphafold3.sif \
    python /app/alphafold/run_alphafold.py \
    --fasta_paths=${SEQUENCE_FILE} \
    --output_dir=${OUTPUT_DIR} \
    --data_dir=/data \
    --model_dir=/models \
    --max_template_date=${MAX_TEMPLATE_DATE} \
    --model_preset=${MODEL_PRESET} \
    --db_preset=${DB_PRESET} \
    --num_multimer_predictions_per_model=1 \
    --recycle_tol=0.0 \
    --num_recycle=${NUM_RECYCLE} \
    --use_gpu_relax=${ENABLE_AMBER_RELAX} \
    --models_to_relax=best

echo "AlphaFold3 prediction completed"

# Extract best model
BEST_MODEL=$(find ${OUTPUT_DIR} -name "*rank_001*.pdb" | head -n 1)
if [ ! -z "${BEST_MODEL}" ]; then
    BEST_MODEL_BASENAME=$(basename ${BEST_MODEL})
    cp ${BEST_MODEL} ${OUTPUT_DIR}/${PROTEIN_NAME}_best_model.pdb
    echo "Best model saved as: ${OUTPUT_DIR}/${PROTEIN_NAME}_best_model.pdb"
else
    echo "No best model found"
fi

echo "Job completed for ${PROTEIN_NAME}"
