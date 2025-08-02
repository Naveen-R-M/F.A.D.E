#!/bin/bash
#SBATCH --job-name=af3_ACE2
#SBATCH --output=/scratch/rajagopalmohanraj.n/F.A.D.E/logs/%x_%j.out
#SBATCH --error=/scratch/rajagopalmohanraj.n/F.A.D.E/logs/%x_%j.err
#SBATCH --partition=gpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:h200:1
#SBATCH --mem=128G
#SBATCH --time=02:00:00

set -euo pipefail
# module load singularity/3.5.3

#— AlphaFold3 container & resource dirs
export AF3_RESOURCES_DIR=/shared/container_repository/AlphaFold
export AF3_IMAGE=${AF3_RESOURCES_DIR}/alphafold3.sif
export DATA_PATH=/shared/container_repository/AlphaFold/database
export MODEL_PATH=/projects/SimBioSys/share/software/AF3/models
export TEMPLATE_DATE=2025-07-15


# Force the correct JAX env vars inside the container.
# Use both prefixes to satisfy Singularity & Apptainer:
export SINGULARITYENV_XLA_CLIENT_MEM_FRACTION=0.95
export APPTAINERENV_XLA_CLIENT_MEM_FRACTION=0.95
export SINGULARITYENV_XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
export APPTAINERENV_XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

#— Load your JSON config & name
CONFIG_FILE="data/inputs/configs/ACE2_alphafold.json"
PROTEIN_NAME="ACE2"

echo "Starting AlphaFold3 for ${PROTEIN_NAME}"
echo "  Config   : ${CONFIG_FILE}"

#— Parse core settings via jq
SEQUENCE_FILE=$(jq -r '.sequence_file'    ${CONFIG_FILE})
MODEL_SEEDS=$(jq -r '.model_seeds | join(",")' ${CONFIG_FILE})
OUTPUT_DIR=$(jq -r '.output_dir'           ${CONFIG_FILE})
MAX_TEMPLATE_DATE=$(jq -r '.max_template_date' ${CONFIG_FILE})
MODEL_PRESET=$(jq -r '.model_preset'       ${CONFIG_FILE})
DB_PRESET=$(jq -r '.db_preset'             ${CONFIG_FILE})
NUM_RECYCLE=$(jq -r '.num_recycle'         ${CONFIG_FILE})
ENABLE_AMBER_RELAX=$(jq -r '.enable_amber_relax' ${CONFIG_FILE})

LIGAND_ID=$(jq -r '.ligand.id // empty' ${CONFIG_FILE})
if [ -n "${LIGAND_ID}" ]; then
    LIGAND_CCD_CODES=$(jq -r '.ligand.ccd_codes // [] | map("\""+.+"\"") | join(",")' ${CONFIG_FILE})
else
    LIGAND_CCD_CODES=""
fi

#— Prep output dir
mkdir -p ${OUTPUT_DIR}

#— Read FASTA
FASTA_CONTENT=$(cat ${SEQUENCE_FILE})
SEQUENCE=$(echo "${FASTA_CONTENT}" | grep -v "^>" | tr -d '\n')
SEQUENCE_ID=$(echo "${FASTA_CONTENT}" | grep "^>" | sed 's/^>//' | head -n 1)

#— Build AlphaFold3 JSON
AF3_INPUT_JSON="${OUTPUT_DIR}/${PROTEIN_NAME}_af3_input.json"
{
  echo '{'
  echo '  "name": "'"${PROTEIN_NAME}"'",'
  echo '  "modelSeeds": ['"${MODEL_SEEDS}"'],'
  echo '  "sequences": ['
  echo '    {'
  echo '      "protein": {'
  echo '        "id": "'"${SEQUENCE_ID}"'",'
  echo '        "sequence": "'"${SEQUENCE}"'"'
  echo '      }'
  echo '    }'
  if [ -n "${LIGAND_ID}" ]; then
    echo '    ,{'
    echo '      "ligand": {'
    echo '        "id": "'"${LIGAND_ID}"'",'
    echo '        "ccdCodes": ['"${LIGAND_CCD_CODES}"']'
    echo '      }'
    echo '    }'
  fi
  echo '  ],'
  echo '  "dialect": "alphafold3",'
  echo '  "version": 1'
  echo '}'
} > "${AF3_INPUT_JSON}"

echo "Generated JSON: ${AF3_INPUT_JSON}"
cat ${AF3_INPUT_JSON}

#— Run the container
echo "Running AlphaFold3..."

singularity exec --nv --cleanenv \
  --bind /projects:/work \
  --bind /scratch:/scratch \
  --bind /shared:/shared \
  --bind "${DATA_PATH}":/data \
  --bind "${MODEL_PATH}":/models \
  --bind /lib64/libcuda.so.1:/usr/local/nvidia/lib64/libcuda.so.1 \
  --bind /lib64/libnvidia-ml.so.1:/usr/local/nvidia/lib64/libnvidia-ml.so.1 \
  ${AF3_IMAGE} \
  python /app/alphafold/run_alphafold.py \
    --json_path="${AF3_INPUT_JSON}" \
    --output_dir="${OUTPUT_DIR}/$(basename "${AF3_INPUT_JSON}" .json)" \
    --logtostderr \
    --db_dir       "/data" \
    --model_dir    "/models" \
    --max_template_date "${MAX_TEMPLATE_DATE}" \
    --model_preset "${MODEL_PRESET}" \
    --db_preset "${DB_PRESET}" \
    --num_recycle "${NUM_RECYCLE}" \
    --enable_amber_relax="${ENABLE_AMBER_RELAX}" \
    --flash_attention_implementation xla

echo "AlphaFold3 prediction completed"

#— Extract & copy best model
BEST_MODEL=$(find ${OUTPUT_DIR} -name "*rank_001*.pdb" | head -n 1 \
             || find ${OUTPUT_DIR} -name "*model_0*.pdb" | head -n 1)
if [ -n "${BEST_MODEL}" ]; then
  cp "${BEST_MODEL}" "${OUTPUT_DIR}/${PROTEIN_NAME}_best_model.pdb"
  echo "Best model → ${OUTPUT_DIR}/${PROTEIN_NAME}_best_model.pdb"
else
  echo "⚠️ No best model found; listing contents:"
  ls -la "${OUTPUT_DIR}"
fi

echo "Job done for ${PROTEIN_NAME}"
