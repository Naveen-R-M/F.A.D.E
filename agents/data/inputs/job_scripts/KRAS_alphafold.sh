#!/bin/bash
#SBATCH --job-name=af3_KRAS
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

# ── AlphaFold 3 container and resource paths ─────────────────────────────
export AF3_RESOURCES_DIR=/shared/container_repository/AlphaFold
export AF3_IMAGE=${AF3_RESOURCES_DIR}/alphafold3.sif
export DATA_PATH=/shared/container_repository/AlphaFold/database
export MODEL_PATH=/projects/SimBioSys/share/software/AF3/models
export TEMPLATE_DATE=2025-07-15        # fallback if config omits it

# Ensure JAX inside container respects GPU memory fraction
export SINGULARITYENV_XLA_CLIENT_MEM_FRACTION=0.95
export APPTAINERENV_XLA_CLIENT_MEM_FRACTION=0.95
export SINGULARITYENV_XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
export APPTAINERENV_XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

# ── Configuration file supplied by the user ──────────────────────────────
CONFIG_FILE="/scratch/rajagopalmohanraj.n/F.A.D.E/agents/data/inputs/configs/KRAS_alphafold.json"
PROTEIN_NAME="KRAS"

echo "Starting AlphaFold 3 for ${PROTEIN_NAME}"
echo "  Config: ${CONFIG_FILE}"

# jq is required; abort with a clear message if missing
command -v jq >/dev/null 2>&1 || { echo >&2 "❌ jq not found — install or module-load it."; exit 1; }

# ── Pull fields we still use from the JSON config ────────────────────────
SEQUENCE_FILE=$(jq -r '.sequence_file'             "${CONFIG_FILE}")
MODEL_SEEDS=$(jq -r '.model_seeds | join(",")'      "${CONFIG_FILE}")
OUTPUT_DIR=$(jq -r '.output_dir'                    "${CONFIG_FILE}")
MAX_TEMPLATE_DATE=$(jq -r '.max_template_date // ""' "${CONFIG_FILE}")
NUM_RECYCLE=$(jq -r '.num_recycle // 10'            "${CONFIG_FILE}")  # default 10 cycles

# Optional ligand section
LIGAND_ID=$(jq -r '.ligand.id // empty' "${CONFIG_FILE}")
if [ -n "${LIGAND_ID}" ]; then
  LIGAND_CCD_CODES=$(jq -r \
    '.ligand.ccd_codes // [] | map("\""+.+"\"") | join(",")' "${CONFIG_FILE}")
fi

mkdir -p "${OUTPUT_DIR}"

# ── Build a valid AlphaFold 3 input JSON ─────────────────────────────────
FASTA_CONTENT=$(cat "${SEQUENCE_FILE}")
SEQUENCE=$(echo "${FASTA_CONTENT}" | grep -v "^>" | tr -d '\n')

# Clean header → uppercase letters only; fall back to PROTEIN_NAME
SEQ_ID_RAW=$(echo "${FASTA_CONTENT}" | grep "^>" | head -n1 | sed 's/^>//')
SEQUENCE_ID=$(echo "${SEQ_ID_RAW}" | tr '[:lower:]' '[:upper:]' | tr -cd 'A-Z')
if [ -z "${SEQUENCE_ID}" ]; then
  SEQUENCE_ID=$(echo "${PROTEIN_NAME}" | tr '[:lower:]' '[:upper:]' | tr -cd 'A-Z')
fi

AF3_INPUT_JSON="${OUTPUT_DIR}/${PROTEIN_NAME}_af3_input.json"
{
  echo '{'
  echo '  "name": "'"${PROTEIN_NAME}"'",'
  echo '  "modelSeeds": ['"${MODEL_SEEDS}"'],'
  echo '  "sequences": ['
  echo '    { "protein": { "id": "'"${SEQUENCE_ID}"'", "sequence": "'"${SEQUENCE}"'" } }'
  if [ -n "${LIGAND_ID:-}" ]; then
    echo '    ,{ "ligand": { "id": "'"${LIGAND_ID}"'", "ccdCodes": ['"${LIGAND_CCD_CODES}"'] } }'
  fi
  echo '  ],'
  echo '  "dialect": "alphafold3",'
  echo '  "version": 1'
  echo '}'
} > "${AF3_INPUT_JSON}"

echo "Generated AF3 JSON:"
cat "${AF3_INPUT_JSON}"

# ── Launch AlphaFold 3 inside the container ──────────────────────────────
singularity exec --nv --cleanenv \
  --bind /projects:/work \
  --bind /scratch:/scratch \
  --bind /shared:/shared \
  --bind "${DATA_PATH}":/data \
  --bind "${MODEL_PATH}":/models \
  --bind /lib64/libcuda.so.1:/usr/local/nvidia/lib64/libcuda.so.1 \
  --bind /lib64/libnvidia-ml.so.1:/usr/local/nvidia/lib64/libnvidia-ml.so.1 \
  "${AF3_IMAGE}" \
  python /app/alphafold/run_alphafold.py \
    --json_path        "${AF3_INPUT_JSON}" \
    --output_dir       "${OUTPUT_DIR}/$(basename "${AF3_INPUT_JSON}" .json)" \
    --logtostderr \
    --db_dir           "/data" \
    --model_dir        "/models" \
    --max_template_date "${MAX_TEMPLATE_DATE:-${TEMPLATE_DATE}}" \
    --num_recycles     "${NUM_RECYCLE}" \
    --flash_attention_implementation xla

echo "AlphaFold 3 prediction completed."

# ── Copy the best-ranked model to a predictable filename ─────────────────
BEST_MODEL=$(find "${OUTPUT_DIR}" -name "*rank_001*.pdb" | head -n1 \
             || find "${OUTPUT_DIR}" -name "*model_0*.pdb" | head -n1)
if [ -n "${BEST_MODEL}" ]; then
  cp "${BEST_MODEL}" "${OUTPUT_DIR}/${PROTEIN_NAME}_best_model.pdb"
  echo "Best model saved → ${OUTPUT_DIR}/${PROTEIN_NAME}_best_model.pdb"
else
  echo "⚠️  No best model found; directory listing:"
  ls -la "${OUTPUT_DIR}"
fi

echo "Job finished for ${PROTEIN_NAME}"
