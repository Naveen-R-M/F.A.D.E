#!/bin/bash
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --job-name=docking
#SBATCH --output=${LOG_DIR}/docking_%j.out
#SBATCH --error=${LOG_DIR}/docking_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER}@northeastern.edu

# Load modules
module load schrodinger/2024-4

# Set environment variables
export PYTHONPATH=${PROJECT_ROOT}:${PYTHONPATH}

# Run Glide docking
$SCHRODINGER/glide \
  -JOBNAME ${JOB_NAME} \
  -HOST localhost \
  -WAIT \
  -LOCAL \
  -NOJOBID \
  ${INPUT_FILE}

# Exit with the command's exit code
exit $?
