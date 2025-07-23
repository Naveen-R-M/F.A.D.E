#!/bin/bash
#SBATCH --partition=${PARTITION}
#SBATCH --constraint=${CONSTRAINT}
#SBATCH --gres=${GRES}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEMORY}
#SBATCH --time=${WALLTIME}
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${LOG_DIR}/${JOB_NAME}_%j.out
#SBATCH --error=${LOG_DIR}/${JOB_NAME}_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER}@northeastern.edu

# Load modules
module load ${MODULES}

# Set environment variables
export PYTHONPATH=${PROJECT_ROOT}:${PYTHONPATH}
export SINGULARITY_BIND="${BIND_PATHS}"

# Run the command
${COMMAND}

# Exit with the command's exit code
exit $?
