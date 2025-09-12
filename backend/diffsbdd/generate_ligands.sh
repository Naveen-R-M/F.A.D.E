#!/bin/bash

#SBATCH --job-name=diffsbdd-lig-gen
#SBATCH --partition=gpu-short
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00

# Generate de novo small molecule ligand candidates
# via Euclidean diffusion provided by DiffSBDD

# NOTE: For running specifically on the Explorer cluster.
# the DiffSBDD repo and its dependencies would have to be somehow
# included so that it could be run from FADE internally (?).
# Here, I provided the full path to where I cloned the DiffSBDD repo,
# where I pulled the model checkpoint and saved it to,
# and where I installed its python dependencies to an environment called diffsbdd.

n_samples=50
ckpt="/scratch/kantorow.j/fade-tests/diffsbdd-tests/models/checkpoints/crossdocked_fullatom_cond.ckpt"

while [[ "$#" -gt 0 ]]; do
	case $1 in
		--pdbfile) pdbfile="$2"; shift ;;
		--ref_ligand) ref_ligand="$2"; shift ;;
		--outfile) outfile="$2"; shift ;;
		--n_samples) n_samples="$2"; shift ;;
		--checkpoint) ckpt="$2"; shift ;;
		*) echo "Unknown parameter passed: $1"; exit 1 ;;
	esac
	shift
done

genscript="/scratch/kantorow.j/fade-tests/diffsbdd-tests/DiffSBDD/generate_ligands.py"

source /projects/SimBioSys/jkant/miniconda3-jkant/miniconda3-jkant/etc/profile.d/conda.sh
conda activate diffsbdd

python $genscript $ckpt --pdbfile $pdbfile --ref_ligand $ref_ligand --outfile $outfile --n_samples $n_samples

exit 0