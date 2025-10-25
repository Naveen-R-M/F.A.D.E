"""
Show the actual DiffSBDD command that will be run.
"""

from pathlib import Path
from fade.config import config

def show_diffsbdd_command():
    """Show the DiffSBDD command structure."""
    
    # Example values
    job_id = "test123"
    pdb_file = "/path/to/protein.pdb"
    pocket_residues = ["A:12", "A:61", "A:116"]
    n_samples = 100
    sanitize = True
    
    module_path = "/projects/SimBioSys/share/software/modulefiles"
    remote_dir = f"/scratch/rajagopalmohanraj.n/F.A.D.E/diffsbdd_inputs/{job_id}"
    remote_pdb = f"{remote_dir}/inputs/protein.pdb"
    output_file = f"generated_ligands_{job_id}.sdf"
    resi_list_str = " ".join(pocket_residues)
    
    command = f"""
cd {remote_dir}
module use {module_path}
module load diffsbdd

sbatch_diffsbdd_generate \\
    {job_id} \\
    checkpoints/crossdocked_fullatom_cond.ckpt \\
    --pdbfile {remote_pdb} \\
    --resi_list {resi_list_str} \\
    --outfile outputs/{output_file} \\
    --n_samples {n_samples} \\
    {"--sanitize" if sanitize else ""}
"""
    
    print("=== DiffSBDD Command Structure ===\n")
    print(command)
    print("\n=== Notes ===")
    print("- sbatch_diffsbdd_generate is a wrapper script that should exist after loading the module")
    print("- The checkpoint path is relative, which might cause issues")
    print("- The --resi_list format should match what DiffSBDD expects")
    print("\n=== Potential Issues ===")
    print("1. Checkpoint file path might need to be absolute")
    print("2. The sbatch_diffsbdd_generate command might not exist")
    print("3. The command syntax might be different")
    print("4. Module might not be setting up environment correctly")
    
    print("\n=== Alternative Command Format ===")
    alt_command = f"""
cd {remote_dir}
module use {module_path}
module load diffsbdd

# Direct python call (if sbatch wrapper doesn't exist)
python /path/to/diffsbdd/generate.py \\
    --pdb {remote_pdb} \\
    --residues {resi_list_str} \\
    --output outputs/{output_file} \\
    --n_samples {n_samples}
"""
    print(alt_command)


if __name__ == "__main__":
    show_diffsbdd_command()
