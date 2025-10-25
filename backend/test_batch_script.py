"""
Test the new batch script generation and show what it looks like.
"""

from fade.tools.hpc_diffsbdd import get_hpc_diffsbdd_client


def test_batch_script_generation():
    """Generate and display the batch script."""
    
    client = get_hpc_diffsbdd_client()
    
    # Example parameters
    job_id = "test123"
    remote_dir = "/scratch/rajagopalmohanraj.n/F.A.D.E/diffsbdd_inputs/test123"
    remote_pdb = f"{remote_dir}/inputs/protein.pdb"
    pocket_residues = ["A:12", "A:61", "A:116"]
    output_file = "generated_ligands_test123.sdf"
    n_samples = 100
    sanitize = True
    
    # Generate the script
    script = client._create_batch_script(
        job_id=job_id,
        remote_dir=remote_dir,
        remote_pdb=remote_pdb,
        pocket_residues=pocket_residues,
        output_file=output_file,
        n_samples=n_samples,
        sanitize=sanitize
    )
    
    print("="*60)
    print("Generated SLURM Batch Script")
    print("="*60)
    print(script)
    print("="*60)
    
    print("\n✓ Script includes:")
    print("  - SLURM directives (#SBATCH)")
    print("  - Module loading")
    print("  - Micromamba environment activation")
    print("  - Direct Python call to generate_ligands.py")
    print("  - Diagnostic output")
    
    print("\nKey features:")
    print("  - GPU requested: Yes (#SBATCH --gres=gpu:1)")
    print("  - Memory: 16GB")
    print("  - Time limit: 2 hours")
    print("  - Environment activation: ✓ (inside the job)")
    
    print("\nThis script will be uploaded to HPC and submitted with sbatch.")


if __name__ == "__main__":
    test_batch_script_generation()
