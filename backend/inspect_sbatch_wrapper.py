"""
Inspect the sbatch_diffsbdd_generate wrapper script to understand
what it does and whether it activates the environment.
"""

from fade.tools.ssh_client import get_hpc_client


def inspect_sbatch_wrapper():
    """Read and display the sbatch_diffsbdd_generate script."""
    
    ssh = get_hpc_client()
    if not ssh.connect():
        print("❌ Failed to connect")
        return
    
    module_path = "/projects/SimBioSys/share/software/modulefiles"
    
    print("="*60)
    print("Inspecting sbatch_diffsbdd_generate")
    print("="*60 + "\n")
    
    # Find the script
    cmd = f'''
module use {module_path}
module load diffsbdd
which sbatch_diffsbdd_generate
'''
    
    stdout, stderr, code = ssh.execute_command(cmd)
    
    if code != 0 or not stdout.strip():
        print("❌ Script not found")
        return
    
    script_path = stdout.strip()
    print(f"Script location: {script_path}\n")
    
    # Read the script
    stdout, stderr, code = ssh.execute_command(f"cat {script_path}")
    
    if code == 0:
        print("Script contents:")
        print("-"*60)
        print(stdout)
        print("-"*60)
        
        # Analyze the script
        print("\n" + "="*60)
        print("Analysis:")
        print("="*60)
        
        if "micromamba activate" in stdout or "conda activate" in stdout:
            print("✓ Script includes environment activation")
        else:
            print("⚠ Script does NOT include environment activation")
            print("  This might be the issue - the job needs environment activation")
        
        if "generate_ligands.py" in stdout:
            print("✓ Script calls generate_ligands.py")
        
        if "#SBATCH" in stdout:
            print("✓ Script is a SLURM batch script")
            
            # Extract SBATCH directives
            print("\n  SLURM directives:")
            for line in stdout.split('\n'):
                if line.strip().startswith("#SBATCH"):
                    print(f"    {line.strip()}")
    else:
        print(f"❌ Failed to read script: {stderr}")
    
    ssh.close()


if __name__ == "__main__":
    inspect_sbatch_wrapper()
