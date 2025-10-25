"""
Quick diagnostic runner for DiffSBDD debugging.
Run this after the fix to see what's actually happening.
"""

import sys
from fade.tools.ssh_client import get_hpc_client
from fade.utils import get_logger

logger = get_logger("diagnostic")


def run_all_diagnostics(slurm_job_id: str = None):
    """Run all diagnostic checks."""
    
    print("\n" + "="*60)
    print("F.A.D.E DiffSBDD Diagnostics")
    print("="*60 + "\n")
    
    ssh = get_hpc_client()
    if not ssh.connect():
        print("❌ Failed to connect to HPC")
        return
    
    print("✓ Connected to HPC\n")
    
    # Check 1: Module system
    print("1. Checking module system...")
    stdout, _, code = ssh.execute_command("module list 2>&1")
    if code == 0:
        print(f"✓ Module system working")
        print(f"   Currently loaded: {stdout.strip() if stdout.strip() else 'none'}")
    else:
        print("❌ Module system not working")
    
    # Check 2: DiffSBDD module
    print("\n2. Checking DiffSBDD module...")
    module_path = "/projects/SimBioSys/share/software/modulefiles"
    stdout, _, code = ssh.execute_command(
        f"module use {module_path} && module load diffsbdd && echo 'SUCCESS' 2>&1"
    )
    if "SUCCESS" in stdout:
        print("✓ DiffSBDD module loads successfully")
    else:
        print("❌ DiffSBDD module failed to load")
        print(f"   Output: {stdout}")
    
    # Check 3: sbatch command
    print("\n3. Checking for sbatch_diffsbdd_generate...")
    stdout, _, code = ssh.execute_command(
        f"module use {module_path} && module load diffsbdd && which sbatch_diffsbdd_generate 2>&1"
    )
    if code == 0 and stdout.strip():
        print(f"✓ Command found: {stdout.strip()}")
    else:
        print("❌ Command not found")
        
        # Try to find it
        print("   Searching for similar commands...")
        stdout, _, _ = ssh.execute_command(
            f"module use {module_path} && module load diffsbdd && compgen -c | grep -i diff"
        )
        if stdout:
            print(f"   Found: {stdout}")
    
    # Check 4: Checkpoint file
    print("\n4. Checking for checkpoint file...")
    stdout, _, code = ssh.execute_command(
        "find /projects/SimBioSys -name '*crossdocked*.ckpt' -type f 2>/dev/null | head -1"
    )
    if stdout.strip():
        print(f"✓ Checkpoint file found: {stdout.strip()}")
    else:
        print("❌ Checkpoint file not found")
    
    # Check 5: Recent jobs (if job ID provided)
    if slurm_job_id:
        print(f"\n5. Checking SLURM job {slurm_job_id}...")
        
        # Check squeue
        stdout, _, code = ssh.execute_command(f"squeue -j {slurm_job_id}")
        if code == 0 and len(stdout.split('\n')) > 1:
            print(f"✓ Job is in queue:\n{stdout}")
        else:
            print("   Job not in squeue, checking sacct...")
            
            # Check sacct
            stdout, _, code = ssh.execute_command(
                f"sacct -j {slurm_job_id} --format=JobID,State,ExitCode,Start,End"
            )
            if stdout.strip():
                print(f"   Job history:\n{stdout}")
            
            # Check for output file
            print(f"\n   Looking for output files...")
            stdout, _, _ = ssh.execute_command(
                f"find /scratch/rajagopalmohanraj.n/F.A.D.E -name 'slurm-{slurm_job_id}.*' 2>/dev/null"
            )
            if stdout.strip():
                files = stdout.strip().split('\n')
                for file in files:
                    print(f"   Found: {file}")
                    # Print last 20 lines
                    out, _, _ = ssh.execute_command(f"tail -20 {file}")
                    print(f"   Last 20 lines:\n{out}")
            else:
                print(f"   No output files found")
    
    # Check 6: Workspace
    print("\n6. Checking workspace...")
    stdout, _, _ = ssh.execute_command(
        "ls -lh /scratch/rajagopalmohanraj.n/F.A.D.E/ 2>&1"
    )
    print(f"   Workspace contents:\n{stdout}")
    
    ssh.close()
    
    print("\n" + "="*60)
    print("Diagnostics Complete")
    print("="*60 + "\n")
    
    print("Next steps:")
    print("1. If DiffSBDD module doesn't load: Contact HPC admin")
    print("2. If sbatch command not found: Need to find correct command name")
    print("3. If checkpoint not found: Need to get checkpoint path from admin")
    print("4. If job failed: Read the error logs above")


if __name__ == "__main__":
    job_id = sys.argv[1] if len(sys.argv) > 1 else None
    
    if job_id:
        print(f"Running diagnostics for job {job_id}...")
    else:
        print("Running general diagnostics...")
        print("(Provide SLURM job ID as argument to check specific job)")
    
    run_all_diagnostics(job_id)
