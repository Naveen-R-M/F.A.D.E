"""
Quick script to check SLURM job status and logs.
Usage: python check_job_status.py <slurm_job_id>
"""

import sys
from fade.tools.ssh_client import get_hpc_client

def check_job_status(slurm_job_id: str):
    """Check status of a SLURM job."""
    ssh_client = get_hpc_client()
    
    if not ssh_client.connect():
        print("Failed to connect to HPC")
        return
    
    print(f"\n=== Checking SLURM Job {slurm_job_id} ===\n")
    
    # Check squeue
    print("1. Checking squeue (running jobs)...")
    stdout, stderr, exit_code = ssh_client.execute_command(
        f"squeue -j {slurm_job_id} -l"
    )
    print(f"Exit code: {exit_code}")
    if stdout:
        print(f"Output:\n{stdout}")
    if stderr:
        print(f"Error:\n{stderr}")
    
    # Check sacct
    print("\n2. Checking sacct (job history)...")
    stdout, stderr, exit_code = ssh_client.execute_command(
        f"sacct -j {slurm_job_id} --format=JobID,JobName,State,ExitCode,Start,End,Elapsed -P"
    )
    print(f"Exit code: {exit_code}")
    if stdout:
        print(f"Output:\n{stdout}")
    if stderr:
        print(f"Error:\n{stderr}")
    
    # Find job directory
    print("\n3. Looking for job output files...")
    stdout, stderr, exit_code = ssh_client.execute_command(
        f"find /scratch/rajagopalmohanraj.n/F.A.D.E -name 'slurm-{slurm_job_id}.out' 2>/dev/null"
    )
    
    if stdout.strip():
        output_file = stdout.strip().split('\n')[0]
        print(f"Found output file: {output_file}")
        
        # Print last 50 lines
        print("\n4. SLURM output (last 50 lines):")
        stdout, stderr, exit_code = ssh_client.execute_command(
            f"tail -50 {output_file}"
        )
        if stdout:
            print(stdout)
    else:
        print("No output file found yet (job may not have started)")
        
        # Check in diffsbdd_inputs directories
        print("\n5. Checking diffsbdd_inputs directories...")
        stdout, stderr, exit_code = ssh_client.execute_command(
            f"ls -ltr /scratch/rajagopalmohanraj.n/F.A.D.E/diffsbdd_inputs/*/slurm-{slurm_job_id}.* 2>/dev/null"
        )
        if stdout:
            print(stdout)
    
    ssh_client.close()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_job_status.py <slurm_job_id>")
        print("Example: python check_job_status.py 2394707")
        sys.exit(1)
    
    job_id = sys.argv[1]
    check_job_status(job_id)
