"""
Check the generate_ligands.py script to understand its expected arguments.
"""

from fade.tools.ssh_client import get_hpc_client


def check_generate_ligands_script():
    """Examine the generate_ligands.py script to verify argument format."""
    
    ssh = get_hpc_client()
    if not ssh.connect():
        print("❌ Failed to connect")
        return
    
    print("="*60)
    print("Checking generate_ligands.py")
    print("="*60 + "\n")
    
    diffsbdd_path = "/projects/SimBioSys/share/software/DiffSBDD"
    script_path = f"{diffsbdd_path}/generate_ligands.py"
    
    # Check if script exists
    print("1. Checking if script exists...")
    stdout, _, code = ssh.execute_command(f"test -f {script_path} && echo EXISTS")
    if "EXISTS" in stdout:
        print(f"✓ Script found: {script_path}")
    else:
        print(f"❌ Script not found: {script_path}")
        return
    
    # Get help/usage information
    print("\n2. Getting script usage information...")
    cmd = f'''
module use /projects/SimBioSys/share/software/modulefiles
module load diffsbdd
eval "$(micromamba shell hook --shell bash)"
micromamba activate {diffsbdd_path}/env/
python {script_path} --help
'''
    
    stdout, stderr, code = ssh.execute_command(cmd, timeout=30)
    
    if code == 0 or stdout:
        print("Script help output:")
        print("-"*60)
        print(stdout if stdout else stderr)
        print("-"*60)
    else:
        print(f"❌ Failed to get help: {stderr}")
        
        # Try to read first 50 lines of the script
        print("\n3. Reading script header to understand arguments...")
        stdout, _, _ = ssh.execute_command(f"head -100 {script_path}")
        if stdout:
            print("Script header:")
            print("-"*60)
            print(stdout)
            print("-"*60)
    
    # Check for checkpoint file
    print("\n4. Checking for checkpoint file...")
    checkpoint_path = f"{diffsbdd_path}/checkpoints/crossdocked_fullatom_cond.ckpt"
    stdout, _, code = ssh.execute_command(f"test -f {checkpoint_path} && echo EXISTS")
    if "EXISTS" in stdout:
        print(f"✓ Checkpoint found: {checkpoint_path}")
    else:
        print(f"❌ Checkpoint not found: {checkpoint_path}")
        
        # Search for checkpoints
        print("   Searching for checkpoints...")
        stdout, _, _ = ssh.execute_command(
            f"find {diffsbdd_path} -name '*.ckpt' -type f 2>/dev/null"
        )
        if stdout:
            print(f"   Found checkpoints:\n{stdout}")
    
    ssh.close()
    
    print("\n" + "="*60)
    print("Check complete")
    print("="*60)


if __name__ == "__main__":
    check_generate_ligands_script()
