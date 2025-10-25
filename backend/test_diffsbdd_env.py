"""
Test if the DiffSBDD environment activation works.
"""

from fade.tools.ssh_client import get_hpc_client
from fade.utils import get_logger

logger = get_logger("test.diffsbdd_env")


def test_diffsbdd_environment():
    """Test if DiffSBDD environment can be activated and used."""
    
    ssh = get_hpc_client()
    if not ssh.connect():
        print("❌ Failed to connect to HPC")
        return
    
    print("✓ Connected to HPC\n")
    
    module_path = "/projects/SimBioSys/share/software/modulefiles"
    env_path = "/projects/SimBioSys/share/software/DiffSBDD/env/"
    
    # Test 1: Check if micromamba is available
    print("1. Checking micromamba...")
    stdout, stderr, code = ssh.execute_command("which micromamba")
    if code == 0 and stdout.strip():
        print(f"✓ micromamba found: {stdout.strip()}")
    else:
        print("❌ micromamba not found")
        print(f"   You may need to: module load micromamba")
        return
    
    # Test 2: Check if environment exists
    print("\n2. Checking DiffSBDD environment...")
    stdout, stderr, code = ssh.execute_command(f"test -d {env_path} && echo EXISTS")
    if "EXISTS" in stdout:
        print(f"✓ Environment directory exists: {env_path}")
    else:
        print(f"❌ Environment directory not found: {env_path}")
        return
    
    # Test 3: Try to activate environment and check torch
    print("\n3. Testing environment activation and torch import...")
    test_cmd = f'''
module use {module_path}
module load diffsbdd
eval "$(micromamba shell hook --shell bash)"
micromamba activate {env_path}
python -c "import torch; print('PyTorch version:', torch.__version__)"
'''
    
    stdout, stderr, code = ssh.execute_command(test_cmd)
    print(f"Exit code: {code}")
    if code == 0 and "PyTorch version" in stdout:
        print(f"✓ Environment activation works!")
        print(f"   {stdout.strip()}")
    else:
        print("❌ Environment activation failed")
        print(f"   stdout: {stdout}")
        print(f"   stderr: {stderr}")
    
    # Test 4: Check what's in the environment
    print("\n4. Checking environment contents...")
    list_cmd = f'''
eval "$(micromamba shell hook --shell bash)"
micromamba activate {env_path}
micromamba list | grep -E "(torch|rdkit|numpy)"
'''
    
    stdout, stderr, code = ssh.execute_command(list_cmd)
    if stdout.strip():
        print(f"✓ Key packages:\n{stdout}")
    
    # Test 5: Check if sbatch_diffsbdd_generate exists and what it does
    print("\n5. Examining sbatch_diffsbdd_generate script...")
    check_cmd = f'''
module use {module_path}
module load diffsbdd
which sbatch_diffsbdd_generate
'''
    
    stdout, stderr, code = ssh.execute_command(check_cmd)
    if code == 0 and stdout.strip():
        script_path = stdout.strip()
        print(f"✓ Script found: {script_path}")
        
        # Read the script to see what it does
        stdout, stderr, code = ssh.execute_command(f"cat {script_path}")
        if code == 0:
            print(f"\n   Script contents:\n{stdout}")
    else:
        print("❌ sbatch_diffsbdd_generate not found")
    
    ssh.close()
    print("\n" + "="*60)
    print("Environment test complete")
    print("="*60)


if __name__ == "__main__":
    test_diffsbdd_environment()
