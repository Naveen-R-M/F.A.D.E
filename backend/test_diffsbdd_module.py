"""
Test DiffSBDD module and command availability on HPC.
"""

from fade.tools.ssh_client import get_hpc_client
from fade.utils import get_logger

logger = get_logger("test.diffsbdd_module")

def test_diffsbdd_module():
    """Test if DiffSBDD module is available and working."""
    ssh_client = get_hpc_client()
    
    if not ssh_client.connect():
        print("❌ Failed to connect to HPC")
        return
    
    print("✓ Connected to HPC\n")
    
    # Test 1: Check if module command works
    print("1. Testing module command...")
    stdout, stderr, exit_code = ssh_client.execute_command("module avail 2>&1 | head -20")
    if exit_code == 0:
        print("✓ Module command works")
    else:
        print(f"❌ Module command failed: {stderr}")
        return
    
    # Test 2: Check if module path exists
    print("\n2. Checking module path...")
    module_path = "/projects/SimBioSys/share/software/modulefiles"
    stdout, stderr, exit_code = ssh_client.execute_command(f"ls {module_path}")
    if exit_code == 0:
        print(f"✓ Module path exists: {module_path}")
        print(f"Contents:\n{stdout[:500]}")
    else:
        print(f"❌ Module path not found: {module_path}")
    
    # Test 3: Check if diffsbdd module exists
    print("\n3. Checking for diffsbdd module...")
    stdout, stderr, exit_code = ssh_client.execute_command(
        f"module use {module_path} && module avail diffsbdd 2>&1"
    )
    print(f"Output:\n{stdout}")
    if "diffsbdd" in stdout.lower():
        print("✓ Found diffsbdd module")
    else:
        print("❌ diffsbdd module not found")
    
    # Test 4: Try loading module and checking commands
    print("\n4. Testing module load and commands...")
    test_cmd = f"""
module use {module_path}
module load diffsbdd
which sbatch_diffsbdd_generate
type sbatch_diffsbdd_generate
"""
    stdout, stderr, exit_code = ssh_client.execute_command(test_cmd)
    print(f"Exit code: {exit_code}")
    print(f"Output:\n{stdout}")
    if stderr:
        print(f"Stderr:\n{stderr}")
    
    # Test 5: Check checkpoint file
    print("\n5. Checking for checkpoint files...")
    stdout, stderr, exit_code = ssh_client.execute_command(
        "find /projects/SimBioSys -name '*crossdocked*.ckpt' 2>/dev/null | head -5"
    )
    if stdout.strip():
        print(f"✓ Found checkpoint files:\n{stdout}")
    else:
        print("❌ No checkpoint files found")
    
    # Test 6: Check recent DiffSBDD jobs
    print("\n6. Checking recent DiffSBDD jobs...")
    stdout, stderr, exit_code = ssh_client.execute_command(
        "sacct -u $(whoami) --format=JobID,JobName,State,Start -n | grep -i diff | tail -10"
    )
    if stdout.strip():
        print(f"Recent DiffSBDD jobs:\n{stdout}")
    else:
        print("No recent DiffSBDD jobs found")
    
    ssh_client.close()
    print("\n=== Test complete ===")


if __name__ == "__main__":
    test_diffsbdd_module()
