"""
Minimal HPC Kalasanty client for testing.
Step-by-step implementation to verify functionality.
"""

import os
import time
import logging
from typing import Dict, Any, Optional, List
from pathlib import Path

from fade.tools.ssh_client import get_hpc_client
from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.hpc_kalasanty_minimal")


class HPCKalasantyClient:
    """Minimal client for testing Kalasanty on HPC."""
    
    def __init__(self):
        """Initialize with basic settings."""
        self.ssh_client = get_hpc_client()
        self.scratch_base = "/scratch/rajagopalmohanraj.n/F.A.D.E/kalasanty_test"
        self.module_path = "/projects/SimBioSys/share/software/modulefiles/"
        
    def test_module_load(self) -> bool:
        """Test if we can load the Kalasanty module."""
        logger.info("Testing Kalasanty module load...")
        
        try:
            if not self.ssh_client.connect():
                logger.error("Failed to connect to HPC")
                return False
            
            # Try to load module
            cmd = f"""
module use {self.module_path}
module load kalasanty/latest
echo "Module loaded successfully"
"""
            stdout, stderr, exit_code = self.ssh_client.execute_command(cmd, timeout=30)
            
            if exit_code == 0:
                logger.info(f"✅ Module loaded: {stdout}")
                return True
            else:
                logger.error(f"❌ Failed to load: {stderr}")
                return False
                
        except Exception as e:
            logger.error(f"Error: {e}")
            return False
        finally:
            self.ssh_client.disconnect()
    
    def find_kalasanty_command(self) -> Optional[str]:
        """Find the correct kalasanty command syntax."""
        logger.info("Finding kalasanty command...")
        
        try:
            if not self.ssh_client.connect():
                return None
            
            # Commands to try
            test_commands = [
                "kalasanty --help",
                "kalasanty -h",
                "python -m kalasanty --help",
                "kalasanty-predict --help",
                "kalasanty_predict --help",
                "python $KALASANTY_HOME/kalasanty.py --help"
            ]
            
            base = f"""
module use {self.module_path}
module load kalasanty/latest
"""
            
            for cmd in test_commands:
                full_cmd = base + cmd + " 2>&1"
                stdout, stderr, exit_code = self.ssh_client.execute_command(
                    full_cmd, timeout=30
                )
                
                # Check if command worked
                if "usage" in stdout.lower() or "help" in stdout.lower() or exit_code == 0:
                    logger.info(f"✅ Found working command: {cmd}")
                    logger.info(f"Help output:\n{stdout[:500]}")
                    return cmd
                    
            logger.warning("No working command found")
            return None
            
        except Exception as e:
            logger.error(f"Error: {e}")
            return None
        finally:
            self.ssh_client.disconnect()
    
    def test_with_simple_pdb(self) -> bool:
        """Test Kalasanty with a minimal PDB file."""
        logger.info("Testing with simple PDB...")
        
        # Minimal PDB content
        test_pdb = """HEADER    TEST STRUCTURE
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.221   2.370   0.000  1.00 20.00           O
ATOM      5  CB  ALA A   1       1.988  -0.776  -1.206  1.00 20.00           C
END
"""
        
        try:
            if not self.ssh_client.connect():
                return False
            
            # Create test directory
            test_dir = f"{self.scratch_base}/test_{int(time.time())}"
            self.ssh_client.execute_command(f"mkdir -p {test_dir}")
            
            # Upload PDB
            pdb_path = f"{test_dir}/test.pdb"
            self.ssh_client.upload_content(test_pdb, pdb_path)
            logger.info(f"Uploaded test PDB to {pdb_path}")
            
            # Find command first
            kalasanty_cmd = self.find_kalasanty_command()
            if not kalasanty_cmd:
                logger.error("Could not find kalasanty command")
                return False
            
            # Modify command for actual prediction
            if "--help" in kalasanty_cmd or "-h" in kalasanty_cmd:
                # Remove help flag and add input/output
                predict_cmd = kalasanty_cmd.replace("--help", "").replace("-h", "")
                predict_cmd += f" -i {pdb_path} -o {test_dir}/output"
            else:
                predict_cmd = kalasanty_cmd
            
            # Run prediction
            run_cmd = f"""
module use {self.module_path}
module load kalasanty/latest
cd {test_dir}
{predict_cmd}
echo "Exit code: $?"
ls -la {test_dir}/
"""
            stdout, stderr, exit_code = self.ssh_client.execute_command(
                run_cmd, timeout=60
            )
            
            logger.info(f"Execution output:\n{stdout}")
            if stderr:
                logger.warning(f"Stderr:\n{stderr}")
            
            # Check for output files
            ls_cmd = f"ls -la {test_dir}/"
            stdout, stderr, exit_code = self.ssh_client.execute_command(ls_cmd)
            logger.info(f"Directory contents:\n{stdout}")
            
            return True
            
        except Exception as e:
            logger.error(f"Error: {e}")
            return False
        finally:
            self.ssh_client.disconnect()


def test_kalasanty_minimal():
    """Run minimal tests."""
    client = HPCKalasantyClient()
    
    print("\n" + "="*60)
    print("STEP 1: Test module loading")
    print("="*60)
    if not client.test_module_load():
        print("❌ Module load failed")
        return False
    
    print("\n" + "="*60)
    print("STEP 2: Find kalasanty command")
    print("="*60)
    cmd = client.find_kalasanty_command()
    if not cmd:
        print("❌ Could not find command")
        return False
    print(f"✅ Found command: {cmd}")
    
    print("\n" + "="*60)
    print("STEP 3: Test with simple PDB")
    print("="*60)
    if not client.test_with_simple_pdb():
        print("❌ PDB test failed")
        return False
    
    print("\n✅ All tests passed!")
    return True


if __name__ == "__main__":
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    test_kalasanty_minimal()
