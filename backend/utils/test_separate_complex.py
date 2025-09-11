#!/usr/bin/env python

"""
Simple test script for separate_complex.py
"""

import os
import subprocess
import sys

def test_separate_complex():
    """Test the separate_complex.py script with example usage."""
    
    print("Testing separate_complex.py script...")
    
    # Check if the script exists
    script_path = os.path.join(os.path.dirname(__file__), 'separate_complex.py')
    
    if not os.path.exists(script_path):
        print(f"Error: Script not found at {script_path}")
        return False
    
    # Test help functionality
    print("\n1. Testing help option...")
    try:
        result = subprocess.run([sys.executable, script_path, '--help'], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            print("✓ Help option works correctly")
            print("Help output preview:")
            print(result.stdout[:200] + "..." if len(result.stdout) > 200 else result.stdout)
        else:
            print("✗ Help option failed")
            return False
    except Exception as e:
        print(f"✗ Error testing help: {e}")
        return False
    
    # Test with missing file (should fail gracefully)
    print("\n2. Testing with non-existent file...")
    try:
        result = subprocess.run([sys.executable, script_path, 'nonexistent.pdb'], 
                              capture_output=True, text=True)
        if result.returncode != 0 and "not found" in result.stdout:
            print("✓ Correctly handles missing input file")
        else:
            print("✗ Should fail with missing file")
            return False
    except Exception as e:
        print(f"✗ Error testing missing file: {e}")
        return False
    
    print("\n✓ Basic functionality tests passed!")
    print("\nTo test with real PDB files, run:")
    print(f"python {script_path} your_complex.pdb")
    print(f"python {script_path} your_complex.pdb --outprefix my_structure")
    print(f"python {script_path} your_complex.pdb --include-cofactors")
    
    return True

if __name__ == "__main__":
    test_separate_complex()
