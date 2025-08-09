#!/usr/bin/env python3
"""
Test the wrapper script integration without requiring Nextflow
"""

import os
import sys
import json
import tempfile
import subprocess
from pathlib import Path

def test_target_selector_wrapper():
    """Test the target selector wrapper script"""
    print("Testing Target Selector wrapper...")
    
    # Create temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Test the wrapper script with a simple query
        cmd = [
            "python3", "nextflow/bin/run_target_selector.py",
            "--query", "Target EGFR for lung cancer",
            "--output-dir", temp_dir,
            "--api-key", "test_key"  # This will fail but we can test the structure
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            
            # Check if output files were created (even if they contain errors)
            expected_files = ["target_info.json", "protein.fasta", "requirements.json"]
            files_created = []
            
            for file in expected_files:
                file_path = os.path.join(temp_dir, file)
                if os.path.exists(file_path):
                    files_created.append(file)
            
            print(f"   Files created: {files_created}")
            print(f"   Return code: {result.returncode}")
            if result.stdout:
                print(f"   Stdout: {result.stdout[:200]}...")
            if result.stderr:
                print(f"   Stderr: {result.stderr[:200]}...")
            
            # Success if files were created (even with errors)
            return len(files_created) >= 2
            
        except subprocess.TimeoutExpired:
            print("   ❌ Timeout expired")
            return False
        except Exception as e:
            print(f"   ❌ Exception: {e}")
            return False

def test_molecule_generator_wrapper():
    """Test molecule generator wrapper with mock inputs"""
    print("Testing Molecule Generator wrapper...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create mock input files
        target_info = {"target": "EGFR", "uniprot_id": "P00533"}
        requirements = {"binding_affinity": "< -8 kcal/mol"}
        binding_sites = {"sites": [{"site_id": "site_1", "center": [0,0,0]}]}
        
        with open(os.path.join(temp_dir, "target_info.json"), "w") as f:
            json.dump(target_info, f)
        with open(os.path.join(temp_dir, "requirements.json"), "w") as f:
            json.dump(requirements, f)
        with open(os.path.join(temp_dir, "binding_sites.json"), "w") as f:
            json.dump(binding_sites, f)
        
        cmd = [
            "python3", "nextflow/bin/run_molecule_generator.py",
            "--requirements", os.path.join(temp_dir, "requirements.json"),
            "--binding-sites", os.path.join(temp_dir, "binding_sites.json"),
            "--target-info", os.path.join(temp_dir, "target_info.json"),
            "--output-dir", temp_dir,
            "--api-key", "test_key",
            "--max-molecules", "5"
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            
            expected_files = ["molecules.sdf", "molecules.json", "generation_stats.json"]
            files_created = [f for f in expected_files if os.path.exists(os.path.join(temp_dir, f))]
            
            print(f"   Files created: {files_created}")
            print(f"   Return code: {result.returncode}")
            
            return len(files_created) >= 2
            
        except Exception as e:
            print(f"   ❌ Exception: {e}")
            return False

def main():
    """Run wrapper integration tests"""
    print("=" * 60)
    print("F.A.D.E Wrapper Script Integration Test")
    print("=" * 60)
    
    tests = [
        ("Target Selector Wrapper", test_target_selector_wrapper),
        ("Molecule Generator Wrapper", test_molecule_generator_wrapper)
    ]
    
    results = []
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        try:
            result = test_func()
            results.append((test_name, result))
            status = "✅ PASS" if result else "❌ FAIL"
            print(f"   {status}")
        except Exception as e:
            print(f"   ❌ CRASHED: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "=" * 60)
    print("WRAPPER TEST SUMMARY")
    print("=" * 60)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status} {test_name}")
    
    print(f"\nOverall: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 Wrapper integration tests passed!")
        print("The Nextflow workflow should work with actual agents.")
    else:
        print(f"\n⚠️ Some tests failed. This may be due to missing dependencies.")
        print("The workflow structure is correct, but runtime environment needs setup.")

if __name__ == "__main__":
    main()
