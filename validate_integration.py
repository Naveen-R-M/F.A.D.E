#!/usr/bin/env python3
"""
Validation script to test F.A.D.E Nextflow integration
"""

import os
import sys
import json
from pathlib import Path

def test_imports():
    """Test that all agent imports work correctly"""
    print("Testing Python agent imports...")
    
    try:
        from agents.target_selector.target_selector import TargetSelector
        print("✅ TargetSelector import successful")
    except Exception as e:
        print(f"❌ TargetSelector import failed: {e}")
        return False
    
    try:
        from agents.structure_predictor.structure_predictor import StructurePredictor
        print("✅ StructurePredictor import successful")
    except Exception as e:
        print(f"❌ StructurePredictor import failed: {e}")
        return False
    
    try:
        from agents.molecule_generator.molecule_generator import MoleculeGenerator
        print("✅ MoleculeGenerator import successful")
    except Exception as e:
        print(f"❌ MoleculeGenerator import failed: {e}")
        return False
    
    try:
        from utils.gemini_client import GeminiClient
        from utils.uniprot_client import UniProtClient
        print("✅ Utility imports successful")
    except Exception as e:
        print(f"❌ Utility imports failed: {e}")
        return False
    
    return True

def test_wrapper_scripts():
    """Test that wrapper scripts exist and are executable"""
    print("\nTesting wrapper scripts...")
    
    wrapper_scripts = [
        "nextflow/bin/run_target_selector.py",
        "nextflow/bin/run_structure_predictor.py", 
        "nextflow/bin/run_binding_site_analysis.py",
        "nextflow/bin/run_molecule_generator.py",
        "nextflow/bin/run_docking.py",
        "nextflow/bin/run_lead_optimization.py",
        "nextflow/bin/run_reporting.py"
    ]
    
    all_good = True
    for script in wrapper_scripts:
        if os.path.exists(script) and os.access(script, os.X_OK):
            print(f"✅ {script} exists and is executable")
        else:
            print(f"❌ {script} missing or not executable")
            all_good = False
    
    return all_good

def test_nextflow_files():
    """Test that Nextflow files are properly structured"""
    print("\nTesting Nextflow files...")
    
    required_files = [
        "nextflow/main.nf",
        "nextflow/nextflow.config",
        "nextflow/modules/target_selection.nf",
        "nextflow/modules/structure_prediction.nf",
        "nextflow/modules/binding_site.nf",
        "nextflow/modules/molecule_generation.nf",
        "nextflow/modules/docking.nf",
        "nextflow/modules/lead_optimization.nf",
        "nextflow/modules/reporting.nf"
    ]
    
    all_good = True
    for file in required_files:
        if os.path.exists(file):
            print(f"✅ {file} exists")
        else:
            print(f"❌ {file} missing")
            all_good = False
    
    return all_good

def test_environment():
    """Test environment configuration"""
    print("\nTesting environment configuration...")
    
    # Check SCRATCH directory
    scratch_dir = os.environ.get('SCRATCH')
    if scratch_dir:
        print(f"✅ SCRATCH directory: {scratch_dir}")
    else:
        print("❌ SCRATCH environment variable not set")
        return False
    
    # Check conda environment
    conda_env = f"{scratch_dir}/conda-envs/fade"
    if os.path.exists(conda_env):
        print(f"✅ Conda environment exists: {conda_env}")
    else:
        print(f"⚠️ Conda environment not found: {conda_env}")
        print("   This is expected if you haven't created it yet")
    
    # Check .env file
    if os.path.exists(".env"):
        print("✅ .env file exists")
        # Check for API key (without revealing it)
        with open(".env", "r") as f:
            env_content = f.read()
            if "GEMINI_API_KEY" in env_content and "=" in env_content:
                print("✅ GEMINI_API_KEY configured in .env")
            else:
                print("⚠️ GEMINI_API_KEY not found in .env")
    else:
        print("⚠️ .env file not found")
    
    return True

def test_agent_functionality():
    """Test basic agent functionality (if environment allows)"""
    print("\nTesting basic agent functionality...")
    
    try:
        # Test TargetSelector with minimal input
        from agents.target_selector.target_selector import TargetSelector
        
        # This will fail without API key, but we can test object creation
        try:
            target_selector = TargetSelector()
            print("✅ TargetSelector object creation successful")
        except Exception as e:
            if "API key" in str(e) or "GEMINI_API_KEY" in str(e):
                print("⚠️ TargetSelector requires API key (expected)")
            else:
                print(f"❌ TargetSelector creation failed: {e}")
                return False
        
        return True
        
    except Exception as e:
        print(f"❌ Agent functionality test failed: {e}")
        return False

def main():
    """Run all validation tests"""
    print("=" * 60)
    print("F.A.D.E Nextflow Integration Validation")
    print("=" * 60)
    
    tests = [
        ("Python Imports", test_imports),
        ("Wrapper Scripts", test_wrapper_scripts),
        ("Nextflow Files", test_nextflow_files),
        ("Environment", test_environment),
        ("Agent Functionality", test_agent_functionality)
    ]
    
    results = []
    for test_name, test_func in tests:
        print(f"\n{'='*20} {test_name} {'='*20}")
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"❌ {test_name} test crashed: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)
    
    passed = 0
    total = len(results)
    
    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status:<8} {test_name}")
        if result:
            passed += 1
    
    print(f"\nOverall: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 All validations passed! The integration is ready to use.")
        print("\nNext steps:")
        print("1. Ensure conda environment is created and activated")
        print("2. Set GEMINI_API_KEY in environment")
        print("3. Run: ./run_nextflow_integrated.sh --dry-run 'Test query'")
        print("4. Run: ./run_nextflow_integrated.sh 'Your actual query'")
    else:
        print(f"\n⚠️ {total - passed} validation(s) failed. Please address the issues above.")
        return False
    
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
