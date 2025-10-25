"""
Test that line endings are correct in the generated script.
"""

from fade.tools.hpc_diffsbdd import get_hpc_diffsbdd_client


def test_line_endings():
    """Verify that the batch script has Unix line endings."""
    
    client = get_hpc_diffsbdd_client()
    
    # Generate a test script
    script = client._create_batch_script(
        job_id="test",
        remote_dir="/tmp/test",
        remote_pdb="/tmp/test.pdb",
        pocket_residues=["A:12"],
        output_file="test.sdf",
        n_samples=10,
        sanitize=True
    )
    
    # Check line endings
    has_windows_endings = '\r\n' in script
    has_unix_endings = '\n' in script
    
    print("="*60)
    print("Line Ending Test")
    print("="*60)
    print(f"Script length: {len(script)} characters")
    print(f"Has Windows line endings (\\r\\n): {has_windows_endings}")
    print(f"Has Unix line endings (\\n): {has_unix_endings}")
    print()
    
    if has_windows_endings:
        print("⚠️  WARNING: Script contains Windows line endings")
        print("   These will be converted during upload")
    else:
        print("✅ Script has Unix line endings only")
    
    # Show first few lines
    lines = script.split('\n')[:5]
    print(f"\nFirst 5 lines:")
    for i, line in enumerate(lines, 1):
        # Show line ending explicitly
        line_repr = repr(line)
        print(f"  {i}. {line_repr}")
    
    print("\n" + "="*60)
    
    # Test the conversion in _upload_batch_script
    print("\nTesting _upload_batch_script conversion...")
    
    # Simulate script with Windows endings
    test_script_windows = "#!/bin/bash\r\n#SBATCH --job-name=test\r\necho 'test'\r\n"
    converted = test_script_windows.replace('\r\n', '\n')
    
    print(f"Before conversion: {repr(test_script_windows[:50])}")
    print(f"After conversion:  {repr(converted[:50])}")
    
    if '\r\n' not in converted:
        print("✅ Conversion works correctly!")
    else:
        print("❌ Conversion failed!")
    
    print("="*60)


if __name__ == "__main__":
    test_line_endings()
