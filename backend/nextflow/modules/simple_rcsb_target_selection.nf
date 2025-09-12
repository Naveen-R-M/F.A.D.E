/*
 * Simple RCSB Target Selection Module - No LLM Dependencies
 * Quick test version without Gemini API calls
 */

process simpleRcsbTargetSelection {
    tag "$query"
    label 'process_medium'
    publishDir "${params.output_dir}/01_target_selection", mode: 'copy'
    
    conda "${params.fade_env_path}"
    
    input:
    val query
    
    output:
    path 'target_info.json', emit: target_info
    path 'protein.fasta', emit: fasta
    path 'structure.pdb', optional: true, emit: structure_pdb
    val query, emit: pdb_id
    path 'requirements.json', emit: requirements
    
    script:
    """
    # Debug info
    echo "=== SIMPLE RCSB TARGET SELECTION ==="
    echo "Query: ${query}"
    echo "Working directory: \\$PWD"
    echo "Project directory: ${projectDir}"
    echo "====================================="
    
    # Set up environment
    export PYTHONPATH="${projectDir}:\\$PYTHONPATH"
    
    # Run the simple RCSB target selector (no API needed)
    python "${projectDir}/bin/run_simple_rcsb_selector.py" \\
        --query "${query}" \\
        --output-dir . \\
        --verbose
    
    # Verify outputs
    echo "=== CHECKING OUTPUTS ==="
    for file in target_info.json protein.fasta requirements.json; do
        if [ -f "\\$file" ]; then
            echo "✓ \\$file created"
            ls -la "\\$file"
        else
            echo "✗ \\$file missing"
            exit 1
        fi
    done
    
    if [ -f "structure.pdb" ]; then
        echo "✓ structure.pdb downloaded"
        head -3 structure.pdb
    else
        echo "✗ structure.pdb missing"
        exit 1
    fi
    
    echo "=== SIMPLE RCSB SELECTION COMPLETE ==="
    """
    
    stub:
    """
    echo '{"target": "KRAS", "pdb_id": "7OM5", "source": "simple_rcsb"}' > target_info.json
    echo '>KRAS|7OM5|KRAS from simple RCSB' > protein.fasta
    echo 'PLACEHOLDER_SEQUENCE' >> protein.fasta
    echo '{"binding_affinity": "< -8 kcal/mol"}' > requirements.json
    echo "HEADER    SIMPLE TEST STRUCTURE" > structure.pdb
    echo "END" >> structure.pdb
    """
}

workflow SIMPLE_TARGET_SELECTION {
    take:
    query
    
    main:
    simpleRcsbTargetSelection(query)
    
    emit:
    target_info = simpleRcsbTargetSelection.out.target_info
    fasta = simpleRcsbTargetSelection.out.fasta
    structure_pdb = simpleRcsbTargetSelection.out.structure_pdb
    pdb_id = simpleRcsbTargetSelection.out.pdb_id
    requirements = simpleRcsbTargetSelection.out.requirements
}
