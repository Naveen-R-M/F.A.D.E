/*
 * RCSB Target Selection Module
 * Processes natural language queries to identify protein targets using RCSB PDB
 */

process rcsbTargetSelection {
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
    def api_key = params.gemini_api_key ?: ""
    def model = params.gemini_model ?: "models/gemini-2.5-flash"
    """
    # Set up environment
    export PYTHONPATH="${projectDir}:\$PYTHONPATH"
    
    # Load .env file if it exists
    if [ -f "${projectDir}/.env" ]; then
        source "${projectDir}/.env"
    fi
    
    # Debug output
    echo "Running RCSB target selection with query: ${query}"
    echo "API key length: \${#GEMINI_API_KEY}"
    echo "Project directory: ${projectDir}"
    echo "Looking for script at: ${projectDir}/bin/run_rcsb_target_selector.py"
    ls -la "${projectDir}/bin/" || echo "bin directory not found"
    
    # Use the API key from parameters
    export GEMINI_API_KEY="${api_key}"
    
    # Run the RCSB target selector agent
    python "${projectDir}/bin/run_rcsb_target_selector.py" \\
        --query "${query}" \\
        --output-dir . \\
        --api-key "${api_key}" \\
        --model "${model}"
    
    # Verify outputs exist and log results
    if [ -f "target_info.json" ]; then
        echo "Target info created:"
        cat target_info.json
    else
        echo "ERROR: target_info.json not created"
        exit 1
    fi
    
    if [ -f "protein.fasta" ]; then
        echo "FASTA file created:"
        head -2 protein.fasta
    else
        echo "ERROR: protein.fasta not created"
        exit 1
    fi
    
    if [ -f "requirements.json" ]; then
        echo "Requirements created:"
        cat requirements.json
    else
        echo "ERROR: requirements.json not created"
        exit 1
    fi
    
    # Check if PDB structure was downloaded
    if [ -f "structure.pdb" ]; then
        echo "PDB structure downloaded:"
        head -5 structure.pdb | grep "^HEADER\\|^TITLE\\|^COMPND"
    else
        echo "ERROR: No PDB structure file found - RCSB target selection failed"
        exit 1
    fi
    """
    
    stub:
    """
    # Create stub outputs for testing
    echo '{"target": "KRAS", "pdb_id": "7OM5", "source": "rcsb_pdb", "confidence": 1.0}' > target_info.json
    echo '>KRAS|7OM5|GTPase KRas from RCSB PDB' > protein.fasta
    echo 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMSCKCVLS' >> protein.fasta
    echo '{"binding_affinity": "< -8 kcal/mol", "selectivity": "high", "blood_brain_barrier": "permeable"}' > requirements.json
    
    # Create a minimal PDB stub
    echo "HEADER    TRANSFERASE/DNA                         20-JUN-21   7OM5" > structure.pdb
    echo "TITLE     KRAS G12D MUTANT BOUND TO INHIBITOR" >> structure.pdb
    echo "COMPND    MOL_ID: 1;" >> structure.pdb
    echo "ATOM      1  CA  GLY A  10      24.384  14.020  18.292  1.00 41.24           C" >> structure.pdb
    echo "END" >> structure.pdb
    """
}

workflow TARGET_SELECTION {
    take:
    query
    
    main:
    rcsbTargetSelection(query)
    
    emit:
    target_info = rcsbTargetSelection.out.target_info
    fasta = rcsbTargetSelection.out.fasta
    structure_pdb = rcsbTargetSelection.out.structure_pdb
    pdb_id = rcsbTargetSelection.out.pdb_id
    requirements = rcsbTargetSelection.out.requirements
}
