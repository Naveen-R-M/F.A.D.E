/*
 * Target Selection Module
 * Processes natural language queries to identify protein targets using the actual TargetSelector agent
 */

process targetSelection {
    tag "$query"
    label 'process_medium'
    publishDir "${params.output_dir}/01_target_selection", mode: 'copy'
    
    conda "${params.fade_env_path}"
    
    input:
    val query
    
    output:
    path 'target_info.json', emit: target_info
    path 'protein.fasta', emit: fasta
    val 'none', emit: pdb_id
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
    
    # Debug output (simplified)
    echo "Running target selection with API key length: \${#GEMINI_API_KEY}"
    echo "Params API key length: ${api_key.length()}"
    
    # Use the API key from parameters (passed from Nextflow config)
    export GEMINI_API_KEY="${api_key}"
    
    # Run the target selector agent
    run_target_selector.py \\
        --query "${query}" \\
        --output-dir . \\
        --api-key "${api_key}" \\
        --model "${model}"
    
    # Verify outputs exist
    if [ ! -f "target_info.json" ]; then
        echo '{"error": "Target selection failed", "target": "unknown"}' > target_info.json
    fi
    
    if [ ! -f "protein.fasta" ]; then
        echo '>unknown|unknown|Unknown protein' > protein.fasta
        echo 'UNKNOWN' >> protein.fasta
    fi
    
    if [ ! -f "requirements.json" ]; then
        echo '{"binding_affinity": "< -8 kcal/mol"}' > requirements.json
    fi
    """
    
    stub:
    """
    echo '{"target": "KRAS", "uniprot_id": "P01116", "pdb_id": "none"}' > target_info.json
    echo '>KRAS|P01116|GTPase KRas' > protein.fasta
    echo 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMSCKCVLS' >> protein.fasta
    echo '{"binding_affinity": "< -8 kcal/mol"}' > requirements.json
    """
}

workflow TARGET_SELECTION {
    take:
    query
    
    main:
    targetSelection(query)
    
    emit:
    target_info = targetSelection.out.target_info
    fasta = targetSelection.out.fasta
    pdb_id = targetSelection.out.pdb_id
    requirements = targetSelection.out.requirements
}
