/*
 * Lead Optimization Module
 * Refines and optimizes top-scoring molecules
 */

process leadOptimization {
    tag "${top_hits.baseName}"
    label 'process_medium'
    publishDir "${params.output_dir}/06_lead_optimization", mode: 'copy'
    
    conda "${params.fade_env_path}"
    
    cpus 4
    memory '8.GB'
    time '2.h'
    
    input:
    path top_hits
    path structure
    
    output:
    path 'optimized_leads.json', emit: optimized_leads
    path 'optimization_stats.json', emit: stats
    
    script:
    def api_key = params.gemini_api_key ?: System.getenv('GEMINI_API_KEY') ?: ""
    def model = params.gemini_model ?: "models/gemini-2.5-flash"
    """
    # Set up environment
    export PYTHONPATH="${projectDir}:\$PYTHONPATH"
    export GEMINI_API_KEY="${api_key}"
    
    # Run lead optimization
    run_lead_optimization.py \\
        --top-hits ${top_hits} \\
        --structure ${structure} \\
        --output-dir . \\
        --api-key "${api_key}" \\
        --model "${model}"
    
    # Verify outputs exist
    if [ ! -f "optimized_leads.json" ]; then
        echo '{"error": "Lead optimization failed", "leads": []}' > optimized_leads.json
    fi
    
    if [ ! -f "optimization_stats.json" ]; then
        echo '{"input_leads": 0, "optimized_leads": 0}' > optimization_stats.json
    fi
    """
    
    stub:
    """
    echo '[{"molecule_id": "OPT_001", "optimized_score": -9.0, "smiles": "CCO"}]' > optimized_leads.json
    echo '{"input_leads": 1, "optimized_leads": 1, "average_improvement": 0.5}' > optimization_stats.json
    """
}

workflow LEAD_OPTIMIZATION {
    take:
    top_hits
    structure
    
    main:
    leadOptimization(top_hits, structure)
    
    emit:
    optimized_leads = leadOptimization.out.optimized_leads
    stats = leadOptimization.out.stats
}
