/*
 * Binding Site Analysis Module
 * Identifies and analyzes potential binding sites in protein structures
 */

process bindingSiteAnalysis {
    tag "${structure.baseName}"
    label 'process_medium'
    publishDir "${params.output_dir}/03_binding_site_analysis", mode: 'copy'
    
    conda "${params.fade_env_path}"
    
    input:
    path structure
    path target_info
    
    output:
    path 'binding_sites.json', emit: binding_sites
    path 'site_analysis.json', emit: analysis
    
    script:
    def api_key = params.gemini_api_key ?: System.getenv('GEMINI_API_KEY') ?: ""
    def model = params.gemini_model ?: "models/gemini-2.5-flash"
    """
    # Set up environment
    export PYTHONPATH="${projectDir}:\$PYTHONPATH"
    export GEMINI_API_KEY="${api_key}"
    
    # Run binding site analysis
    run_binding_site_analysis.py \\
        --structure ${structure} \\
        --target-info ${target_info} \\
        --output-dir . \\
        --api-key "${api_key}" \\
        --model "${model}"
    
    # Create site analysis summary
    if [ -f "binding_sites.json" ]; then
        python3 -c "
import json
with open('binding_sites.json', 'r') as f:
    sites = json.load(f)
    
analysis = {
    'num_sites': len(sites.get('sites', [])),
    'method': sites.get('method', 'unknown'),
    'confidence': 'high' if len(sites.get('sites', [])) > 0 else 'low'
}

with open('site_analysis.json', 'w') as f:
    json.dump(analysis, f, indent=2)
"
    else
        echo '{"error": "Binding site analysis failed"}' > site_analysis.json
    fi
    """
    
    stub:
    """
    echo '{"sites": [{"site_id": "site_1", "center": [0,0,0], "radius": 10.0}]}' > binding_sites.json
    echo '{"num_sites": 1, "confidence": "medium"}' > site_analysis.json
    """
}

workflow BINDING_SITE_ANALYSIS {
    take:
    structure
    target_info
    
    main:
    bindingSiteAnalysis(structure, target_info)
    
    emit:
    binding_sites = bindingSiteAnalysis.out.binding_sites
    analysis = bindingSiteAnalysis.out.analysis
}
