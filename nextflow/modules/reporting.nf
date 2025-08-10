/*
 * Reporting Module
 * Generates comprehensive final reports from all pipeline results
 */

process finalReporting {
    tag "final_report"
    label 'process_low'
    publishDir "${params.output_dir}/07_final_report", mode: 'copy'
    
    conda "${params.fade_env_path}"
    
    input:
    path target_info
    path structure
    path docking_results
    path optimized_leads
    
    output:
    path 'final_report.json', emit: report
    path 'analysis_summary.md', emit: summary
    path 'top_candidates.sdf', emit: top_candidates
    
    script:
    def api_key = params.gemini_api_key ?: System.getenv('GEMINI_API_KEY') ?: ""
    def model = params.gemini_model ?: "models/gemini-2.5-flash"
    """
    # Set up environment
    export PYTHONPATH="${projectDir}:\$PYTHONPATH"
    export GEMINI_API_KEY="${api_key}"
    
    # Run final reporting
    run_reporting.py \\
        --target-info ${target_info} \\
        --structure ${structure} \\
        --docking-results ${docking_results} \\
        --optimized-leads ${optimized_leads} \\
        --output-dir . \\
        --api-key "${api_key}" \\
        --model "${model}"
    
    # Verify outputs exist
    if [ ! -f "final_report.json" ]; then
        echo '{"error": "Report generation failed", "status": "failed"}' > final_report.json
    fi
    
    if [ ! -f "analysis_summary.md" ]; then
        echo "# F.A.D.E Report Generation Failed" > analysis_summary.md
        echo "An error occurred during report generation." >> analysis_summary.md
    fi
    
    if [ ! -f "top_candidates.sdf" ]; then
        echo "# No candidates available" > top_candidates.sdf
    fi
    """
    
    stub:
    """
    echo '{"status": "success", "top_candidates": 1}' > final_report.json
    echo "# F.A.D.E Drug Discovery Report" > analysis_summary.md
    echo "## Summary: 1 candidate identified" >> analysis_summary.md
    echo "MOL_001" > top_candidates.sdf
    echo "End of SDF" >> top_candidates.sdf
    """
}

workflow REPORTING {
    take:
    target_info
    structure
    docking_results
    optimized_leads
    
    main:
    finalReporting(target_info, structure, docking_results, optimized_leads)
    
    emit:
    report = finalReporting.out.report
    summary = finalReporting.out.summary
    top_candidates = finalReporting.out.top_candidates
}
