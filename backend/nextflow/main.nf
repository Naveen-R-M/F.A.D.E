#!/usr/bin/env nextflow

/*
 * F.A.D.E (Fully Agentic Drug Engine) - Minimal Nextflow Workflow
 * Starting with RCSB Target Selection only
 */

nextflow.enable.dsl=2

// Import sub-workflows
include { TARGET_SELECTION } from './modules/rcsb_target_selection'

// Workflow parameters with defaults
params.query = "Design a drug for EGFR targeting lung cancer"
params.output_dir = "${launchDir}/results_${new Date().format('yyyyMMdd_HHmmss')}"

// Environment and API configuration
params.fade_env_path = "${System.getenv('SCRATCH')}/conda-envs/fade"
params.gemini_api_key = System.getenv('GEMINI_API_KEY') ?: ""
params.gemini_model = System.getenv('GEMINI_MODEL') ?: "models/gemini-2.5-flash"

// Validation
if (!params.query) {
    error "Error: --query parameter is required"
}

if (!params.gemini_api_key) {
    log.warn "Warning: GEMINI_API_KEY not set. Target selection may fail."
}

// Print workflow header
log.info """
╔══════════════════════════════════════════════════════════════════╗
║                  F.A.D.E MINIMAL PIPELINE v1.0                  ║
║                    (RCSB Target Selection Only)                 ║
╚══════════════════════════════════════════════════════════════════╝
Query            :  ${params.query}
Output directory :  ${params.output_dir}
Conda env        :  ${params.fade_env_path}
Work Directory   :  ${projectDir}
======================================================================
"""

// Main workflow
workflow FADE {
    
    // Create input channel from query
    query_ch = Channel.value(params.query)
    
    // Step 1: Target Selection (RCSB-based)
    TARGET_SELECTION(query_ch)
}

// Entry point
workflow {
    FADE()
}

// Workflow completion handler
workflow.onComplete {
    def status = workflow.success ? 'SUCCESS' : 'FAILED'
    def duration = workflow.duration
    
    log.info """
    =====================================
    F.A.D.E Pipeline Completed!
    Status: ${status}
    Duration: ${duration}
    Results: ${params.output_dir}
    Query: "${params.query}"
    =====================================
    """
    
    if (workflow.success) {
        // Create completion summary
        def summaryFile = new File("${params.output_dir}/pipeline_summary.txt")
        summaryFile.text = """
F.A.D.E Pipeline Execution Summary
=================================
Query: ${params.query}
Status: ${status}
Duration: ${duration}
Start time: ${workflow.start}
Completion time: ${workflow.complete}
Results directory: ${params.output_dir}

Key outputs:
- 01_target_selection/target_info.json
- 01_target_selection/structure.pdb (from RCSB)
- 01_target_selection/protein.fasta
- 01_target_selection/requirements.json

To view target info: cat ${params.output_dir}/01_target_selection/target_info.json
"""
        log.info "Pipeline summary written to: ${params.output_dir}/pipeline_summary.txt"
    } else {
        log.error "Pipeline execution failed. Check the logs for details."
        log.error "Error message: ${workflow.errorMessage}"
    }
}

// Error handler
workflow.onError {
    log.error "Pipeline error: ${workflow.errorMessage}"
    log.error "Error report will be saved to: ${params.output_dir}/error_report.txt"
    
    // Create error report file
    def errorFile = new File("${params.output_dir}/error_report.txt")
    errorFile.text = """
F.A.D.E Pipeline Error Report
============================
Timestamp: ${new Date()}
Query: ${params.query}
Error: ${workflow.errorMessage}
Duration before failure: ${workflow.duration}

Troubleshooting steps:
1. Check that the conda environment is properly activated: ${params.fade_env_path}
2. Verify GEMINI_API_KEY is set in environment
3. Check individual process logs in work/ directory

Command to resume from failure:
nextflow run main.nf --query "${params.query}" --resume
"""
}
