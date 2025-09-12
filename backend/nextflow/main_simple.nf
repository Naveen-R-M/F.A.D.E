#!/usr/bin/env nextflow

/*
 * F.A.D.E (Fully Agentic Drug Engine) - Simple Test Workflow
 * Testing RCSB functionality without LLM dependencies
 */

nextflow.enable.dsl=2

// Import simple workflow
include { SIMPLE_TARGET_SELECTION } from './modules/simple_rcsb_target_selection'

// Workflow parameters
params.query = "Design a drug for EGFR targeting lung cancer"
params.output_dir = "${launchDir}/results_${new Date().format('yyyyMMdd_HHmmss')}"
params.fade_env_path = "${System.getenv('SCRATCH') ?: '/scratch/rajagopalmohanraj.n'}/conda-envs/fade"

// Validation
if (!params.query) {
    error "Error: --query parameter is required"
}

// Print workflow header
log.info """
╔══════════════════════════════════════════════════════════════════╗
║                F.A.D.E SIMPLE TEST PIPELINE v1.0                ║
║                   (No LLM - Direct RCSB Only)                   ║
╚══════════════════════════════════════════════════════════════════╝
Query            :  ${params.query}
Output directory :  ${params.output_dir}
Work Directory   :  ${projectDir}
======================================================================
"""

// Main workflow
workflow FADE_SIMPLE {
    
    // Create input channel from query
    query_ch = Channel.value(params.query)
    
    // Step 1: Simple Target Selection (No LLM)
    SIMPLE_TARGET_SELECTION(query_ch)
}

// Entry point
workflow {
    FADE_SIMPLE()
}

// Workflow completion handler
workflow.onComplete {
    def status = workflow.success ? 'SUCCESS' : 'FAILED'
    def duration = workflow.duration
    
    log.info """
    =====================================
    Simple RCSB Test Completed!
    Status: ${status}
    Duration: ${duration}
    Results: ${params.output_dir}
    =====================================
    """
}
