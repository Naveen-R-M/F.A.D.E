#!/usr/bin/env nextflow

/*
 * F.A.D.E (Fully Agentic Drug Engine) - Nextflow Workflow
 * Main workflow orchestrating all drug discovery pipeline stages
 */

nextflow.enable.dsl=2

// Import sub-workflows
include { TARGET_SELECTION } from './modules/target_selection'
include { STRUCTURE_PREDICTION } from './modules/structure_prediction'
include { BINDING_SITE_ANALYSIS } from './modules/binding_site'
include { MOLECULE_GENERATION } from './modules/molecule_generation'
include { DOCKING } from './modules/docking'
include { LEAD_OPTIMIZATION } from './modules/lead_optimization'
include { REPORTING } from './modules/reporting'

// Workflow parameters with defaults
params.query = "Design a drug for EGFR targeting lung cancer"
params.output_dir = "${launchDir}/results"
params.max_molecules = 1000
params.max_iterations = 3
params.use_alphafold = true
params.docking_method = "vina"
params.parallel_jobs = 10
params.resume_from = null

// Resource parameters
params.max_cpus = 16
params.max_memory = '128.GB'
params.max_time = '240.h'

// Container/Module paths
params.alphafold_container = "/shared/containers/alphafold3.sif"
params.schrodinger_module = "schrodinger/2024-4"
params.rdkit_env = "rdkit-env"

// Print workflow header
log.info """
╔══════════════════════════════════════════════════════════════════╗
║                  F.A.D.E NEXTFLOW PIPELINE v1.0                  ║
╚══════════════════════════════════════════════════════════════════╝
Query           : ${params.query}
Output directory: ${params.output_dir}
Max molecules   : ${params.max_molecules}
Max iterations  : ${params.max_iterations}
Parallel jobs   : ${params.parallel_jobs}
=====================================
"""

// Main workflow
workflow FADE {
    
    // Create input channel from query
    query_ch = Channel.value(params.query)
    
    // Step 1: Target Selection
    // Identifies protein target from natural language query
    TARGET_SELECTION(query_ch)
    
    // Step 2: Structure Prediction/Retrieval
    // Either fetches from PDB or predicts using AlphaFold3
    STRUCTURE_PREDICTION(
        TARGET_SELECTION.out.target_info,
        TARGET_SELECTION.out.fasta,
        TARGET_SELECTION.out.pdb_id
    )
    
    // Step 3: Binding Site Analysis
    // Analyzes pockets and binding sites
    BINDING_SITE_ANALYSIS(
        STRUCTURE_PREDICTION.out.structure,
        TARGET_SELECTION.out.target_info
    )
    
    // Step 4: Molecule Generation
    // Generate molecules (simplified - single iteration for now)
    iteration_ch = Channel.value(0)
    MOLECULE_GENERATION(
        TARGET_SELECTION.out.requirements,
        BINDING_SITE_ANALYSIS.out.binding_sites,
        iteration_ch
    )
    
    // Step 5: Molecular Docking
    // Performs high-throughput virtual screening
    DOCKING(
        MOLECULE_GENERATION.out.molecules,
        STRUCTURE_PREDICTION.out.prepared_receptor,
        BINDING_SITE_ANALYSIS.out.binding_sites
    )
    
    // Step 6: Lead Optimization
    // Optimizes top hits
    LEAD_OPTIMIZATION(
        DOCKING.out.top_hits,
        STRUCTURE_PREDICTION.out.structure
    )
    
    // Step 7: Generate Report
    // Creates comprehensive report
    REPORTING(
        TARGET_SELECTION.out.target_info,
        STRUCTURE_PREDICTION.out.structure,
        DOCKING.out.all_results,
        LEAD_OPTIMIZATION.out.optimized_leads
    )
}

// Entry point
workflow {
    FADE()
}

// Workflow completion handler
workflow.onComplete {
    log.info """
    =====================================
    Pipeline completed!
    Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration: ${workflow.duration}
    Results: ${params.output_dir}
    =====================================
    """
    
    if (!workflow.success) {
        log.error "Pipeline execution failed. Check the logs for details."
    }
}

// Error handler
workflow.onError {
    log.error "Pipeline error: ${workflow.errorMessage}"
    log.error "Error report saved to: ${params.output_dir}/error_report.txt"
}
