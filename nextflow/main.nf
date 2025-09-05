#!/usr/bin/env nextflow

/*
 * F.A.D.E (Fully Agentic Drug Engine) - Nextflow Workflow
 * Main workflow orchestrating all drug discovery pipeline stages using actual Python agents
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
params.output_dir = "${launchDir}/results_${new Date().format('yyyyMMdd_HHmmss')}"
params.max_molecules = 100
params.max_iterations = 1
params.use_alphafold = true
params.docking_method = "vina"
params.parallel_jobs = 4

// Resource parameters
params.max_cpus = 16
params.max_memory = '64.GB'
params.max_time = '240.h'

// Environment and API configuration
params.fade_env_path = "${System.getenv('SCRATCH')}/conda-envs/fade"
params.gemini_api_key = System.getenv('GEMINI_API_KEY') ?: ""
params.gemini_model = System.getenv('GEMINI_MODEL') ?: "models/gemini-2.5-flash"

// Container/Module paths
params.alphafold_container = "/shared/container_repository/AlphaFold/alphafold3.sif"
params.schrodinger_module = "schrodinger/2024-4"

// Validation
if (!params.query) {
    error "Error: --query parameter is required"
}

if (!params.gemini_api_key) {
    log.warn "Warning: GEMINI_API_KEY not set. Some agents may fail."
}

// Print workflow header
log.info """
╔══════════════════════════════════════════════════════════════════╗
║                  F.A.D.E NEXTFLOW PIPELINE v2.0                  ║
║                    (Integrated Python Agents)                    ║
╚══════════════════════════════════════════════════════════════════╝
Query            :  ${params.query}
Output directory :  ${params.output_dir}
Max molecules    :  ${params.max_molecules}
Max iterations   :  ${params.max_iterations}
Docking method   :  ${params.docking_method}
Conda env        :  ${params.fade_env_path}
Work Directory   :  ${projectDir}
======================================================================
"""

// Main workflow
workflow FADE {
    
    // Create input channel from query
    query_ch = Channel.value(params.query)
    
    // Step 1: Target Selection
    // Uses TargetSelector agent to identify protein target from natural language query
    TARGET_SELECTION(query_ch)
    
    // Step 2: Structure Prediction/Retrieval
    // Uses StructurePredictor agent with AlphaFold3 or PDB retrieval
    STRUCTURE_PREDICTION(
        TARGET_SELECTION.out.target_info,
        TARGET_SELECTION.out.fasta,
        TARGET_SELECTION.out.pdb_id
    )
    
    // Step 3: Binding Site Analysis
    // Uses BindingSiteDetector to analyze pockets and binding sites
    BINDING_SITE_ANALYSIS(
        STRUCTURE_PREDICTION.out.structure,
        TARGET_SELECTION.out.target_info
    )
    
    // Step 4: Molecule Generation
    // Uses MoleculeGenerator agent with LLM guidance and RDKit validation
    iteration_ch = Channel.value(0)
    MOLECULE_GENERATION(
        TARGET_SELECTION.out.requirements,
        BINDING_SITE_ANALYSIS.out.binding_sites,
        TARGET_SELECTION.out.target_info,
        iteration_ch
    )
    
    // Step 5: Molecular Docking
    // Performs high-throughput virtual screening with AutoDock Vina or Glide
    DOCKING(
        MOLECULE_GENERATION.out.molecules,
        STRUCTURE_PREDICTION.out.prepared_receptor,
        BINDING_SITE_ANALYSIS.out.binding_sites
    )
    
    // Step 6: Lead Optimization
    // Uses Refiner agent concepts to optimize top hits
    LEAD_OPTIMIZATION(
        DOCKING.out.top_hits,
        STRUCTURE_PREDICTION.out.structure
    )
    
    // Step 7: Generate Final Report
    // Creates comprehensive natural language and technical reports
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
- 02_structure_prediction/structure.pdb
- 03_binding_site_analysis/binding_sites.json
- 04_molecule_generation/molecules.sdf
- 05_docking/docking_results.json
- 06_lead_optimization/optimized_leads.json
- 07_final_report/analysis_summary.md

To view results: cat ${params.output_dir}/07_final_report/analysis_summary.md
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
3. Ensure all required modules are available on the cluster
4. Check individual process logs in work/ directory
5. Try running with --with-trace for detailed execution tracking

Command to resume from failure:
./run_nextflow_integrated.sh --resume "${params.query}"
"""
}
