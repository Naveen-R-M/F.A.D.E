/*
 * Molecule Generation Module
 * Uses MoleculeGenerator agent to create potential drug candidates
 */

process moleculeGeneration {
    tag "iteration_${iteration}"
    label 'process_high'
    publishDir "${params.output_dir}/04_molecule_generation", mode: 'copy'
    
    conda "${params.fade_env_path}"
    
    cpus 4
    memory '16.GB'
    time '2.h'
    
    input:
    path requirements
    path binding_sites
    path target_info
    val iteration
    
    output:
    path 'molecules.sdf', emit: molecules
    path 'molecules.json', emit: molecules_json
    path 'generation_stats.json', emit: stats
    
    script:
    def api_key = params.gemini_api_key ?: System.getenv('GEMINI_API_KEY') ?: ""
    def model = params.gemini_model ?: "models/gemini-2.5-flash"
    def max_molecules = params.max_molecules ?: 100
    """
    # Set up environment
    export PYTHONPATH="${projectDir}:\$PYTHONPATH"
    export GEMINI_API_KEY="${api_key}"
    
    # Load RDKit module if available
    module load RDKit/2023.09.5 || echo "RDKit module not available, using conda environment"
    
    # Run molecule generation
    run_molecule_generator.py \\
        --requirements ${requirements} \\
        --binding-sites ${binding_sites} \\
        --target-info ${target_info} \\
        --output-dir . \\
        --api-key "${api_key}" \\
        --model "${model}" \\
        --max-molecules ${max_molecules} \\
        --iteration ${iteration}
    
    # Verify outputs exist
    if [ ! -f "molecules.sdf" ]; then
        echo "# Molecule generation failed" > molecules.sdf
    fi
    
    if [ ! -f "molecules.json" ]; then
        echo '{"error": "Molecule generation failed", "molecules": []}' > molecules.json
    fi
    
    if [ ! -f "generation_stats.json" ]; then
        echo '{"total_generated": 0, "valid_molecules": 0}' > generation_stats.json
    fi
    """
    
    stub:
    """
    echo "MOL_001" > molecules.sdf
    echo "  F.A.D.E Generated" >> molecules.sdf
    echo "  SMILES: CCO" >> molecules.sdf
    echo '\$\$\$\$' >> molecules.sdf
    echo '[{"id": "MOL_001", "smiles": "CCO"}]' > molecules.json
    echo '{"total_generated": 1, "valid_molecules": 1}' > generation_stats.json
    """
}

workflow MOLECULE_GENERATION {
    take:
    requirements
    binding_sites
    target_info
    iteration
    
    main:
    moleculeGeneration(requirements, binding_sites, target_info, iteration)
    
    emit:
    molecules = moleculeGeneration.out.molecules
    molecules_json = moleculeGeneration.out.molecules_json
    stats = moleculeGeneration.out.stats
}
