/*
 * Molecular Docking Module
 * Performs virtual screening
 */

process prepareLigands {
    tag "ligand_prep"
    label 'process_medium'
    publishDir "${params.output_dir}/05_docking/prepared", mode: 'copy'
    
    input:
    path molecules
    
    output:
    path "ligands_prepared.pdbqt", emit: ligands
    
    script:
    """
    # Create dummy prepared ligands
    echo "REMARK Prepared ligands" > ligands_prepared.pdbqt
    echo "ATOM      1  C   LIG A   1       0.000   0.000   0.000" >> ligands_prepared.pdbqt
    """
    
    stub:
    """
    touch ligands_prepared.pdbqt
    """
}

process runDocking {
    tag "docking"
    label 'docking'
    publishDir "${params.output_dir}/05_docking/results", mode: 'copy'
    
    input:
    path receptor
    path ligands
    path binding_sites
    
    output:
    path "docking_results.csv", emit: results
    
    script:
    """
    # Create dummy docking results
    echo "ligand,docking_score" > docking_results.csv
    echo "molecule_1,-7.8" >> docking_results.csv
    echo "molecule_2,-6.5" >> docking_results.csv
    """
    
    stub:
    """
    echo "ligand,score" > docking_results.csv
    """
}

process aggregateDockingResults {
    tag "aggregation"
    label 'process_low'
    publishDir "${params.output_dir}/05_docking", mode: 'copy'
    
    input:
    path results
    
    output:
    path "all_docking_results.csv", emit: all_results
    path "top_hits.csv", emit: top_hits
    
    script:
    """
    cp ${results} all_docking_results.csv
    head -n 2 all_docking_results.csv > top_hits.csv
    """
    
    stub:
    """
    echo "ligand,score" > all_docking_results.csv
    echo "ligand,score" > top_hits.csv
    """
}

workflow DOCKING {
    take:
    molecules
    receptor
    binding_sites
    
    main:
    prepareLigands(molecules)
    runDocking(receptor, prepareLigands.out.ligands, binding_sites)
    aggregateDockingResults(runDocking.out.results)
    
    emit:
    all_results = aggregateDockingResults.out.all_results
    top_hits = aggregateDockingResults.out.top_hits
}
