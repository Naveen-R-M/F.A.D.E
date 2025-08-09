/*
 * Target Selection Module
 * Processes natural language queries to identify protein targets
 */

process targetSelection {
    tag "$query"
    label 'process_medium'
    publishDir "${params.output_dir}/01_target_selection", mode: 'copy'
    
    container 'python:3.9'
    
    input:
    val query
    
    output:
    path 'target_info.json', emit: target_info
    path 'protein.fasta', emit: fasta
    val 'none', emit: pdb_id
    path 'requirements.json', emit: requirements
    
    script:
    """
    # Create stub outputs for testing (since Python agents aren't set up yet)
    cat > target_info.json << 'JSONEOF'
    {
        "target": "KRAS",
        "uniprot_id": "P01116",
        "pdb_id": "none",
        "description": "GTPase KRas"
    }
    JSONEOF
    
    cat > protein.fasta << 'FASTAEOF'
    >KRAS|P01116|GTPase KRas
    MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMSCKCVLS
    FASTAEOF
    
    cat > requirements.json << 'JSONEOF'
    {
        "binding_affinity": "< -8 kcal/mol",
        "selectivity": "high",
        "toxicity": "low",
        "bbb_permeability": "good"
    }
    JSONEOF
    
    echo "none" > pdb_id.txt
    """
    
    stub:
    """
    echo '{"target": "KRAS", "uniprot_id": "P01116", "pdb_id": "none"}' > target_info.json
    echo '>KRAS|P01116|GTPase KRas' > protein.fasta
    echo 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMSCKCVLS' >> protein.fasta
    echo '{"binding_affinity": "< -8 kcal/mol"}' > requirements.json
    echo "none" > pdb_id.txt
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
