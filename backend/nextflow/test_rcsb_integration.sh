#!/bin/bash

echo "🧪 TESTING RCSB-INTEGRATED STRUCTURE PREDICTION"
echo "================================================"

# Navigate to test directory
cd /scratch/rajagopalmohanraj.n/F.A.D.E/nextflow/bin

# Create test input files (simulating target_selection output)
mkdir -p rcsb_integration_test
cd rcsb_integration_test

# Create target_info.json (Human KRAS - should be in PDB)
cat > target_info.json << 'EOF'
{
  "target": "KRAS",
  "uniprot_id": "P01116",
  "pdb_id": "none",
  "description": "GTPase KRas",
  "mutation": "G12D",
  "mutations": [
    {
      "original_residue": "G", 
      "position": 12,
      "mutated_residue": "D"
    }
  ],
  "organism": "Homo sapiens",
  "success": true
}
EOF

# Create protein.fasta (Human KRAS - more likely to be in PDB)
cat > protein.fasta << 'EOF'
>P01116|KRAS|GTPase KRas [Homo sapiens]
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMSCKCVLS
EOF

echo "📄 Created test input files:"
echo "   - target_info.json ($(stat -c%s target_info.json 2>/dev/null || echo "unknown") bytes)"
echo "   - protein.fasta ($(stat -c%s protein.fasta 2>/dev/null || echo "unknown") bytes)"
echo ""

echo "🔍 Test 1: RCSB Integration Test"
echo "================================"
echo "Running with RCSB lookup enabled (should find KRAS structure)..."

python ../run_structure_predictor.py \
    --target-info target_info.json \
    --fasta-file protein.fasta \
    --output-dir structure_test \
    --use-alphafold

echo ""
echo "📋 Test 1 Results:"
echo "=================="
ls -la *.pdb *.json 2>/dev/null || echo "No output files found"

echo ""
echo "📄 structure_analysis.json:"
if [ -f "structure_analysis.json" ]; then
    cat structure_analysis.json | python -m json.tool 2>/dev/null || cat structure_analysis.json
else
    echo "File not found"
fi

echo ""
echo "📄 structure.pdb header:"
if [ -f "structure.pdb" ]; then
    head -10 structure.pdb
else
    echo "File not found"
fi

echo ""
echo "✅ RCSB Integration Test Completed"
echo ""
echo "Expected results:"
echo "- Should find experimental KRAS structure from PDB"
echo "- Structure retrieval should be fast (under 5 minutes)"
echo "- Should skip AlphaFold3 completely"
echo "- structure_analysis.json should show method: 'experimental_rcsb'"
