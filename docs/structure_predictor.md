# Structure Predictor Agent

## Overview

The Structure Predictor Agent is a key component of the F.A.D.E (Fully Agentic Drug Engine) pipeline. This agent handles the generation, validation, and analysis of 3D protein structures from sequence data, and identifies potential binding sites for drug targeting.

## Features

- **AlphaFold3 Integration**: Interfaces with AlphaFold3 for state-of-the-art protein structure prediction
- **Agentic Error Recovery**: Automatically handles errors and retries operations
- **Structure Validation**: Validates predicted structures using multiple metrics
- **Binding Site Detection**: Identifies potential binding sites for drug targeting
- **Docking Preparation**: Prepares structures for molecular docking

## Components

### Structure Predictor Agent

The main agent class that orchestrates the structure prediction process. It:
1. Takes protein sequences from the Target Selector
2. Runs structure prediction using AlphaFold3
3. Validates the predicted structures
4. Identifies potential binding sites
5. Prepares the structures for docking

### PDB Processor

Handles parsing, analysis, and preparation of PDB files:
- Extracts key information from PDB files
- Prepares structures for docking
- Identifies binding sites based on geometric analysis

### Structure Validator

Validates protein structures and provides quality metrics:
- Checks for missing residues
- Evaluates geometry
- Extracts confidence scores (pLDDT, pTM)
- Provides human-readable interpretations of validation scores

### Binding Site Detector

Identifies potential binding sites in protein structures:
- Processes known binding sites
- Detects new binding sites using geometric analysis
- Uses LLM to analyze and rank binding sites

## Usage Example

```python
from agents.structure_predictor import StructurePredictor

# Initialize the agent
structure_predictor = StructurePredictor()

# Prepare input data
input_data = {
    "sequences": {
        "KRAS_G12D": {
            "sequence": "MTEYKLVVVGADGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM",
            "binding_sites": [
                {
                    "name": "GTP Binding Pocket",
                    "type": "nucleotide",
                    "residues": [10, 11, 12, 13, 14, 15, 16, 17]
                }
            ],
            "mutations": [
                {
                    "original_residue": "G",
                    "position": 12,
                    "mutated_residue": "D"
                }
            ]
        }
    },
    "job_configs": {
        "KRAS_G12D_alphafold_job": "path/to/job_script.sh"
    }
}

# Run the structure prediction process
results = structure_predictor.process(input_data)

# Access results
for target_name, structure_info in results["structures"].items():
    print(f"Structure for {target_name}:")
    print(f"  PDB file: {structure_info['pdb_file']}")
    print(f"  Confidence: {structure_info['confidence_scores']['overall']}")

for target_name, binding_sites in results["binding_sites"].items():
    print(f"Binding sites for {target_name}:")
    for i, site in enumerate(binding_sites):
        print(f"  Site {i+1}: {site['name']} (score: {site['score']})")
```

## Dependencies

- BioPython: For PDB file parsing and manipulation
- NumPy: For numerical operations
- SciPy: For spatial calculations
- Google Generative AI: For LLM-based analysis
- AlphaFold3: For structure prediction (via Singularity container)
- SLURM: For job scheduling on HPC clusters

## Future Improvements

- Integration with additional structure prediction tools (ESMFold, RoseTTAFold)
- More sophisticated binding site detection algorithms
- Integration with molecular dynamics simulations
- Support for multi-chain and protein-protein complexes
