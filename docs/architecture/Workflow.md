# F.A.D.E: Fully Agentic Drug Engine

## Project Overview

F.A.D.E is an agentic AI pipeline for discovering drug-like molecules targeting specific proteins. The system transforms natural language queries into a computational drug discovery workflow, using multiple specialized agents to handle different aspects of the discovery process.

## System Components

- **Target Selector Agent**: Processes natural language to identify protein targets
- **Structure Predictor Agent**: Generates 3D protein structures
- **Molecule Generator Agent**: Creates potential drug candidates
- **Evaluator Agent**: Assesses drug-like properties
- **Docking Agent**: Simulates protein-ligand binding
- **Refiner Agent**: Improves molecules based on evaluation results
- **Memory Manager**: Tracks history of attempts and results

## Detailed Workflow Diagrams

### 1. Initial Query Processing Phase

```
User Query (Natural Language)
        │
        ▼
┌─────────────────┐
│       LLM       │  Parse query and extract parameters
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ Target Selector │  Identify protein target and requirements
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│   UniProt API   │  Retrieve protein sequence
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ Config Generator│  Create config files and job scripts
└────────┬────────┘
        │
        ▼
Job Configuration Files
```

### 2. Structure Prediction Phase

```
Protein Sequence (FASTA)
        │
        ▼
┌─────────────────┐
│  SLURM Scheduler│  Queue GPU job
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│   AlphaFold3    │  3D structure prediction
│   Container     │  (Running on GPU node)
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ Structure       │  Analyze quality, identify binding sites
│ Analysis (PyMOL)│
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│  PDB Processing │  Prepare structure for docking
│  (BioPython)    │
└────────┬────────┘
        │
        ▼
Protein Structure (PDB format)
```

### 3. Molecule Generation Phase

```
Target Information + Binding Site Data
        │
        ▼
┌─────────────────┐
│       LLM       │  Generate molecule ideas based on context
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ SMILES Generator│  Create valid SMILES strings
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│     RDKit       │  Validate molecules and generate 3D structures
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ Property        │  Calculate molecular properties
│ Calculator      │  (LogP, MW, TPSA, etc.)
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ Lipinski Filter │  Apply Rule of Five and other filters
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ ML-based        │  Predict BBB permeability, toxicity
│ Property Models │
└────────┬────────┘
        │
        ▼
Candidate Molecules (SDF/MOL2 format)
```

### 4. Docking Simulation Phase

```
Protein Structure + Candidate Molecules
        │
        ▼
┌─────────────────┐
│ Receptor Prep   │  Prepare protein for docking
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ Ligand Prep     │  Prepare molecules for docking
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│  SLURM Scheduler│  Queue docking job
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ Schrodinger     │  Perform molecular docking
│ Glide           │  (Running on CPU nodes)
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ Pose Analysis   │  Analyze binding poses
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ Interaction     │  Map protein-ligand interactions
│ Mapping         │
└────────┬────────┘
        │
        ▼
Docking Results (Scores + Poses)
```

### 5. Refinement Loop Phase

```
Docking Results + Property Data
        │
        ▼
┌─────────────────┐
│ Candidate       │  Rank and select molecules for refinement
│ Ranking         │
└────────┬────────┘
        │
        ▼
        ┌───────────────────────────────┐
        │                               │
        ▼                               │
┌─────────────────┐                     │
│ LLM-guided      │  Suggest molecular modifications  │
│ Refinement      │                     │
└────────┬────────┘                     │
        │                               │
        ▼                               │
┌─────────────────┐                     │
│ Structure       │  Implement suggested changes     │
│ Modification    │                     │
└────────┬────────┘                     │
        │                               │
        ▼                               │
┌─────────────────┐                     │
│ RDKit Property  │  Recalculate properties          │
│ Calculation     │                     │
└────────┬────────┘                     │
        │                               │
        ▼                               │
┌─────────────────┐                     │
│ Re-docking      │  Dock refined molecules          │
│ (Glide)         │                     │
└────────┬────────┘                     │
        │                               │
        ▼                               │
┌─────────────────┐                     │
│ Improvement     │  Evaluate if improvements made   │
│ Assessment      │                     │
└────────┬────────┘                     │
        │                               │
        └───────────────────────────────┘
        │          (Loop 3-5 times)
        ▼
Refined Candidates (Final Selection)
```

### 6. Result Analysis & Presentation Phase

```
Final Candidate Molecules + All Data
        │
        ▼
┌─────────────────┐
│ Data Aggregation│  Compile all results
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ Comparative     │  Compare candidates against each other
│ Analysis        │
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ Technical Data  │  Structured data of all findings
│ Compilation     │
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ JSON/CSV        │  Format results in structured format
│ Formatter       │
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│      LLM        │  Convert technical data to natural language
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ Result          │  Format and organize the final report
│ Presentation    │
└────────┬────────┘
        │
        ▼
Natural Language Results + Supporting Data
```

### 7. End-to-End Complete Flow

```
User Query (Natural Language)
        │
        ▼
┌─────────────────┐
│ Target Selector │  LLM parses query to extract parameters
│     Agent       │
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│   Structure     │  AlphaFold3 generates 3D protein structure
│ Predictor Agent │
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│   Molecule      │  LLM + RDKit generate initial molecules
│ Generator Agent │
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│   Evaluator     │  RDKit calculates properties and filters
│     Agent       │
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│    Docking      │  Schrodinger Glide performs docking
│     Agent       │
└────────┬────────┘
        │
        │◄───────────────────┐
        ▼                    │
┌─────────────────┐          │
│    Refiner      │  LLM suggests improvements
│     Agent       │          │
└────────┬────────┘          │
        │                    │
        ▼                    │
┌─────────────────┐          │
│  Re-evaluation  │  Repeat evaluation and docking
│     Cycle       │          │
└────────┬────────┘          │
        │                    │
        └────────────────────┘
        │     (3-5 cycles)
        ▼
┌─────────────────┐
│ Result Analysis │  Compile and analyze all data
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ Technical Data  │  JSON/CSV structured data
└────────┬────────┘
        │
        ▼
┌─────────────────┐
│ LLM Translation │  Convert to natural language
└────────┬────────┘
        │
        ▼
User-Friendly Results (Natural Language)
```

### 8. Data Transformation Flow

```
Natural Language Query
        │
        ▼
Target Parameters (Structured Data)
        │
        ▼
Protein Sequence (FASTA)
        │
        ▼
3D Protein Structure (PDB)
        │
        ▼
Molecule Ideas (Conceptual)
        │
        ▼
SMILES Strings (Chemical Notation)
        │
        ▼
3D Molecules (SDF/MOL2)
        │
        ▼
Property Data (Numerical Values)
        │
        ▼
Docking Results (Scores + 3D Coordinates)
        │
        ▼
Binding Interactions (Structural Data)
        │
        ▼
Refinement Suggestions (Structured Modifications)
        │
        ▼
Final Candidates (Complete Data Package)
        │
        ▼
Technical Results (JSON/CSV)
        │
        ▼
Natural Language Explanation
```

## Implementation Timeline

- **Week 1-2**: Environment setup, initial agent implementation
- **Week 3-4**: Complete individual agent development
- **Week 5-6**: Agent integration, workflow orchestration
- **Week 7-8**: Testing, optimization, and documentation

## Sample Input and Output

### Sample Input (Natural Language)

```
I'm looking for potential drug candidates that could target the KRAS G12D mutant protein, which is implicated in pancreatic cancer. The molecules should be able to cross the blood-brain barrier, have low toxicity, and follow Lipinski's rule of five for drug-likeness. I'm particularly interested in molecules that might bind to the GTP binding pocket of KRAS to inhibit its activity. Please generate 3-5 promising candidates.
```

### Sample Output (Natural Language)

```
I've completed the drug discovery process for KRAS G12D mutant protein and found 3 promising candidates that meet your criteria. Here's what I discovered:

CANDIDATE 1: "Compound KR-371"
This molecule has excellent properties for a potential drug. It has a molecular weight of 312.3 daltons and a LogP value of 2.1, making it well within Lipinski's rules. Our analysis shows it should readily cross the blood-brain barrier and demonstrates low predicted toxicity in preliminary models.

Most importantly, it showed strong binding to the GTP pocket of KRAS G12D with a binding score of -9.1 kcal/mol. The molecule forms key hydrogen bonds with residues Asp57 and Ser39, while also making a π-stacking interaction with Tyr32 in the binding pocket.

CANDIDATE 2: "Compound KR-225"
This alternative candidate has a slightly lower molecular weight of 290.2 daltons with a LogP of 1.8. It also crosses the blood-brain barrier effectively and shows low toxicity profiles.

Its binding to KRAS G12D is strong (-8.7 kcal/mol) but approaches the target differently, interacting primarily with the switch II region near the GTP binding pocket. This could provide an alternative mechanism of inhibition.

CANDIDATE 3: "Compound KR-493"
Our third candidate is structurally distinct from the first two, offering a different scaffold for optimization. With a molecular weight of 328.4 and LogP of 2.6, it still maintains drug-like properties while having the potential to be further modified.

Its binding score of -8.3 kcal/mol is slightly lower than the other candidates, but it forms unique interactions with Gly12D itself - the mutated residue of interest - which could provide specificity for the mutant over the wild-type protein.

All three candidates were generated through an iterative process that evaluated 43 total molecules across 7 design cycles. The 3D structure of KRAS G12D was predicted using AlphaFold3 with a high confidence score (pLDDT > 90), and docking was performed using Schrodinger's Glide software.

Would you like me to provide the detailed molecular structures for these compounds, or would you prefer to explore any specific aspect of these candidates in more depth?
```

## Technical Resources Used

### Cluster Resources
- Northeastern University HPC Cluster
- GPU Partition (V100, A100, H200 GPUs)
- SLURM Job Scheduler
- SCRATCH Space for data storage

### Software Tools
- AlphaFold3 (Structure Prediction)
- ESMFold (Alternative Structure Prediction)
- RDKit (Molecule Handling and Property Calculation)
- Schrodinger Glide (Molecular Docking)
- AutoDock Vina (Alternative Docking)
- LangGraph (Agent Orchestration)
- Claude/GPT API (LLM Integration)

### Environment Configuration
- Python 3.9+
- CUDA 12.3.0
- Miniconda3/24.11.1

## Expected Runtime

- **Initial Processing**: 5-10 minutes
- **Structure Prediction**: 1-4 hours (+ queue time)
- **Initial Molecule Generation**: 45-75 minutes
- **Docking**: 1-2 hours (+ queue time)
- **Refinement Cycles**: 2-6 hours (3-5 cycles)
- **Result Compilation & Delivery**: 30-60 minutes

**Total Expected Time**: 5-14 hours (Typically 6-8 hours for standard runs)
