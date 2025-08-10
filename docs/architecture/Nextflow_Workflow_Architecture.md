# F.A.D.E Nextflow Workflow Architecture

## Table of Contents
1. [Introduction to Nextflow](#introduction-to-nextflow)
2. [Benefits for F.A.D.E Project](#benefits-for-fade-project)
3. [High-Level Architecture](#high-level-architecture)
4. [Detailed Workflow Architecture](#detailed-workflow-architecture)
5. [Data Flow Architecture](#data-flow-architecture)
6. [Resource Management](#resource-management)
7. [Error Handling & Recovery](#error-handling--recovery)
8. [Monitoring & Reporting](#monitoring--reporting)
9. [File Organization](#file-organization)
10. [Integration Strategy](#integration-strategy)
11. [Example Implementation](#example-implementation)

## Introduction to Nextflow

Nextflow is a powerful workflow management system and domain-specific language (DSL) designed for writing and executing data-intensive computational pipelines, particularly in bioinformatics and scientific computing.

### Key Features

#### 1. **Workflow Orchestration**
- Allows you to chain together different processes/tools into complex pipelines
- Automatically manages task dependencies and execution order
- Handles data flow between pipeline steps

#### 2. **Parallel & Distributed Computing**
- Automatically parallelizes independent tasks
- Can scale from laptops to HPC clusters and cloud platforms
- Supports various executors: local, SLURM, PBS, AWS Batch, Google Cloud, Kubernetes

#### 3. **Portability & Reproducibility**
- Integrates with container technologies (Docker, Singularity, Podman)
- Ensures pipelines run consistently across different computing environments
- Built-in support for Conda environments

#### 4. **Resume Capability**
- Can resume failed pipelines from the last successful step
- Caches intermediate results to avoid re-computation
- Tracks execution history and provenance

#### 5. **Key Concepts**
- **Processes**: Individual computational steps (like running a program or script)
- **Channels**: Data flow connections between processes
- **Workflows**: Complete pipelines combining multiple processes
- **Operators**: Methods to manipulate and transform data channels

### Common Use Cases in Drug Discovery
- RNA-seq analysis pipelines
- Genome assembly and annotation
- Variant calling workflows
- Complex multi-step simulations
- High-throughput virtual screening
- Lead optimization workflows

## Benefits for F.A.D.E Project

### 1. **Workflow Orchestration & Agent Coordination**
Your F.A.D.E project currently has multiple agents that need to work in sequence:
- Target Selector Agent → Structure Predictor → Molecule Generator → Evaluator → Docking → Refiner

**Nextflow would help by:**
- Defining each agent as a process with clear inputs/outputs
- Automatically managing dependencies between agents
- Handling data flow between stages without manual file management
- Providing a visual workflow that clearly shows the pipeline structure

### 2. **HPC Integration & Resource Management**
Your project already uses SLURM (`run_fade_slurm.sh`) and needs different resources for different stages:
- AlphaFold3 needs GPU nodes
- Docking might need high-memory nodes
- Molecule generation needs CPU resources

**Nextflow would:**
- Automatically submit jobs to SLURM with appropriate resource requirements
- Handle queue management and job scheduling
- Retry failed jobs automatically
- Scale resources based on workload

### 3. **Parallel Processing & Scalability**
Your pipeline could process multiple molecules or targets in parallel:

**Nextflow enables:**
- Parallel execution of independent tasks (e.g., evaluating 100 molecules simultaneously)
- Dynamic parallelization based on available resources
- Efficient use of HPC resources by running multiple docking simulations in parallel

### 4. **Resume Capability & Fault Tolerance**
Drug discovery pipelines can run for days and failures are common:

**Nextflow provides:**
- Automatic checkpointing of completed stages
- Resume from last successful step if pipeline fails
- No need to restart from scratch after failures
- Caching of intermediate results

### 5. **Containerization & Reproducibility**
Your project uses AlphaFold3 container and various Python environments:

**Nextflow offers:**
- Native support for Singularity/Docker containers
- Each process can use different containers/environments
- Ensures reproducible results across different HPC systems
- Simplifies dependency management

## High-Level Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                     NEXTFLOW ORCHESTRATION LAYER                 │
│  ┌─────────────────────────────────────────────────────────┐    │
│  │  Main Workflow Controller (main.nf)                     │    │
│  │  • Parameter validation                                 │    │
│  │  • Workflow branching logic                            │    │
│  │  • Resource allocation strategy                        │    │
│  └─────────────────────────────────────────────────────────┘    │
│                                                                   │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐         │
│  │ Config       │  │ Executors    │  │ Containers   │         │
│  │ • params     │  │ • SLURM      │  │ • AlphaFold3 │         │
│  │ • resources  │  │ • Local      │  │ • RDKit      │         │
│  │ • profiles   │  │ • AWS Batch  │  │ • Schrodinger│         │
│  └──────────────┘  └──────────────┘  └──────────────┘         │
└─────────────────────────────────────────────────────────────────┘
                                ↓
┌─────────────────────────────────────────────────────────────────┐
│                        WORKFLOW PROCESSES                        │
└─────────────────────────────────────────────────────────────────┘
```

## Detailed Workflow Architecture

```
USER QUERY
    │
    ▼
┌───────────────────────────────────────────────────────────────┐
│                    INITIALIZATION PHASE                        │
├───────────────────────────────────────────────────────────────┤
│  Process: validateInput                                        │
│  • Parse natural language query                                │
│  • Validate required parameters                                │
│  • Set workflow configuration                                  │
│  Resources: 1 CPU, 2GB RAM                                     │
│  Container: python:3.9                                         │
└────────────────┬──────────────────────────────────────────────┘
                 │
                 ▼
┌───────────────────────────────────────────────────────────────┐
│                  TARGET IDENTIFICATION PHASE                   │
├───────────────────────────────────────────────────────────────┤
│  Process: identifyTarget                                       │
│  • LLM parsing of drug requirements                           │
│  • UniProt API queries                                        │
│  • Sequence retrieval                                         │
│  Resources: 2 CPUs, 4GB RAM                                   │
│  Container: fade-base                                         │
│  Output: Channel<target_info.json, protein.fasta>            │
└────────────────┬──────────────────────────────────────────────┘
                 │
                 ├─────────────[Fork: Structure Available?]
                 │                          │
                 ▼                          ▼
        [NO: Need Prediction]      [YES: Fetch from PDB]
                 │                          │
┌────────────────▼──────────────┐  ┌───────▼────────────────┐
│   STRUCTURE PREDICTION PHASE  │  │  STRUCTURE RETRIEVAL    │
├────────────────────────────────┤  ├─────────────────────────┤
│  Process: predictStructure     │  │  Process: fetchPDB      │
│  • AlphaFold3 execution       │  │  • Download from RCSB   │
│  • Quality assessment         │  │  • Structure validation │
│  Resources: 1 GPU, 32GB RAM   │  │  Resources: 1 CPU, 2GB  │
│  Container: alphafold3.sif    │  │  Container: pymol       │
│  Queue: gpu                   │  │  Queue: shared          │
│  Time: 8 hours               │  │  Time: 10 minutes      │
└────────────────┬──────────────┘  └───────┬────────────────┘
                 │                          │
                 └──────────┬───────────────┘
                           │
                           ▼
┌───────────────────────────────────────────────────────────────┐
│                 BINDING SITE ANALYSIS PHASE                    │
├───────────────────────────────────────────────────────────────┤
│  Process: analyzeBindingSites                                  │
│  • Cavity detection (CASTp/P2Rank)                            │
│  • Active site identification                                  │
│  • Pocket characterization                                     │
│  Resources: 4 CPUs, 8GB RAM                                   │
│  Container: structure-tools                                    │
│  Output: Channel<binding_sites.json, prepared_receptor.pdb>   │
└────────────────┬──────────────────────────────────────────────┘
                 │
                 ▼
┌───────────────────────────────────────────────────────────────┐
│              MOLECULE GENERATION PHASE (PARALLEL)              │
├───────────────────────────────────────────────────────────────┤
│  Process: generateMoleculesBatch                               │
│  • Splits into N parallel jobs                                 │
│  • Each job generates M molecules                              │
│  ┌─────────────────────────────────────────────────────┐     │
│  │  Subworkflow: moleculeGeneration                     │     │
│  │  ┌─────────────────┐  ┌─────────────────┐          │     │
│  │  │ LLM Generation  │→│ SMILES Validation│          │     │
│  │  └─────────────────┘  └─────────────────┘          │     │
│  │           ↓                    ↓                    │     │
│  │  ┌─────────────────┐  ┌─────────────────┐          │     │
│  │  │ 3D Conformer    │→│ Property Filter  │          │     │
│  │  └─────────────────┘  └─────────────────┘          │     │
│  └─────────────────────────────────────────────────────┘     │
│  Resources: 4 CPUs, 8GB RAM per job                          │
│  Container: rdkit-env                                         │
│  Parallelism: 10-100 concurrent jobs                         │
│  Output: Channel<molecules_*.sdf>.flatten()                  │
└────────────────┬──────────────────────────────────────────────┘
                 │
                 ▼
┌───────────────────────────────────────────────────────────────┐
│           INITIAL SCREENING PHASE (PARALLEL)                   │
├───────────────────────────────────────────────────────────────┤
│  Process: rapidScreen                                          │
│  • Pharmacophore matching                                      │
│  • Shape-based screening                                       │
│  • Quick ADMET predictions                                     │
│  Resources: 2 CPUs, 4GB RAM per molecule                      │
│  Container: screening-tools                                    │
│  Parallelism: 100 concurrent                                  │
│  Output: Channel<filtered_molecules.sdf>                      │
└────────────────┬──────────────────────────────────────────────┘
                 │
                 ▼
┌───────────────────────────────────────────────────────────────┐
│              DOCKING PHASE (PARALLEL)                          │
├───────────────────────────────────────────────────────────────┤
│  Process: performDocking                                       │
│  ┌─────────────────────────────────────────────────────┐     │
│  │  Dynamic Resource Allocation:                        │     │
│  │  • If molecules < 100: Glide SP (4 CPU, 8GB)       │     │
│  │  • If molecules < 1000: Glide HTVS (2 CPU, 4GB)    │     │
│  │  • If molecules > 1000: Vina (1 CPU, 2GB)          │     │
│  └─────────────────────────────────────────────────────┘     │
│  Container: schrodinger-suite                                 │
│  Queue: docking-priority                                      │
│  Parallelism: 50 concurrent                                   │
│  Output: Channel<docking_results_*.csv>                      │
└────────────────┬──────────────────────────────────────────────┘
                 │
                 ▼
┌───────────────────────────────────────────────────────────────┐
│                  SCORING & RANKING PHASE                       │
├───────────────────────────────────────────────────────────────┤
│  Process: scoreAndRank                                         │
│  • Aggregate docking scores                                    │
│  • MM-GBSA rescoring (top 10%)                                │
│  • Consensus scoring                                           │
│  Resources: 8 CPUs, 16GB RAM                                  │
│  Container: scoring-tools                                      │
│  Output: Channel<ranked_molecules.csv, top_hits.sdf>         │
└────────────────┬──────────────────────────────────────────────┘
                 │
                 ├─────────[Branch: Satisfactory Results?]
                 │                     │
                 ▼                     ▼
           [NO: Refine]          [YES: Finalize]
                 │                     │
┌────────────────▼──────────────┐     │
│     REFINEMENT PHASE          │     │
├────────────────────────────────┤     │
│  Process: refineMolecules      │     │
│  • Analyze failed molecules   │     │
│  • LLM-guided modifications   │     │
│  • Fragment replacement       │     │
│  • Scaffold hopping          │     │
│  Resources: 4 CPUs, 8GB RAM   │     │
│  Container: fade-refiner      │     │
│  Output: Channel<refined.sdf> │     │
└────────────────┬──────────────┘     │
                 │                     │
                 └──[Loop back to]────→┘
                    Molecule Generation
                    (Max 3 iterations)
                           │
                           ▼
┌───────────────────────────────────────────────────────────────┐
│                    LEAD OPTIMIZATION PHASE                     │
├───────────────────────────────────────────────────────────────┤
│  Process: optimizeLeads                                        │
│  • Free energy perturbation (top 10)                         │
│  • ADMET optimization                                         │
│  • Synthetic accessibility                                    │
│  Resources: 8 CPUs, 32GB RAM per molecule                    │
│  Container: schrodinger-fep                                   │
│  Queue: high-memory                                           │
│  Time: 24 hours                                              │
└────────────────┬──────────────────────────────────────────────┘
                 │
                 ▼
┌───────────────────────────────────────────────────────────────┐
│                    REPORTING PHASE                             │
├───────────────────────────────────────────────────────────────┤
│  Process: generateReport                                       │
│  • Compile all results                                        │
│  • Generate visualizations                                    │
│  • Create PyMOL sessions                                      │
│  • Produce summary statistics                                 │
│  Resources: 2 CPUs, 4GB RAM                                  │
│  Container: reporting-tools                                   │
│  Output: Channel<final_report.html, results.xlsx>            │
└───────────────────────────────────────────────────────────────┘
```

## Data Flow Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                     NEXTFLOW CHANNELS                        │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  query_ch ────────→ target_ch ─────→ structure_ch          │
│                           ↓                ↓                │
│                     sequence_ch      binding_ch             │
│                                           ↓                 │
│                                     molecules_ch            │
│                                           ↓                 │
│                            ┌──────────────┴──────────┐      │
│                            ↓                         ↓      │
│                      screening_ch              docking_ch   │
│                            ↓                         ↓      │
│                            └──────────┬──────────────┘      │
│                                       ↓                     │
│                                  ranking_ch                 │
│                                       ↓                     │
│                                refinement_ch ←──────┐       │
│                                       ↓              │       │
│                                   leads_ch           │       │
│                                       ↓              │       │
│                                  report_ch           │       │
│                                                      │       │
│  Feedback Loop: ─────────────────────────────────────┘       │
│                                                              │
└──────────────────────────────────────────────────────────────┘
```

### Channel Definitions

- **query_ch**: Initial user query in natural language
- **target_ch**: Identified protein target information
- **sequence_ch**: Protein sequence data (FASTA format)
- **structure_ch**: 3D protein structure (PDB format)
- **binding_ch**: Binding site information and prepared receptor
- **molecules_ch**: Generated molecule candidates (SDF format)
- **screening_ch**: Molecules that passed initial screening
- **docking_ch**: Docking results with scores
- **ranking_ch**: Ranked molecules based on all metrics
- **refinement_ch**: Refined molecules for next iteration
- **leads_ch**: Final optimized lead compounds
- **report_ch**: Comprehensive results and visualizations

## Resource Management

```
┌─────────────────────────────────────────────────────────────┐
│                 NEXTFLOW RESOURCE MANAGER                    │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  ┌──────────────────────────────────────────┐               │
│  │         Queue Assignment Logic           │               │
│  ├──────────────────────────────────────────┤               │
│  │  • GPU Tasks → gpu-queue (limited)       │               │
│  │  • Docking → docking-queue (dedicated)   │               │
│  │  • Small tasks → shared-queue            │               │
│  │  • Memory intensive → highmem-queue      │               │
│  └──────────────────────────────────────────┘               │
│                                                              │
│  ┌──────────────────────────────────────────┐               │
│  │      Dynamic Scaling Configuration       │               │
│  ├──────────────────────────────────────────┤               │
│  │  molecules_count = molecules_ch.count()  │               │
│  │                                          │               │
│  │  if (molecules_count < 100) {            │               │
│  │    cpus = 4; memory = '8GB'             │               │
│  │    executor.queueSize = 10              │               │
│  │  } else if (molecules_count < 1000) {    │               │
│  │    cpus = 2; memory = '4GB'             │               │
│  │    executor.queueSize = 50              │               │
│  │  } else {                                │               │
│  │    cpus = 1; memory = '2GB'             │               │
│  │    executor.queueSize = 100             │               │
│  │  }                                       │               │
│  └──────────────────────────────────────────┘               │
│                                                              │
└──────────────────────────────────────────────────────────────┘
```

### Resource Profiles by Process

| Process | CPUs | Memory | GPU | Queue | Time Limit |
|---------|------|--------|-----|-------|------------|
| validateInput | 1 | 2GB | - | shared | 10m |
| identifyTarget | 2 | 4GB | - | shared | 30m |
| predictStructure | 4 | 32GB | 1 | gpu | 8h |
| fetchPDB | 1 | 2GB | - | shared | 10m |
| analyzeBindingSites | 4 | 8GB | - | shared | 1h |
| generateMolecules | 4 | 8GB | - | shared | 2h |
| rapidScreen | 2 | 4GB | - | shared | 30m |
| performDocking | 4 | 8GB | - | docking | 4h |
| scoreAndRank | 8 | 16GB | - | shared | 1h |
| refineMolecules | 4 | 8GB | - | shared | 2h |
| optimizeLeads | 8 | 32GB | - | highmem | 24h |
| generateReport | 2 | 4GB | - | shared | 30m |

## Error Handling & Recovery

```
┌─────────────────────────────────────────────────────────────┐
│              NEXTFLOW ERROR HANDLING LAYER                   │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  ┌──────────────────────────────────────────┐               │
│  │         Process-Level Retry Logic         │               │
│  ├──────────────────────────────────────────┤               │
│  │  errorStrategy = 'retry'                  │               │
│  │  maxRetries = 3                          │               │
│  │  errorBackOff = 'exponential'            │               │
│  │                                          │               │
│  │  Retry Conditions:                       │               │
│  │  • Exit 137: Memory exceeded → +50% RAM  │               │
│  │  • Exit 140: Time limit → +2 hours      │               │
│  │  • Exit 1: General failure → retry      │               │
│  └──────────────────────────────────────────┘               │
│                                                              │
│  ┌──────────────────────────────────────────┐               │
│  │       Checkpoint & Resume System          │               │
│  ├──────────────────────────────────────────┤               │
│  │  Work Directory Structure:                │               │
│  │  work/                                    │               │
│  │  ├── 2a/3f5b... (completed)             │               │
│  │  ├── 4c/8d9e... (completed)             │               │
│  │  ├── 6e/1a2c... (failed) ← retry here   │               │
│  │  └── cache.db                           │               │
│  └──────────────────────────────────────────┘               │
│                                                              │
└──────────────────────────────────────────────────────────────┘
```

### Error Recovery Strategies

1. **Automatic Retries**
   - Network failures: Immediate retry
   - Resource exhaustion: Retry with increased resources
   - Queue timeout: Retry with extended time limit

2. **Checkpoint Recovery**
   - All completed processes are cached
   - Resume command: `nextflow run main.nf -resume`
   - No recomputation of successful steps

3. **Failure Notifications**
   - Email alerts on critical failures
   - Slack integration for team notifications
   - Detailed error logs in `work/` directory

## Monitoring & Reporting

```
┌─────────────────────────────────────────────────────────────┐
│                  NEXTFLOW MONITORING LAYER                   │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  ┌──────────────────────────────────────────┐               │
│  │         Real-time Progress Tracking       │               │
│  ├──────────────────────────────────────────┤               │
│  │  Tower Dashboard Integration:             │               │
│  │  • Live process status                   │               │
│  │  • Resource utilization graphs           │               │
│  │  • Queue wait times                      │               │
│  │  • Success/failure rates                 │               │
│  └──────────────────────────────────────────┘               │
│                                                              │
│  ┌──────────────────────────────────────────┐               │
│  │           Execution Reports               │               │
│  ├──────────────────────────────────────────┤               │
│  │  reports/                                 │               │
│  │  ├── timeline.html (Gantt chart)         │               │
│  │  ├── report.html (summary)               │               │
│  │  ├── trace.txt (detailed metrics)        │               │
│  │  └── dag.dot (workflow graph)            │               │
│  └──────────────────────────────────────────┘               │
│                                                              │
└──────────────────────────────────────────────────────────────┘
```

### Monitoring Features

1. **Real-time Monitoring**
   - Console output with progress bars
   - Web interface via Nextflow Tower
   - Resource usage tracking
   - Queue status monitoring

2. **Performance Metrics**
   - CPU utilization per process
   - Memory consumption patterns
   - I/O statistics
   - Network usage for API calls

3. **Execution Reports**
   - HTML timeline showing process execution
   - Resource usage summaries
   - Bottleneck identification
   - Cost analysis for cloud deployments

## File Organization

```
F.A.D.E-nextflow/
├── main.nf                 # Main workflow
├── nextflow.config         # Configuration
├── modules/
│   ├── local/
│   │   ├── target_selector.nf
│   │   ├── structure_predictor.nf
│   │   ├── molecule_generator.nf
│   │   ├── evaluator.nf
│   │   ├── docking.nf
│   │   └── refiner.nf
│   └── nf-core/           # Community modules
│       ├── rdkit/
│       └── pymol/
├── subworkflows/
│   ├── structure_preparation.nf
│   ├── lead_optimization.nf
│   └── iterative_refinement.nf
├── conf/
│   ├── base.config        # Base configuration
│   ├── resources.config   # Resource definitions
│   ├── containers.config  # Container specs
│   └── profiles/
│       ├── northeastern.config
│       ├── aws.config
│       └── local.config
├── bin/                   # Helper scripts
│   └── (existing Python agents)
└── assets/
    └── schemas/
        └── parameters.json
```

### Directory Descriptions

- **main.nf**: Entry point defining the complete workflow
- **modules/**: Individual process definitions
- **subworkflows/**: Reusable workflow components
- **conf/**: Configuration files for different environments
- **bin/**: Executable scripts called by processes
- **assets/**: Static files and schemas

## Integration Strategy

### Phase 1: Hybrid Approach (Months 1-2)
```
┌──────────────────────────────────────────────────────────────┐
│              HYBRID ARCHITECTURE APPROACH                     │
├──────────────────────────────────────────────────────────────┤
│                                                               │
│  Nextflow Process ←→ Existing Python Agent                   │
│  ─────────────────────────────────────────                  │
│                                                               │
│  process runTargetSelector {                                 │
│      script:                                                 │
│      """                                                      │
│      python ${projectDir}/agents/target_selector.py \\      │
│          --query "${params.query}" \\                       │
│          --output target_info.json                          │
│      """                                                      │
│  }                                                           │
│                                                               │
│  Benefits:                                                   │
│  • Preserve existing agent logic                            │
│  • Gradual migration path                                   │
│  • Immediate workflow benefits                              │
│  • No code rewriting required initially                     │
│                                                               │
└──────────────────────────────────────────────────────────────┘
```

### Phase 2: Module Refactoring (Months 3-4)
- Convert Python agents to Nextflow modules
- Implement proper channel communication
- Add container definitions for each module
- Create unit tests for modules

### Phase 3: Full Integration (Months 5-6)
- Complete workflow orchestration
- Advanced features implementation
- Performance optimization
- Production deployment

## Example Implementation

### Basic Nextflow Workflow for F.A.D.E

```nextflow
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define input parameters
params.query = "Find molecules targeting KRAS G12D with good BBB permeability"
params.output_dir = "results"
params.skip_structure = false
params.num_molecules = 100
params.max_iterations = 3

// Import modules
include { targetSelection } from './modules/local/target_selector'
include { structurePrediction } from './modules/local/structure_predictor'
include { generateMolecules } from './modules/local/molecule_generator'
include { performDocking } from './modules/local/docking'
include { scoreAndRank } from './modules/local/scoring'
include { refineMolecules } from './modules/local/refiner'
include { generateReport } from './modules/local/reporting'

// Define the main workflow
workflow {
    // Create input channel from query
    query_ch = Channel.value(params.query)
    
    // Step 1: Target Selection
    target_results = targetSelection(query_ch)
    
    // Step 2: Structure Prediction or Retrieval
    if (!params.skip_structure) {
        structure = structurePrediction(target_results.fasta)
    } else {
        structure = fetchPDBStructure(target_results.pdb_id)
    }
    
    // Step 3: Binding Site Analysis
    binding_sites = analyzeBindingSites(structure)
    
    // Step 4: Iterative Molecule Generation and Optimization
    iteration_ch = Channel.from(1..params.max_iterations)
    
    molecules_ch = Channel.empty()
    
    iteration_ch.each { iteration ->
        // Generate molecules
        new_molecules = generateMolecules(
            target_results.target_info,
            binding_sites,
            iteration
        )
        
        // Screen molecules
        screened = rapidScreen(new_molecules)
        
        // Perform docking
        docking_results = performDocking(
            structure,
            screened
        ).flatten()
        
        // Score and rank
        ranked = scoreAndRank(docking_results.collect())
        
        // Check if refinement needed
        if (iteration < params.max_iterations) {
            refined = refineMolecules(ranked.failed_molecules)
            molecules_ch = molecules_ch.mix(refined)
        }
        
        molecules_ch = molecules_ch.mix(ranked.successful_molecules)
    }
    
    // Step 5: Lead Optimization
    leads = optimizeLeads(molecules_ch.collect())
    
    // Step 6: Generate Final Report
    generateReport(
        target_results.target_info,
        structure,
        leads
    )
}

// Workflow entry point
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Duration: $workflow.duration"
}
```

### Module Example: Target Selector

```nextflow
process targetSelection {
    tag "$query"
    publishDir "${params.output_dir}/target", mode: 'copy'
    
    container 'fade-base:latest'
    
    input:
    val query
    
    output:
    path 'target_info.json', emit: target_info
    path 'protein.fasta', emit: fasta
    val pdb_id, emit: pdb_id
    
    script:
    """
    python ${projectDir}/agents/target_selector.py \\
        --query "$query" \\
        --output target_info.json \\
        --fasta protein.fasta
    
    # Extract PDB ID if available
    pdb_id=\$(jq -r '.pdb_id // "none"' target_info.json)
    echo \$pdb_id > pdb_id.txt
    """
}
```

### Configuration Example

```groovy
// nextflow.config
profiles {
    northeastern {
        process {
            executor = 'slurm'
            queue = 'shared'
            
            withLabel: gpu {
                queue = 'gpu'
                clusterOptions = '--gres=gpu:1'
            }
            
            withLabel: highmem {
                queue = 'highmem'
                memory = '64 GB'
            }
        }
        
        singularity {
            enabled = true
            autoMounts = true
            cacheDir = "$SCRATCH/singularity_cache"
        }
    }
    
    local {
        process.executor = 'local'
        docker.enabled = true
    }
    
    aws {
        process.executor = 'awsbatch'
        aws.region = 'us-east-1'
        aws.batch.cliPath = '/usr/local/bin/aws'
    }
}

params {
    // Default parameters
    max_cpus = 16
    max_memory = '128.GB'
    max_time = '240.h'
    
    // F.A.D.E specific parameters
    gemini_api_key = "$GEMINI_API_KEY"
    alphafold_container = '/shared/containers/alphafold3.sif'
    schrodinger_path = '/opt/schrodinger2024-4'
}
```

## Running the Workflow

### Basic Execution
```bash
# Run with default parameters
nextflow run main.nf

# Run with custom query
nextflow run main.nf --query "Find EGFR inhibitors with low toxicity"

# Resume from previous run
nextflow run main.nf -resume

# Run with specific profile
nextflow run main.nf -profile northeastern

# Run with resource limits
nextflow run main.nf --max_cpus 8 --max_memory 32.GB
```

### Advanced Options
```bash
# Generate execution report
nextflow run main.nf -with-report report.html

# Create timeline
nextflow run main.nf -with-timeline timeline.html

# Enable trace
nextflow run main.nf -with-trace

# Use Nextflow Tower for monitoring
nextflow run main.nf -with-tower
```

## Benefits Summary

1. **Reduced Complexity**
   - Replace complex shell scripts with declarative workflow
   - Automatic dependency management
   - Clear separation of concerns

2. **Improved Reliability**
   - Automatic error handling and retries
   - Checkpoint and resume capability
   - Robust resource management

3. **Enhanced Scalability**
   - Easy scaling from 10 to 10,000 molecules
   - Efficient parallel processing
   - Dynamic resource allocation

4. **Better Monitoring**
   - Real-time progress tracking
   - Detailed execution reports
   - Resource usage analytics

5. **Increased Reproducibility**
   - Container integration
   - Version control for workflows
   - Consistent execution across platforms

## Next Steps

1. **Immediate Actions**
   - Install Nextflow on the HPC cluster
   - Create basic workflow skeleton
   - Test with simple target selection process

2. **Short-term Goals (1-2 months)**
   - Implement hybrid approach for all agents
   - Set up monitoring with Nextflow Tower
   - Create development and production profiles

3. **Medium-term Goals (3-4 months)**
   - Refactor agents as Nextflow modules
   - Implement parallel processing for molecule generation
   - Add comprehensive error handling

4. **Long-term Goals (5-6 months)**
   - Full production deployment
   - Integration with cloud resources
   - Community contribution to nf-core

## Conclusion

Integrating Nextflow into the F.A.D.E project will transform it from a collection of scripts into a robust, scalable, and maintainable drug discovery pipeline. The benefits include improved reliability, better resource utilization, and the ability to easily scale computational drug discovery efforts.

The gradual integration approach ensures that existing work is preserved while gaining immediate benefits from Nextflow's workflow management capabilities. This architecture provides a solid foundation for future enhancements and scaling of the F.A.D.E platform.

---

*Document Version: 1.0*  
*Last Updated: August 2025*  
*Author: F.A.D.E Development Team*
