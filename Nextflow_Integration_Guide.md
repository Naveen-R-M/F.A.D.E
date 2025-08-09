# F.A.D.E Nextflow Integration Guide

## Overview

The F.A.D.E project now includes a fully integrated Nextflow workflow that utilizes the actual Python agents instead of placeholder processes. This integration provides robust workflow management, better resource allocation, and improved fault tolerance.

## Architecture Changes

### Before Integration
- Standalone Python agents (`main.py`)
- Manual job submission and monitoring
- Limited parallelization
- Basic error handling

### After Integration
- **Nextflow workflow orchestration** with automatic dependency management
- **HPC-optimized resource allocation** with SLURM integration
- **Parallel execution** of independent tasks
- **Resume capability** from failed steps
- **Comprehensive monitoring** and reporting

## File Structure

```
F.A.D.E/
├── nextflow/
│   ├── main.nf                           # Main workflow definition
│   ├── nextflow.config                   # Configuration and resource allocation
│   ├── modules/                          # Individual workflow modules
│   │   ├── target_selection.nf
│   │   ├── structure_prediction.nf
│   │   ├── binding_site.nf
│   │   ├── molecule_generation.nf
│   │   ├── docking.nf
│   │   ├── lead_optimization.nf
│   │   └── reporting.nf
│   └── bin/                              # Python wrapper scripts
│       ├── run_target_selector.py
│       ├── run_structure_predictor.py
│       ├── run_binding_site_analysis.py
│       ├── run_molecule_generator.py
│       ├── run_docking.py
│       ├── run_lead_optimization.py
│       └── run_reporting.py
├── run_nextflow_integrated.sh            # Main execution script
├── test_integrated_workflow.sh           # Integration test suite
└── [existing Python agents]              # Original agent implementations
```

## Key Integration Features

### 1. Actual Agent Integration
- **Real Python agents** called via Nextflow wrapper scripts
- **Proper error handling** and fallback mechanisms
- **Data validation** between pipeline stages
- **Full integration** with existing utilities (Gemini client, UniProt client, etc.)

### 2. HPC Resource Management
```groovy
process structurePrediction {
    queue = 'gpu'
    clusterOptions = '--gres=gpu:1'
    cpus = 8
    memory = '32.GB'
    time = '4.h'
    module = ['AlphaFold3/3.0.0']
}
```

### 3. Environment Management
- **Conda environment integration** for Python dependencies
- **Module loading** for HPC software (AlphaFold3, Schrodinger, AutoDock Vina)
- **Environment variable handling** for API keys and paths

### 4. Data Flow Management
```
Query → Target Info → Structure → Binding Sites → Molecules → Docking → Optimization → Report
```

## Usage

### Quick Start
```bash
# Basic usage
./run_nextflow_integrated.sh "Find molecules targeting KRAS G12D with good BBB permeability"

# With custom parameters
./run_nextflow_integrated.sh -m 200 -i 3 -d glide "Design EGFR inhibitors for lung cancer"
```

### Testing
```bash
# Run integration tests
./test_integrated_workflow.sh

# Dry run with stub processes
./run_nextflow_integrated.sh --dry-run "Test query"

# Debug mode with detailed tracing
./run_nextflow_integrated.sh --debug --trace "Debug query"
```

### Advanced Usage
```bash
# Resume failed pipeline
./run_nextflow_integrated.sh --resume "Original query"

# Custom output directory
./run_nextflow_integrated.sh -o /path/to/results "Query"

# Local testing (no SLURM)
./run_nextflow_integrated.sh -p local "Query for local testing"
```

## Configuration

### Environment Setup
1. **Conda environment**: Ensure `${SCRATCH}/conda-envs/fade` exists with all dependencies
2. **API keys**: Set `GEMINI_API_KEY` in environment or `.env` file
3. **Modules**: Verify access to required HPC modules (AlphaFold3, Schrodinger, etc.)

### Nextflow Profiles
- **`northeastern`** (default): SLURM execution with HPC modules
- **`local`**: Local execution for testing
- **`debug`**: Enhanced logging and monitoring
- **`stub`**: Stub processes for workflow validation

### Resource Allocation
```groovy
// CPU-intensive processes
process_high: 8 CPUs, 32GB RAM, 8h
// Standard processes  
process_medium: 4 CPUs, 8GB RAM, 2h
// Lightweight processes
process_low: 1 CPU, 2GB RAM, 30m
```

## Agent Integration Details

### 1. Target Selector Agent
- **Input**: Natural language query
- **Process**: Uses `TargetSelector` Python class with Gemini API
- **Output**: `target_info.json`, `protein.fasta`, `requirements.json`

### 2. Structure Predictor Agent
- **Input**: Target info and FASTA sequence
- **Process**: Uses `StructurePredictor` with AlphaFold3 container
- **Output**: `structure.pdb`, `prepared_receptor.pdb`
- **Resources**: GPU partition, 8 CPUs, 32GB RAM

### 3. Molecule Generator Agent
- **Input**: Requirements and binding sites
- **Process**: Uses `MoleculeGenerator` with LLM guidance and RDKit
- **Output**: `molecules.sdf`, `molecules.json`
- **Resources**: CPU partition, 4 CPUs, 16GB RAM

### 4. Docking Process
- **Input**: Molecules and prepared receptor
- **Process**: AutoDock Vina or Schrodinger Glide
- **Output**: `docking_results.json`, `top_hits.json`
- **Resources**: CPU partition, 8 CPUs, 16GB RAM

## Monitoring and Debugging

### Pipeline Monitoring
```bash
# Check running jobs
squeue -u $USER

# Monitor Nextflow execution
tail -f .nextflow.log

# Check specific process logs
ls work/*/
```

### Generated Reports
- **Timeline**: `${output_dir}/timeline.html`
- **Execution Report**: `${output_dir}/execution_report.html`
- **DAG Visualization**: `${output_dir}/workflow_dag.svg`
- **Trace Log**: `${output_dir}/pipeline_trace.txt`

### Troubleshooting
1. **Environment Issues**: Check conda environment path and activation
2. **API Failures**: Verify GEMINI_API_KEY is set correctly
3. **Module Loading**: Ensure HPC modules are available
4. **Resource Limits**: Check SLURM queue limits and job timeouts
5. **Python Imports**: Verify PYTHONPATH includes F.A.D.E root directory

## Performance Optimization

### Parallel Execution
- Multiple molecules processed simultaneously during docking
- Independent binding site analysis can run in parallel
- Optimization cycles can be parallelized

### Resource Efficiency
- **GPU resources** allocated only for structure prediction
- **High-memory nodes** used only when needed
- **Automatic cleanup** of intermediate files (optional)

### Caching
- Nextflow automatically caches completed processes
- Resume from failures without re-running successful steps
- Intermediate results preserved for inspection

## Migration from Python-only Version

### Advantages of Nextflow Integration
1. **Better resource management** on HPC clusters
2. **Automatic retry** and error recovery
3. **Parallel execution** of independent tasks
4. **Resume capability** for long-running workflows
5. **Built-in monitoring** and reporting
6. **Scalable execution** across multiple compute nodes

### Maintaining Compatibility
- Original Python agents remain unchanged
- `main.py` still works for standalone execution
- Existing monitoring scripts continue to function
- All utilities and libraries preserved

## Future Enhancements

### Planned Improvements
1. **Dynamic resource allocation** based on molecule count
2. **Multi-target workflows** for comparative studies
3. **Advanced optimization** with multiple iterations
4. **Cloud execution** support (AWS, Google Cloud)
5. **Real-time monitoring** dashboard

### Extension Points
- Additional docking software integration
- Machine learning model integration
- Custom property calculators
- Advanced visualization tools

## Testing and Validation

### Test Suite
```bash
# Run all integration tests
./test_integrated_workflow.sh

# Test specific components
./run_nextflow_integrated.sh --dry-run "Test query"
```

### Validation Checklist
- [ ] Conda environment properly configured
- [ ] API keys available
- [ ] HPC modules accessible
- [ ] Wrapper scripts executable
- [ ] Python imports working
- [ ] Nextflow configuration valid
- [ ] Resource allocations appropriate

## Best Practices

### Query Formulation
- Be specific about target proteins
- Include desired molecular properties
- Specify therapeutic context
- Mention selectivity requirements

### Resource Management
- Use appropriate profiles for different environments
- Monitor disk usage in scratch space
- Clean up work directories periodically
- Use resume functionality for failed runs

### Debugging
- Start with dry runs for new queries
- Use debug mode for troubleshooting
- Check individual process logs
- Verify environment setup before full runs

This integration represents a significant advancement in the F.A.D.E project, providing a robust, scalable foundation for AI-powered drug discovery workflows.
