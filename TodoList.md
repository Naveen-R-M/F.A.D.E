# F.A.D.E Implementation Todo List

## Initial Project Setup
- ✅ Create conda environment (fade)
- ✅ Set up project structure
- ✅ Create README.md
- ✅ Create .gitignore
- ✅ Create .env.example
- ✅ Create .env file with API keys

## Step 1: Target Selector Agent (Current Focus)
- ✅ Create base agent class (BaseAgent)
- ✅ Implement Gemini client utility
- ✅ Implement UniProt client utility
- ✅ Implement Config Generator utility
- ✅ Create job templates for structure prediction and docking
- ✅ Implement TargetSelector agent
- ✅ Create main script to demonstrate usage
- ✅ Write unit tests for TargetSelector agent
- ✅ Install required packages for TargetSelector agent
- ✅ Test TargetSelector agent with example queries

## Step 1B: Agentic Target Selector Enhancements
- ✅ Design agentic architecture (see docs/design/agentic_target_selector.md)
- ✅ Implement agentic mixin for base agent
- ✅ Create error analysis component for API errors
- ✅ Implement query reformulation for failed API calls
- ✅ Add sequence validation for scientific accuracy
- ✅ Create adaptive search strategy for protein retrieval
- ✅ Integrate LLM-powered decision making
- ✅ Update unit tests for agentic components
- ✅ Benchmark agentic vs. non-agentic performance

## Step 2: Structure Predictor Agent (Next)
- [ ] Implement StructurePredictor agent
- [ ] Implement PDB file parser/processor
- [ ] Create structure validation utilities
- [ ] Implement binding site detection
- [ ] Write unit tests for StructurePredictor agent
- [ ] Test StructurePredictor agent with example proteins

## Step 3: Molecule Generator Agent
- [ ] Implement MoleculeGenerator agent
- [ ] Integrate with RDKit
- [ ] Implement molecule property calculators
- [ ] Create molecule filtering utilities
- [ ] Write unit tests for MoleculeGenerator agent
- [ ] Test MoleculeGenerator agent with example targets

## Step 4: Evaluator Agent
- [ ] Implement Evaluator agent
- [ ] Create property prediction models/interfaces
- [ ] Implement scoring functions
- [ ] Write unit tests for Evaluator agent
- [ ] Test Evaluator agent with example molecules

## Step 5: Docking Agent
- [ ] Implement Docking agent
- [ ] Create interfaces for docking software
- [ ] Implement pose analysis utilities
- [ ] Write unit tests for Docking agent
- [ ] Test Docking agent with example receptor-ligand pairs

## Step 6: Refiner Agent
- [ ] Implement Refiner agent
- [ ] Create molecule modification utilities
- [ ] Implement learning from docking results
- [ ] Write unit tests for Refiner agent
- [ ] Test Refiner agent with example molecules

## Step 7: Integration and Workflow
- [ ] Create workflow orchestration
- [ ] Implement memory manager for tracking history
- [ ] Create agent communication interfaces
- [ ] Implement parallel execution capabilities
- [ ] Write unit tests for workflow orchestration
- [ ] Test full workflow with example queries

## Step 8: Documentation and Examples
- [ ] Write detailed documentation for each agent
- [ ] Create example queries and expected outputs
- [ ] Write installation and setup guide
- [ ] Create usage tutorials
- [ ] Document API and interfaces

## Step 9: Deployment and CI/CD
- [ ] Set up continuous integration
- [ ] Create containerized version
- [ ] Write deployment scripts for HPC
- [ ] Create user-friendly CLI
- [ ] Implement logging and monitoring
