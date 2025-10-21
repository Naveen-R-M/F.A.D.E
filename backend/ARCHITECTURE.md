# F.A.D.E Backend - LangGraph Architecture

## Clean Architecture ✅

This is a pure **LangGraph** implementation - no custom agent classes!

## Directory Structure

```
fade/
├── agents/           # LangGraph nodes (functions, not classes)
│   ├── research/     # Target research nodes
│   ├── structure/    # Structure prediction nodes (Boltz-2)
│   ├── targeting/    # Pocket detection nodes (SILCS)
│   ├── invention/    # Molecule generation nodes (DiffSBDD)
│   ├── screening/    # Screening nodes (Boltz-2)
│   └── analysis/     # Analysis and reporting nodes
│
├── state/            # LangGraph state definition
│   └── langgraph_state.py  # DrugDiscoveryState TypedDict
│
├── tools/            # External API clients
│   ├── uniprot_api.py   # UniProt protein database
│   ├── chembl_api.py    # ChEMBL compound database
│   └── rcsb_api.py      # RCSB PDB structure database
│
├── workflows/        # LangGraph workflow definitions
│   └── drug_discovery.py  # Main StateGraph workflow
│
├── utils/            # Utilities
│   └── logging.py    # Logging configuration
│
└── config.py         # Configuration management
```

## Key Design Principles

### 1. **Pure LangGraph Implementation**
- All agents are **node functions**, not classes
- State flows through nodes via `DrugDiscoveryState`
- Conditional routing with `add_conditional_edges`
- Compiled graph with `workflow.compile()`

### 2. **Node Pattern**
```python
def node_name(state: DrugDiscoveryState) -> Dict[str, Any]:
    """LangGraph node that processes state and returns updates."""
    # Process state
    # Return state updates (not full state)
    return {"field_to_update": new_value}
```

### 3. **Workflow Definition**
```python
workflow = StateGraph(DrugDiscoveryState)
workflow.add_node("node_name", node_function)
workflow.add_conditional_edges(...)
app = workflow.compile()
```

## Current Implementation Status

### ✅ Completed
- **Research Module**: `target_research_node`, `validate_target_node`
- **API Tools**: UniProt, ChEMBL, RCSB clients
- **Core Workflow**: LangGraph StateGraph setup
- **State Management**: DrugDiscoveryState definition

### 🚧 To Implement (as LangGraph nodes)
- **Structure Module**: Boltz-2 integration for structure prediction
- **Targeting Module**: SILCS integration for pocket detection
- **Invention Module**: DiffSBDD for molecule generation
- **Screening Module**: Boltz-2 for binding affinity prediction
- **Analysis Module**: Final report generation

## Running the Pipeline

```bash
cd backend
python main.py "Find inhibitors for KRAS G12C"
```

## Key Technologies
- **Orchestration**: LangGraph (not custom agents!)
- **LLM**: LangChain with configurable models (GPT-4, Claude, Gemini)
- **Structure Prediction**: Boltz-2 (not AlphaFold3)
- **Pocket Detection**: SILCS
- **Molecule Generation**: DiffSBDD
- **Screening**: Boltz-2 (not traditional docking)

## Important Notes
- **No custom BaseAgent class** - everything is LangGraph nodes
- **No agent classes** - only node functions
- **State updates** - nodes return partial updates, not full state
- **Message history** - using LangGraph's `add_messages` pattern
