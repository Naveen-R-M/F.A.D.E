# F.A.D.E Minimal Implementation Todo List

## Current Status: RCSB Target Selection Only

### ✅ COMPLETED
- **RCSB Target Selection**: Complete RCSB PDB-based target identification
  - LLM-powered query parsing
  - Intelligent structure search and ranking
  - PDB structure download
  - FASTA sequence extraction
  - Requirements parsing

### 🔄 CURRENT PIPELINE
```
Query → RCSB Target Selection → END
```

**Outputs:**
- `target_info.json`: Target and PDB metadata
- `structure.pdb`: Downloaded PDB structure
- `protein.fasta`: Protein sequence
- `requirements.json`: Parsed drug requirements

## 📋 TODO: BUILD AS WE GO

### Next Steps (To be added based on instructions)
- [ ] **Step 2**: TBD - Next component to be specified
- [ ] **Step 3**: TBD - Next component to be specified  
- [ ] **Step 4**: TBD - Next component to be specified
- [ ] **Step 5**: TBD - Next component to be specified
- [ ] **Step 6**: TBD - Next component to be specified

### 🗂️ BACKED UP COMPONENTS
All previous pipeline components have been moved to `nextflow/modules/backup/` and `nextflow/bin/backup/`:
- Structure Prediction (AlphaFold3-based)
- Binding Site Analysis
- Molecule Generation  
- Docking
- Lead Optimization
- Reporting

These can be restored and adapted as needed.

## 🚀 CURRENT TESTING

To test the minimal pipeline:

```bash
cd nextflow
export GEMINI_API_KEY="your_api_key"
nextflow run main.nf --query "Find drugs targeting EGFR for lung cancer"
```

## 📊 PIPELINE STATUS

- **Active Modules**: 1 (RCSB Target Selection)
- **Backed Up Modules**: 7 (Available for restoration)  
- **Pipeline Complexity**: Minimal (Single step)
- **Resource Requirements**: Low (CPU only, no GPU)
- **Dependencies**: RCSB PDB API, Gemini API

---

**Ready for next component based on your instructions.**
