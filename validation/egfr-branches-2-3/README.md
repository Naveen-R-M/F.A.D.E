# EGFR Branch 2 and Branch 3 validation

Production FADE LangGraph runs on Explorer, using private repaired environments and fpocket-only detection. Code baseline: `6da27f3`. Neither run validates the shared installations.

Both used EGFR's **331-residue 7SI1 chain A construct**, not full-length EGFR. RCSB describes 7SI1 as an apo EGFR kinase-domain structure (1.6 Å resolution). Requested scale: 5 generated molecules, at most 2 affinity screens.

## Branch 2: completed

- Application job `11262418`, runtime 19m18s, run ID `bbe6aac4`.
- Explicit existing structure: PDB `7SI1`; no Boltz-2 structure prediction.
- Selected pocket 11: druggability approximately 0.507, volume approximately 342 Å³, 7 residues.
- 5 generated; 4 passed filtering; 2 selected and screened; 0 exceeded P(binder)=0.5.

| Molecule | SMILES | P(binder) | Raw Boltz affinity value |
| --- | --- | ---: | ---: |
| `diffsbdd_36906525_1` | `BrNBr` | 0.1285 | 0.9474 |
| `diffsbdd_36906525_3` | `O=C(O)[C@H]1CCCC(=C2CC=C(C(F)(F)F)CC2)C1` | 0.0682 | 1.7002 |

**Quality caveat:** `BrNBr` is not a credible drug lead despite passing these filters. Completion demonstrates software execution, not medicinal chemistry quality. Pocket correspondence to EGFR's ATP-binding site has not been validated.

Start with `branch2/final_report.md` and `branch2/screening_results.json`. The report's derived affinity quantities are application interpretations; predictions are not experimental binding measurements. Original state, log, input structure, selected pocket, SDF and predicted complexes are included.

## Branch 3: structure succeeded; stopped at pocket cutoff

- Initial application job `11262420` failed during initialization with SQLite `locking protocol` while Branch 2 ran concurrently. Its log is included separately.
- Rerun alone: application job `11264041`, runtime 11m37s, run ID `3cabff54`.
- Forced Boltz-2 prediction from the matching sequence; prediction succeeded.
- fpocket detected 19 pockets; **none passed the unchanged 0.5 druggability cutoff**.
- Best: pocket 15, score 0.420. Its coordinates are included for inspection only; it was **not selected for generation**.
- No DiffSBDD generation, affinity screens or final report for Branch 3.

`branch3/` includes state, logs, predicted PDB/mmCIF, Boltz confidence scores, fpocket metrics, and the below-cutoff pocket. The graph's failure-state snapshot does not retain all completed-node fields; the structure files and log document the successful prediction.

## Reproduction and provenance

`reproduction/` contains copies of the original driver and SLURM script. They call the production graph with explicit target inputs and private env overrides. They contain absolute paths from the original cluster run and **require path adjustment** before reuse. The original RCSB metadata JSONs are not included; retrieve entry 7SI1 and polymer entity 1 from the RCSB API if reproducing the driver.

Artifacts retain their original embedded paths and formatting. The cutoff was not relaxed; no additional run was performed for this upload. Further work should check ATP-site pocket correspondence and strengthen filtering before expanding generation. Results-only branch; no merge into main intended.
