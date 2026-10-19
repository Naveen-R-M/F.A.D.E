# CRBP1 third-arm validation outputs

Results from the production FADE CLI/LangGraph run, SLURM job `11230565` (21 min 7 sec), using code at `6da27f3`.

## Route and outcome

CRBP1 FASTA (135 amino acids) → Boltz-2 structure → fpocket → DiffSBDD → RDKit filtering → Boltz-2 affinity → final report.

- Selected pocket: `pocket_1`, fpocket druggability approximately 0.93, volume approximately 934 Å³.
- Generated 5 molecules; **all 5 passed filtering**, with 2 selected under the screening cap (see `application-state.json`, `filter_statistics`).
- Molecule 3: P(binder) 0.1569, affinity value 1.8536.
- Molecule 5: P(binder) 0.2438, affinity value 2.4013.
- Neither exceeded the 0.5 binder-probability threshold. These are model predictions, not experimentally validated inhibitors.
- Pipeline completed with no recorded error and a final report.

## Files

| File | Contents |
| --- | --- |
| `final_report.md` | Human-readable report; start here |
| `screening_results.json` | Screening scores, SMILES, ranking and descriptors |
| `crbp1.fasta` | Input sequence |
| `predicted_structure.pdb` | Boltz-2 prediction converted to PDB |
| `pocket_1.pdb` | Selected pocket |
| `generated_molecules.sdf` | Five generated molecules |
| `complex_molecule_3.cif` | Predicted complex for `diffsbdd_81107de4_3` |
| `complex_molecule_5.cif` | Predicted complex for `diffsbdd_81107de4_5` |
| `application-state.json` | Complete recorded workflow state |
| `run.log` | Application job log, including regression tests |

Artifacts are copied without rewriting their embedded original cluster paths. Those paths describe the original run and are not portable download links.

## Environment qualification

This run used **private repaired environments** on Explorer, with `FADE_LOCAL_HPC=1`, `DIFFSBDD_PYTHON` and `BOLTZ2_ENV_SETUP` overrides. It does **not** validate the shared environments. The later shared-environment repair attempt was blocked by project storage quota and is separate from these results. The app runtime also used a private Python 3.10 environment. App-job regression tests: 8 passed.

CLI (from `backend/`, with the run's environment variables configured):

```bash
python main.py 'Find inhibitors for CRBP1 from sequence' \
  --fasta crbp1.fasta --force-prediction --pocket-mode fpocket \
  --no-rich --state-output application-state.json
```

Generation and screening caps were 5 and 2 respectively. Output branch only; no merge into `main` is intended.
