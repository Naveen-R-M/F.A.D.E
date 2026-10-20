# F.A.D.E Drug Discovery Report: EGFR kinase domain (7SI1 chain A construct)

*Generated: 2026-10-20 13:09:58*

**Query:** Generate candidate ligands for EGFR kinase domain, branch 2

## Pipeline Summary

- **Target:** EGFR kinase domain (7SI1 chain A construct)
- **Structure source:** PDB
- **Selected pocket:** pocket_11 (druggability 0.51)
- **Molecules screened:** 2
- **Likely binders:** 0
- **Screening method:** Boltz-2 binding-affinity prediction

## Key Insights

- Screened 2 candidate ligand(s) against EGFR kinase domain (7SI1 chain A construct) using Boltz-2 binding-affinity prediction.
- 0 candidate(s) exceeded the binder-probability threshold (0.5).
- Structure source for this run: PDB (branch determines how the pocket was found, screening is identical).
- Top candidate diffsbdd_36906525_1 has predicted affinity value 0.947 (log10 IC50, lower is stronger) and P(binder)=0.128.

## Top Candidates

Ranked by Boltz-2 predicted affinity value (lower log10 IC50 = stronger); binder probability breaks ties.

| Rank | Molecule | P(binder) | Affinity value (log10 IC50) | Pred. IC50 (µM) | QED |
|------|----------|-----------|------------------------------|-----------------|-----|
| 1 | diffsbdd_36906525_1 | 0.128 | 0.947 | 8.86 | 0.54 |
| 2 | diffsbdd_36906525_3 | 0.068 | 1.700 | 50.1 | 0.73 |

### Top candidate structures (SMILES)

- **diffsbdd_36906525_1**: `BrNBr`
- **diffsbdd_36906525_3**: `O=C(O)[C@H]1CCCC(=C2CC=C(C(F)(F)F)CC2)C1`

---

*Affinity values are Boltz-2 predictions (log10 IC50); the derived kcal/mol and IC50 columns are informational conversions and are assay-dependent. Use predicted affinity and binder probability for ranking, and validate top candidates experimentally.*