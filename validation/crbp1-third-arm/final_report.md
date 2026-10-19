# F.A.D.E Drug Discovery Report: sp|P09455|RET1_HUMAN Cellular Retinol-binding protein 1 (CRBP1) GN=RBP1

*Generated: 2026-10-19 14:33:43*

**Query:** Find inhibitors for CRBP1 from sequence

## Pipeline Summary

- **Target:** sp|P09455|RET1_HUMAN Cellular Retinol-binding protein 1 (CRBP1) GN=RBP1
- **Structure source:** Boltz2
- **Selected pocket:** pocket_1 (druggability 0.93)
- **Molecules screened:** 2
- **Likely binders:** 0
- **Screening method:** Boltz-2 binding-affinity prediction

## Key Insights

- Screened 2 candidate ligand(s) against sp|P09455|RET1_HUMAN Cellular Retinol-binding protein 1 (CRBP1) GN=RBP1 using Boltz-2 binding-affinity prediction.
- 0 candidate(s) exceeded the binder-probability threshold (0.5).
- Structure source for this run: Boltz2 (branch determines how the pocket was found, screening is identical).
- Top candidate diffsbdd_81107de4_3 has predicted affinity value 1.854 (log10 IC50, lower is stronger) and P(binder)=0.157.

## Top Candidates

Ranked by Boltz-2 predicted affinity value (lower log10 IC50 = stronger); binder probability breaks ties.

| Rank | Molecule | P(binder) | Affinity value (log10 IC50) | Pred. IC50 (µM) | QED |
|------|----------|-----------|------------------------------|-----------------|-----|
| 1 | diffsbdd_81107de4_3 | 0.157 | 1.854 | 71.4 | 0.67 |
| 2 | diffsbdd_81107de4_5 | 0.244 | 2.401 | 252 | 0.83 |

### Top candidate structures (SMILES)

- **diffsbdd_81107de4_3**: `Oc1ccccc1OC(O)(F)F`
- **diffsbdd_81107de4_5**: `CCN(C(=O)O)C(=O)c1ccc(C)cc1C`

---

*Affinity values are Boltz-2 predictions (log10 IC50); the derived kcal/mol and IC50 columns are informational conversions and are assay-dependent. Use predicted affinity and binder probability for ranking, and validate top candidates experimentally.*