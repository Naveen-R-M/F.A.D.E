# ATP-site correspondence check

Analysis of the **existing** EGFR outputs, using holo 4G5J chain A ligand 0WN (residue 1101) as the reference. No fpocket rerun, model rerun, threshold adjustment or molecule generation was performed.

## Main findings

The reference ATP site comprises 24 residues having a heavy atom within 5 Å of 0WN. All 24 have sequence-aligned counterparts in both EGFR structures.

| Structure and pocket | Druggability | Reference ATP-site residues represented | Fraction of pocket residues belonging to reference site |
|---|---:|---:|---:|
| 7SI1 pocket 11 (used for generation) | 0.507 | 1/24 | 1/7 (14.3%) |
| Predicted pocket 15 (proposed relaxed-cutoff candidate) | 0.420 | 0/24 | 0% |
| 7SI1 pocket 18 (largest residue overlap) | 0.067 | 9/24 | 100% |
| 7SI1 pocket 16 | 0.051 | 6/24 | 100% |
| Predicted pocket 13 (largest residue overlap) | 0.024 | 10/24 | 90.9% |

**These results do not support identifying predicted pocket 15 as the ATP site or relaxing the threshold to select it for an ATP-site experiment.** Branch 2 completed computationally, but its selected pocket 11 likewise has very limited correspondence to the reference ATP site. This does not rule out other biological roles for these cavities; none are established here.

The low-scoring overlapping pockets represent partial coverage of reference-site residues, not independently validated usable ATP cavities. Their druggability scores remain below the original 0.5 threshold. No threshold was changed.

## Alignment caveat

Global sequence alignments yielded 290 paired Cα atoms for 7SI1 and 307 for predicted EGFR, with 100% identity among paired residues. All-atom-pair Cα rigid-fit RMSDs were **6.197 Å** and **9.338 Å**, respectively. Reference-site Cα RMSDs after these global fits were **2.012 Å** and **4.558 Å**. Separate site-only fits yielded **0.527 Å** and **0.517 Å**.

Consequently, whole-structure fits do not optimally align the local site. Spatial distances and ligand coverage in `metrics.md` and `validation.json` are explicitly **global-fit descriptive measurements**, not definitive site-local pose assessments. Sequence-mapped residue overlap, the primary evidence for the conclusions above, does not depend on the rigid fit. Local-site similarity also means these results do not establish that the predicted ATP-site backbone itself is poor.

This analysis does not explain the manuscript's different 7SI1 pocket scores or prove that either preparation reproduces its pocket detection. Pocket IDs are not conserved anatomical labels. The original manuscript input preparation, fpocket version/settings and output files are needed to reconcile that discrepancy.

## Files and provenance

- `metrics.md`: measurements for every existing pocket.
- `validation.json`: reference residues, residue intersections, alignments, transformations and measurements.
- `4G5J.pdb`: reference downloaded from RCSB.
- `*_aligned.pdb`: structures and selected pockets in the reference coordinate frame for inspection.
- `validate.py`, `validate.sbatch`, `validation-11266530.log`: calculation and execution record (Biopython 1.88, NumPy 2.2.6). Scripts retain original absolute cluster paths; adjust for reuse.
- `input_sha256.json`: checksums of all files consumed, keyed by original cluster path.
- `existing_fpocket/`: existing input structures, info files, pocket atom files and sphere coordinates used by the calculation, copied for reproducibility. This is not new fpocket output.

The script checks transform direction/RMSD consistency and rigid-distance preservation. No affinity pose validation or ligand generation was performed in this analysis.
