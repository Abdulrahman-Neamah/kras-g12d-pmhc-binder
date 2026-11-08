# Challenge 5 — Boltz-2 co-folding screen, ipSAE analysis

60 predictions: 20 designs × {G12D, WT, apo}, plus 2 designs × 2 extra seeds.
Boltz-2 (2.0.3), 1 diffusion sample, 3 recycling steps, trifast kernel off.
Chains A = HLA-A*11:01, B = beta-2m, C = peptide, Z = binder.
ipSAE from DunbrackLab/IPSAE, pae_cutoff 10, dist_cutoff 10.
Values below are the C->Z direction (peptide aligned, binder scored) for the peptide
interface and Z->A for the HLA interface; both directions are in ipsae_scores.csv.

## Headline

Of 20 designs, only 5 have any confidently predicted interface to the G12D peptide
(ipSAE > 0.15). Of those 5, only 2 discriminate against the wild-type peptide, and both
are on the same RFD4-Proteina scaffold (d010).

| design | n seeds | G12D pep | WT pep | selectivity | G12D HLA | apo HLA | apo selectivity |
|---|---|---|---|---|---|---|---|
| mpnn_d010_14 | 1 | 0.864 | 0.491 | +0.373 | 0.924 | 0.760 | +0.164 |
| **mpnn_d010_13** | **3** | **0.833** | **0.046** | **+0.787** | **0.913** | **0.000** | **+0.913** |
| rfd4_d015 | 1 | 0.511 | 0.538 | -0.027 | 0.622 | n/a | n/a |
| rfd4_d031 | 1 | 0.237 | 0.203 | +0.034 | 0.390 | n/a | n/a |
| mpnn_d031_03 | 1 | 0.224 | n/a | n/a | 0.405 | n/a | n/a |

The remaining 15 designs score at or near 0.000 on every target, including G12D.

## mpnn_d010_13, three seeds

| seed | G12D pep | WT pep | apo HLA | backbone gap (G12D) |
|---|---|---|---|---|
| default | 0.820 | 0.000 | 0.000 | 3.73 A |
| 101 | 0.858 | 0.000 | 0.000 | 3.78 A |
| 202 | 0.820 | 0.137 | 0.000 | 3.74 A |

G12D binding is seed-stable (SD 0.018). Both counterfactuals collapse to zero in 2 of 3
seeds and stay far below G12D in the third. No interpenetration in any run.

## ipSAE disagrees with ipTM, and the disagreement matters

Boltz's ipTM is a learned prediction head; ipSAE is computed from the PAE matrix with a
strict per-pair cutoff. They conflict here:

| design / target | ipTM (binder-peptide) | ipSAE (binder-peptide) |
|---|---|---|
| mpnn_d010_13 vs WT | 0.52 - 0.73 | 0.000 - 0.137 |
| mpnn_d002_09 vs G12D | 0.685 | 0.123 |
| mpnn_d025_17 vs G12D | 0.633 | 0.011 |
| rfd4_d031 vs G12D | 0.735 | 0.237 |

An ipSAE of exactly 0 means no residue pair across the interface had PAE below 10. The
chains are still in contact (min backbone gap 3.95 A for d010_13 vs WT), so Boltz places
the binder against the groove but has no confidence in where. On that reading ipTM was
crediting designs for overall fold confidence rather than for a real interface.

This conflict should be stated in the submission rather than hidden: the ranking depends
on which metric is trusted, and we trust the PAE-derived one.

## What the redesign did

ProteinMPNN with a positive-charge bias at the electrostatic band inverted the parent's
preference. Parent rfd4_d010: G12D ipSAE 0.000, WT 0.000 (no interface either way, and by
ipTM it mildly preferred WT, 0.247 vs 0.445). Redesigned mpnn_d010_13: G12D 0.833, WT 0.046.

The scaffold alone does not explain it. The sibling mpnn_d010_14, same scaffold, same MPNN
run, binds G12D slightly harder (0.864) but keeps substantial WT (0.491) and apo (0.760)
interface. Sequence, not backbone, decides the selectivity.

## Caveats

- One diffusion sample per prediction. Only d010_13 and d031_01 were run at 3 seeds.
- Boltz co-folding confidence is not a binding affinity and has not been calibrated
  against this target class.
- mpnn_d031_01 looked selective at one seed (+0.42 ipTM gap) and reversed at another
  (-0.17). Single-seed gaps on this system are not reliable.
- 2 WT complexes were flagged interpenetrating (mpnn_d025_17 at 2.20 A, rfd4_d025 at
  1.91 A). Both are non-binders by ipSAE anyway.
