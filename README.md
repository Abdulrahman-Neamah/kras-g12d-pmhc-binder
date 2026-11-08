# kras-g12d-pmhc-binder

De novo binder design against the **KRAS G12D neoantigen peptide VVVGADGVGK presented on
HLA-A\*11:01**, built for Challenge 5 of the Anthropic × Adaptyv protein design competition.

The hard part is not binding. It is *discriminating*: the design must bind the mutant pMHC,
not the wild-type peptide VVVGAGGVGK on the same MHC, and not the empty HLA/β2m
heterodimer. All three are measured under identical conditions.

**Result: one design, `mpnn_d010_13`, discriminates on both counter-targets across three
independent seeds.** Its parent scaffold had no detectable interface at all — a
charge-biased ProteinMPNN redesign created one.

| | ipSAE (binder–peptide) |
|---|---|
| G12D | **0.833** |
| wild-type peptide | 0.046 |
| apo HLA/β2m (binder–HLA) | **0.000** |

---

## The central finding: a metric disagreement that changes the answer

Boltz-2 reports ipTM from a learned prediction head. ipSAE is computed directly from the
PAE matrix with a strict per-residue-pair cutoff. On this system they disagree sharply:

| design / target | ipTM | ipSAE |
|---|---|---|
| mpnn_d010_13 vs wild type | 0.52 – 0.73 | 0.000 – 0.137 |
| mpnn_d002_09 vs G12D | 0.685 | 0.123 |
| mpnn_d025_17 vs G12D | 0.633 | 0.019 |

An ipSAE of exactly 0 means **no residue pair across the interface had PAE below 10**. The
chains are still in contact — 3.95 Å minimum backbone separation for `d010_13` against wild
type — so the model places the binder against the groove but has no confidence in where.

By ipTM, the initial screen of 20 designs looked like a usable field. By ipSAE, **15 of
those 20 have no confident interface to any target, including the one they were designed
for.** They are not weak binders; they are non-binders that ipTM flattered. Across all 27
designs eventually screened, the same 15 remain without an interface — every design added
later came from the one scaffold that worked.

This is not a tie-break. Ranking on ipTM alone would have missed `mpnn_d010_16` entirely
(ipTM wild-type gap +0.10, ipSAE +0.45) and promoted several designs with no interface.

## The epitope was measured, not assumed

The mutation sits at peptide position 6. The intuitive move is to contact the mutant
aspartate directly. Measured on 7OW6 (JDIa41b1 TCR on HLA-A\*11:01/VVVGADGVGK), that is the
wrong move:

| quantity | value |
|---|---|
| Asp6 side chain → nearest TCR atom | 5.73 Å |
| TCR atoms within 4.0 / 5.0 / 6.0 Å | 0 / 0 / 1 |
| Asp6 side-chain SASA buried **by the TCR** | **0.0 Å²** |
| Asp6 side-chain SASA buried **by the HLA** | 76.4 Å² |
| Asp6 OD2 → HLA Arg114 NH2 | 2.95 Å (bidentate salt bridge) |

The receptor never touches it. The mutant carboxylate is buried by the MHC and read
*indirectly*, through the electrostatic field it projects onto the pMHC surface. The same
holds in 7PB2, 6ULN, 6ULR, 6UON, 7SU9 and 9WJ6.

Two further reasons to leave it alone: in the free G12D pMHC structure (7OW4, 1.81 Å)
peptide residues 5 and 6 are **not modelled in any of the four copies** — the state a de
novo binder first meets has no ordered Asp6 — and contacting it means competing with a
2.95 Å salt bridge, destroying the feature you are trying to read.

`analysis/01_position6_readout.py` reproduces these measurements.

## What the design step actually did

ProteinMPNN was run with a **per-residue positive-charge bias** on residues whose CB lies
within 9.0 Å of the point where the G12D-minus-wild-type electrostatic potential difference
is largest, with existing charge-reader residues held fixed.

The motivation is a measured failure mode: unbiased ProteinMPNN places more acidic than
basic residues at that field optimum and destroys existing readers 92% of the time.

It worked, and the effect is not subtle:

| | ipSAE G12D | ipSAE wild type |
|---|---|---|
| parent `rfd4_design_010` | 0.000 | 0.000 |
| redesigned `mpnn_d010_13` | 0.833 | 0.046 |

The scaffold's approach geometry already resembled the natural receptor's — 160.4 Å² of
peptide buried against the TCR's 169 Å², footprint cosine similarity 0.872 — but it had no
interface until the charge bias was applied.

## Replication caught two false positives

Five designs were run at multiple seeds. Two would otherwise have entered the submission:

| design | per-seed wild-type gap | verdict |
|---|---|---|
| mpnn_d031_01 | +0.42, −0.17, 0.00 | **sign flips** — not a discriminator |
| mpnn_d010_16 | +0.45, +0.53, *failed to fold* | **unstable** — no G12D complex at 1 of 3 seeds |
| mpnn_d010_13 | +0.82, +0.86, +0.68 | holds |

Single-seed discrimination gaps on this system are not reliable.

## Nine variants, one scaffold, and what that taught

Eight of the nine `d010` variants bind G12D strongly (ipSAE 0.786–0.883). The ninth,
`_16`, averages 0.538 only because one of its three seeds failed to form the complex at
all; its two working seeds sit in the same band. The variance between the variants is
otherwise entirely in the counter-targets:

- 1 discriminates cleanly on both (`_13`)
- 2 specialise in one counter-target each (`_25` on wild type, `_21` on apo)
- 6 discriminate on neither

**The scaffold solves binding; the sequence solves selectivity.** Which also means
selectivity could not be predicted in advance from the backbone — it had to be screened for.

## A clash check, because the last round needed one

In Challenge 4 the top-ranked design turned out to be threaded through its target: 0.41 Å
between the closest backbone atom pair. ipTM, ipSAE and pDockQ2 all scored it well, because
**none of them contains a steric term**.

`scripts/10_interchain_backbone_clash.py` is a few dozen lines and would have caught it.
Every design here was checked. Two wild-type complexes failed (2.20 Å, 1.91 Å); no
submitted design clashes against G12D.

## Limitations

- One diffusion sample per prediction; only five designs replicated across seeds.
- Boltz-2 co-folding confidence is **not** binding affinity and has not been calibrated
  against pMHC targets.
- Seven of the ten submitted designs are one scaffold family at 70–78% identity. If the
  backbone is wrong, the submission fails as a block.
- Whether the assay measures wild-type binding at all, or only G12D binding, could not be
  determined from the published materials.

## Repository contents

```
analysis/
  01_position6_readout.py      Asp6 contact, SASA and salt-bridge measurements on 7OW5/7OW6
  02_tcr_footprint.py          per-position TCR burial on the peptide
  ipsae_analysis.md            full validation writeup
scripts/
  10_interchain_backbone_clash.py   steric check no confidence metric performs
  build_submission.py          ranks by the worse of the two discrimination gaps
notebooks/
  boltz2_validation.ipynb      three-target co-folding, resumable, with ipSAE scoring
results/
  design_metrics.csv           pLDDT, ipSAE, hydropathy, charge, clash gap per design
  ipsae_scores.csv             ipSAE and pDockQ2, both asymmetric directions, all runs
  boltz2_raw_results.csv       raw Boltz-2 output, 84 predictions
```

### Running these

The analysis scripts need 7OW5.pdb and 7OW6.pdb from the RCSB PDB:
`STRUCT_DIR=/path/to/pdbs python3 analysis/01_position6_readout.py`.
`build_submission.py` needs the design-sequence CSVs, which are withheld (see above);
point `SEQ_DIR` at your own. The clash check takes its inputs as arguments and needs
nothing else. Dependencies: Biopython, NumPy, SciPy.

Ranking uses the **worse** of the two discrimination gaps, not the average: a design has to
clear both goals, and averaging lets a strong apo gap hide a failed wild-type gap.

## Design sequences

Withheld pending competition results, as with the earlier challenges in this series. They
will be added here, with the lab measurements, once published.

## Methods

RFD4-Proteina backbone generation (Tamarind Bio) → motif and composition screening →
charge-biased ProteinMPNN redesign (`--use_soluble_model --omit_AAs CX --sampling_temp 0.1
0.2 --num_seq_per_target 32`) → Boltz-2 2.0.3 co-folding against all three targets
(`--diffusion_samples 1 --recycling_steps 3`) → ipSAE (PAE cutoff 10, distance cutoff 10,
[DunbrackLab/IPSAE](https://github.com/DunbrackLab/IPSAE)) → backbone clash check.

Chains throughout: A = HLA-A\*11:01, B = β2-microglobulin, C = peptide, Z = binder.
