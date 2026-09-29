# FoldX stability of 1XYN ↔ 1KS5 region-swap mutants

This folder estimates how much each region-swap mutant destabilizes its host protein. It uses FoldX 5.1 `BuildModel` on the crystal structures of *T. reesei* xylanase XYN I (PDB **1XYN**, GH11) and *A. niger* endoglucanase EglA / Eg1A (PDB **1KS5**, GH12). The FoldX jobs run on the Hellbender cluster (SLURM). Input preparation and analysis run on any machine with Python 3.

Each swap replaces a segment of one protein with the equivalent segment of the other. Both segments have the same length, so each mutant is modelled as a set of simultaneous point substitutions on the host crystal backbone.

Input structures, sequences and results are **not** in the repository (see `.gitignore`). Section 1 explains how to recreate them.

---

## Contents

| File | Purpose |
|---|---|
| `make_foldx_inputs.py` | Cleans the crystal PDBs and turns each swap sequence into a FoldX mutation list |
| `hpc/env.sh` | FoldX path and number of runs per mutant. **Edit for your cluster** |
| `hpc/repair.sbatch` | `RepairPDB` on 1XYN and 1KS5 (2-task array) |
| `hpc/buildmodel.sbatch` | `BuildModel`, one array task per mutant |
| `hpc/summary.sbatch` | Collects ΔΔG of every run into `foldx_summary.tsv` |
| `hpc/submit.sh` | Submits repair → BuildModel → summary with SLURM dependencies |
| `compare_mutants.py` | Breaks down and compares the results across mutants |
| `af3_segment_check.py` | Optional. Uses AlphaFold3 models to flag swaps whose segment is predicted to refold |

Requirements:
- FoldX 5.1, with an academic licence from <https://foldxsuite.crg.eu>. Do not commit the FoldX binary or `rotabase.txt`; the licence doesn't allow redistribution.
- SLURM.
- Python 3 (standard library only) for input preparation and analysis.
- `af3_segment_check.py` also needs `gemmi` and `numpy`.

---

## 1. Make the input files

### 1.1 Crystal structures → `pdb/`

```bash
cd foldx_swaps
mkdir -p pdb
curl -fL -o pdb/1XYN.pdb https://files.rcsb.org/download/1XYN.pdb
curl -fL -o pdb/1KS5.pdb https://files.rcsb.org/download/1KS5.pdb
```

The file names must be exactly `pdb/1XYN.pdb` and `pdb/1KS5.pdb`. Both are single chains (chain A), numbered 1–178 and 1–223, with no missing residues.

### 1.2 Swap sequences → FASTA

One record per mutant, giving the full-length mutant sequence:

```
>1KS5_M02 | 1KS5:30-35 replaced_by 1XYN:14-19
QTMCSQYDSASSPPYSVNQNLWGEYQGTGQVSYSPDKLSSSGASWHTEWTWSGGEG...
>1XYN_M02 | 1XYN:14-19 replaced_by 1KS5:30-35
ASINYDQNYQTGGSQCVYVSNTGFSVNWNTQDDFVVGVGWTTGSSAPINFGGSFSV...
```

Rules the script checks:
- **Name** (the text before `|`): `<HOST>_<SWAP>`. `HOST` must be `1XYN` or `1KS5`. The two directions of one swap share the same `SWAP` tag (for example `M02`); that's how reciprocal pairs are matched.
- **Sequence:** exactly as long as the host crystal structure (178 or 223 residues). Swaps that change the length can't be modelled with `BuildModel`.
- **Region in the header** (`HOST:start-end`): optional and only used for a warning. Mutations are always read from the sequence itself, by comparing it with the crystal sequence.

The original set is 10 mutants (M02, M03, M05, M08 and M13, each in both directions) in `1KS5_1XYN_region_swap_10_sequences.fasta`.

### 1.3 Build the FoldX inputs

```bash
python3 make_foldx_inputs.py path/to/1KS5_1XYN_region_swap_10_sequences.fasta
```

If you leave out the path, the script reads `../1KS5_1XYN_region_swap_10_sequences.fasta`.

The script writes:

| Output | Content |
|---|---|
| `repair/1XYN.pdb`, `repair/1KS5.pdb` | Chain A protein atoms only. Waters, the 1XYN Ca²⁺ ion (Asp81/His83, outside every swap region) and alternate positions are removed |
| `mutants/<name>/individual_list.txt` | One FoldX mutation line, e.g. `SA30Q,QA31V,CA32S,VA33Y,YA34S,VA35P;` |
| `mutants.tsv` | One row per mutant: name, host, number of substitutions, changed span, mutations. **The row order sets the SLURM array index** |

The script prints each mutant's changed span and stops with an error if:
- a length differs from the host;
- a residue is non-standard;
- the PDB has gaps or insertion codes.

It also checks each reciprocal pair: the segment inserted into one host must equal the other protein's wild-type sequence over the matching span.

> **Known header mismatch.** The M03 headers say `1KS5:44-51` ↔ `1XYN:23-30` (8 residues). The M03 sequences actually swap **1KS5:42-52 ↔ 1XYN:22-32** (11 residues). The script prints a warning for this and uses the sequences. The other four swaps match their headers.

Expected result:

| Swap | 1KS5 host (span, substitutions) | 1XYN host (span, substitutions) |
|---|---|---|
| M02 | 30–35, 6 | 14–19, 6 |
| M03 | 42–52, 10 | 22–32, 10 |
| M05 | 81–89, 9 | 48–56, 9 |
| M08 | 132–139, 6 | 90–97, 6 |
| M13 | 176–187, 12 | 138–149, 12 |

Some counts are lower than the span length because positions where the two proteins already share a residue need no substitution.

---

## 2. Run FoldX on Hellbender

### 2.1 Configure

Edit `hpc/env.sh`:

```bash
FOLDX_DIR=/home/yjx9t/data/EG1/foldx     # folder holding the FoldX binary (and rotabase.txt, if your build needs it)
FOLDX=$FOLDX_DIR/foldx_20261231          # the FoldX 5.1 executable
NRUNS=5                                  # BuildModel runs per mutant
```

If `$FOLDX_DIR/rotabase.txt` exists, the jobs link it into each working directory automatically.

The SLURM settings use partition `requeue`, as in the other Hellbender scripts. Change `#SBATCH --partition` in the `.sbatch` files if needed.

### 2.2 Copy and submit

```bash
# from your computer: copy the folder, including the inputs made in step 1
rsync -av foldx_swaps/ <user>@<hellbender>:/home/yjx9t/data/EG1/foldx_swaps/

# on Hellbender
/home/yjx9t/data/EG1/foldx/foldx_20261231 --help | head    # check the binary runs
bash /home/yjx9t/data/EG1/foldx_swaps/hpc/submit.sh          # add --account=xxx if required
```

`submit.sh` runs from `foldx_swaps/` (whichever folder you call it from) and submits three dependent jobs:

1. **`foldx_repair`** (2 tasks, about 1 minute each): `RepairPDB` makes `repair/<PDB>_Repair.pdb`.
2. **`foldx_build`** (one task per row of `mutants.tsv`, about 5–10 minutes each with 5 runs): `BuildModel --numberOfRuns=$NRUNS` in `mutants/<name>/`. It starts only if both repairs succeeded.
3. **`foldx_summary`**: writes `foldx_summary.tsv`, one row per mutant per run.

Tasks that already have output are skipped, so you can resubmit after a failure. Logs go to `hpc/logs/`.

### 2.3 Check before copying back

```bash
cd /home/yjx9t/data/EG1/foldx_swaps
column -t foldx_summary.tsv       # 10 mutants x 5 runs = 50 rows, none marked MISSING
grep -il "error\|fatal" hpc/logs/*.out repair/*.log mutants/*/buildmodel.log
```

A `MISSING` row means that mutant's BuildModel failed. The reason is in `mutants/<name>/buildmodel.log`.

### 2.4 Copy results back

```bash
rsync -av --exclude 'rotabase.txt' \
  <user>@<hellbender>:/home/yjx9t/data/EG1/foldx_swaps/ foldx_swaps/
```

What BuildModel writes in each `mutants/<name>/`:

| File | Content |
|---|---|
| `Dif_<HOST>_Repair.fxout` | **ΔΔG for each run**, split into energy terms. This is the main result |
| `Raw_<HOST>_Repair.fxout` | Total energies of the mutant and wild-type models before subtraction |
| `<HOST>_Repair_1_<r>.pdb` | Mutant model from run *r* (0 to NRUNS−1) |
| `WT_<HOST>_Repair_1_<r>.pdb` | Wild-type model rebuilt in the same run |
| `Average_<HOST>_Repair.fxout` | Average over runs. Not used: runs are reported individually |

---

## 3. Analyse

### 3.1 Optional: flag swaps that AF3 predicts to refold

```bash
/opt/anaconda3/envs/vina/bin/python af3_segment_check.py     # any Python with gemmi + numpy
```

FoldX keeps the host backbone fixed, so its ΔΔG only makes sense if the swapped segment keeps roughly the host's shape. This script uses the AF3 models of the mutants only as a **warning flag**; they are never scored.

For each mutant, it finds the AF3 job whose sequence matches the FASTA exactly. For each of that job's 5 models, it superimposes everything except the segment ±3 residues onto the host crystal, then measures how far the segment has moved (Cα RMSD over the segment ±2 residues) and its mean pLDDT.

| Flag | Rule | Meaning |
|---|---|---|
| `native_backbone` | Moved ≤ 2 Å and pLDDT ≥ 70 | The fixed-backbone ΔΔG is usable |
| `low_AF3_confidence` | pLDDT < 70 in any model | AF3 keeps the host shape but isn't confident; use with care |
| `AF3_refold` | Moved > 2 Å in every model | The segment is predicted to refold. The ΔΔG only means "doesn't fit the native backbone"; **don't trust its value** |

For reference, AF3's wild-type models differ from the crystals by only 0.4–0.8 Å over the same segments.

Paths the script expects, relative to `foldx_swaps/`:
- `../1KS5_1XYN_region_swap_10_sequences.fasta`
- AF3 jobs as `../fold_*/` or `../SwapMutants_AF3_JY/fold_*/`, either as folders or `.zip` files, each containing `*_model_[0-4].cif`
- AF3 wild-type models `../models/1xyn_apo_model.pdb` and `../models/1ks5_WT_AF3_model_0.pdb`
- `repair/<HOST>_Repair.pdb` and `mutants.tsv`

It writes `af3_segment_check.tsv`. If this file is missing, `compare_mutants.py` still runs, just without the flags.

### 3.2 Compare mutants

```bash
python3 compare_mutants.py
```

This prints three tables and writes `foldx_compare_runs.tsv`, with one row per mutant per run:

1. **Ranges for each mutant** (min–max over runs, never averaged) of:
   - **total:** FoldX ΔΔG.
   - **steric:** `Van der Waals clashes` + `torsional clash`.
   - **packing:** total − steric − disulfide.
   - **per substitution:** total ÷ number of substitutions.
   - **disulfide term.**

   It also gives a tier, the main cause (driver) and the backbone flag.
2. **Reciprocal pairs:** which host tolerates the swapped segment better. Marked "not comparable" if either direction is `AF3_refold`.
3. **Pairs that can't be told apart:** mutants on the same host whose run ranges are separated by 1 kcal/mol or less, roughly FoldX's own error. `AF3_refold` mutants are left out.

---

## 4. Interpreting the numbers

**What ΔΔG is.** In each run, FoldX puts all of a mutant's substitutions into the repaired crystal structure at once and repacks the nearby side chains. It rebuilds the wild type (`WT_*.pdb`) the same way in the same run. It reports **ΔΔG = E(mutant) − E(WT)** in kcal/mol; positive values are destabilizing. The values are not normalized: the penalties of all substitutions add up, so larger swaps tend to score higher. The per-substitution column adjusts for this.

**Energy terms.** The `Dif` file splits ΔΔG into terms that add up exactly to `total energy`. The one exception is `backbone clash`: FoldX reports it but leaves it out of the total, so the analysis leaves it out too.

| Group | Terms | Interpretation |
|---|---|---|
| **steric** | Van der Waals clashes, torsional clash | Atoms forced too close together on the fixed backbone. A real protein could shift its backbone to relieve this, so FoldX tends to **overestimate** it. The strain at that position is still real |
| **packing** | Everything else: H-bonds, polar and hydrophobic solvation, van der Waals attraction, side-chain and main-chain entropy | Chemical mismatch (buried polar group, exposed hydrophobic group, lost H-bond) that moving the backbone would not fix. The **more trustworthy** part |
| **disulfide** | disulfide | Non-zero only when a Cys in a disulfide is replaced (C32S in 1KS5_M02 breaks Cys4–Cys32) |

**Tiers.** FoldX is calibrated on single mutations (error about 1 kcal/mol), and the error grows with the number of substitutions. Read the totals in tiers rather than as exact values:

| Max ΔΔG over runs (kcal/mol) | Tier |
|---|---|
| < 2 | tolerated |
| 2–7 | moderate |
| 7–15 | strong |
| > 15 | incompatible with the host backbone as it is. The exact number isn't meaningful |

Most globular proteins are only about 5–15 kcal/mol more stable than their unfolded state, so values far above that mean "doesn't fit", not a literal destabilization.

**Limitations.**
- The backbone is fixed; mutants that refold are flagged, not modelled (section 3.1).
- ΔΔG estimates folding stability only, not expression, activity or ligand binding.
- Treat the results as a ranking and a way to find problem residues. Confirm them with experimental data (expression, Tm) or MD.

---

## Troubleshooting

| Symptom | Where to look / what to do |
|---|---|
| `FATAL: ... not found or not executable` | Fix `FOLDX_DIR` or `FOLDX` in `hpc/env.sh`, or run `chmod +x` on the binary |
| Repair task fails | Read `repair/<PDB>_repair.log`. If it complains about rotamers, put `rotabase.txt` in `$FOLDX_DIR` |
| BuildModel task fails or a `MISSING` row | Read `mutants/<name>/buildmodel.log`; resubmit with `bash hpc/submit.sh`, which skips finished tasks |
| `make_foldx_inputs.py`: length error | That swap changes the length; `BuildModel` can't model it |
| `compare_mutants.py`: file not found | Results not copied back yet; `mutants/<name>/Dif_*.fxout` is needed |
