# 1KS4 Endoglucanase MD Pipeline

Force field: AMBER ff14SB + TIP3P, built via `tleap`, converted to GROMACS
format and run with `gmx_mpi` (same convention as the XYN WT/Q4 systems —
ff14SB parameters, GROMACS engine, not AMBER's pmemd/sander).

Structure: PDB 1KS4 (A. niger GH12 endoglucanase, EglA, EC 3.2.1.4). Native disulfide
Cys4-Cys32 (crystal SG-SG 2.03 Å).

## Directory layout

```
1KS4_GROMACS/
  *.mdp                       shared min/nvt/equi/prod parameters
  tleap_1KS4.in-style scripts  (per-case, lives one level up in codes/)
  01_amber_to_gromacs.sh       AMBER prmtop/inpcrd -> GROMACS gro/top
  make_posre_backbone.sh       backbone-only restraint (for NPT)
  insert_posre_include.sh      wires posre*.itp into system.top
  submit_replicates.sh         sets up repN/ dirs and submits jobs
  run_1KS4.batch                min -> nvt -> equi -> prod for one replicate
  continue_1KS4.batch          resumes production after 24h wall-time cutoff
  run_analysis.sh              rmsd/rmsf/gyrate/dssp, tagged per case+rep
  <CASE_DIR>/                  e.g. WT/, E7/ -- one per mutant
    system.prmtop, system.inpcrd   (from tleap)
    system.gro, system.top         (from 01_amber_to_gromacs.sh)
    posre.itp, posre_backbone.itp
    rep1/, rep2/, rep3/            (from submit_replicates.sh)
  unused_charmm36m_route/      retired pdb2gmx/CHARMM36m attempt, kept for reference
```

Scripts are run from `1KS4_GROMACS/` with `<CASE_DIR>` as an argument
(e.g. `./01_amber_to_gromacs.sh WT`), not from inside the case directory.

## Step-by-step build (per case)

### 1. Strip ions/waters from the source structure

Remove crystallographic Pd2+ ions and waters, keep protein ATOM records only,
to get a clean starting structure (`input.pdb`).

### 2. Add hydrogens

```bash
module load amber/22
reduce -BUILD input.pdb > 1ks4_ready.pdb 2> reduce_info.log
```
Note: the positional-argument form (`reduce -BUILD input.pdb > out`) is
required — stdin redirection (`reduce -BUILD < input.pdb`) silently produced
only a usage message on this AMBER/22 build.

`reduce` also distinguishes free-thiol vs. disulfide-bonded cysteines: bonded
Cys residues come out with no `HG` atom. Confirmed for residues 4 and 32.

### 3. Fix histidine tautomers

`tleap` defaults bare `HIS` to the `HIE` template regardless of which
hydrogen `reduce` actually placed. Check every His residue for which
hydrogen it has (HE2 -> HIE, HD1 -> HID, both -> HIP) and rename any that
don't match the HIE default. For 1KS4 WT: residue 46 has HE2 (matches
default, no change needed logically, but was still renamed explicitly for
clarity), residue 108 has HD1 (needed HID).

```bash
sed 's/HIS A  46/HIE A  46/; s/HIS A 108/HID A 108/' 1ks4_ready.pdb > 1ks4_ready.pdb.fixed2
```
Always confirm column spacing against the actual PDB (`grep -n`) before
running — get the literal match string right rather than guessing, and
verify with `diff`/`wc -c` (file size should be unchanged) plus an
atom-level `awk` check afterward.

### 4. Fix disulfide-bonded cysteine naming

`bond mol.X.SG mol.Y.SG` in the tleap script creates connectivity but does
**not** retroactively reclassify the residue template — tleap picks the
template (CYS = free thiol vs. CYX = disulfide-bonded) at `loadpdb` time
based on the literal residue name in the PDB. Both must say `CYX` for
residues actually in a disulfide bond, or tleap throws SH-SH/HS-SH-SH
bond/angle/torsion parameter errors partway through.

```bash
sed 's/CYS A   4/CYX A   4/; s/CYS A  32/CYX A  32/' 1ks4_ready.pdb > 1ks4_ready.pdb.fixed3
```
Watch for incidental matches in `LINK`/`SITE` header records (documentary
only, harmless, but will inflate a naive `grep -c` count) — confirm with
`grep -n` which lines are actually `ATOM` records before trusting a count.

### 5. Build the AMBER topology

```bash
tleap -f tleap_1KS4.in
```
`tleap_1KS4.in`:
```
source leaprc.protein.ff14SB
source leaprc.water.tip3p
mol = loadpdb 1ks4_ready.pdb
bond mol.4.SG mol.32.SG
solvateOct mol TIP3PBOX 10.0
addIons mol Na+ 0
addIons mol Cl- 0
saveamberparm mol system.prmtop system.inpcrd
quit
```
Check the log: `grep -iE "error|fatal" tleap_1KS4.log` should be empty, and
the log should report `Errors = 0`. Verify the disulfide geometry survived
correctly:
```bash
ambpdb -p system.prmtop -c system.inpcrd | grep CYX
```
(For WT: SG-SG distance recomputed as 2.026 Å, matching the crystal's 2.03 Å.)

### 6. Convert AMBER topology to GROMACS format

```bash
cd 1KS4_GROMACS
./01_amber_to_gromacs.sh <CASE_DIR>
```
Internally uses ParmEd via `amber.python` (bundled with the AMBER module,
no separate acpype dependency) to convert `system.prmtop`/`system.inpcrd`
directly to `system.top`/`system.gro`, preserving all ff14SB parameters
including the CYX disulfide. Also generates the heavy-atom restraint file
(`posre.itp`, since ParmEd's conversion — unlike `pdb2gmx` — doesn't produce
one automatically) and wires it into `system.top` under `#ifdef POSRES`.

Checkpoint:
```bash
ls -la <CASE_DIR>/system.gro <CASE_DIR>/system.top <CASE_DIR>/posre.itp
grep -A2 "#ifdef POSRES" <CASE_DIR>/system.top
```

### 7. Build the backbone restraint (for NPT)

```bash
./make_posre_backbone.sh <CASE_DIR>
```
Generates `posre_backbone.itp` and wires it in under `#ifdef POSRES_BB`.

Checkpoint: `system.top` should now show **both** `#ifdef` blocks
(`POSRES`->`posre.itp`, `POSRES_BB`->`posre_backbone.itp`).

### 8. Launch replicates

```bash
./submit_replicates.sh <CASE_DIR> 3
```
Creates `<CASE_DIR>/rep1,2,3/`, copies `system.gro`/`system.top`/restraint
files into each, and submits `run_1KS4.batch` per replicate. Checks that
both restraint files exist before submitting.

Checkpoint: `squeue -u $USER` should show 3 jobs.

## Simulation protocol (mdp files)

| Stage | File | Length | Restraint | Notes |
|---|---|---|---|---|
| Minimization | `min.mdp` | — | none | steepest descent, `DispCorr = EnerPres` |
| NVT | `nvt.mdp` | 500 ps | heavy-atom (`-DPOSRES`) | `gen_vel = yes`, `gen_seed = -1` (independent velocities per replicate), `pcoupl = no` |
| NPT equilibration | `equi.mdp` | 2 ns | backbone only (`-DPOSRES_BB`) | `pcoupl = C-rescale` (gentler than Parrinello-Rahman for equilibration) |
| Production | `prod.mdp` | 200 ns | none | `pcoupl = Parrinello-Rahman`, `nstxout-compressed = 5000` (10 ps/frame), `DispCorr = EnerPres` |

Each replicate runs its own full min -> nvt -> equi -> prod chain
(`run_1KS4.batch`) so that `gen_seed = -1` gives genuinely independent
trajectories, not just re-launches from the same equilibrated state.

Production restart after the 24h SLURM wall-time limit:
```bash
cd <CASE_DIR>/repN
sbatch ../../continue_1KS4.batch    # -cpi prod_1.cpt -append
```

## Analysis

```bash
./run_analysis.sh <CASE_DIR> <REP>
```
Runs `gmx_mpi trjconv`/`rmsd`/`rmsf`/`gyrate`/`do_dssp` selecting groups by
**name** ("Protein", "Backbone") rather than index, since group numbering
isn't guaranteed to match between a `pdb2gmx`-built topology and this
ParmEd-converted one. Output files are tagged
`1KS4_<CASE_DIR>_rep<N>_*.xvg` to avoid collisions with other projects'
identically-named output files.

To characterize a specific salt bridge or contact of interest (e.g. near a
mutation site) without biasing the simulation, measure it post-hoc on the
production trajectory rather than restraining it during the run:
```bash
gmx_mpi mindist -f prod_1.xtc -s prod_1.tpr -n index.ndx -od saltbridge_dist.xvg
```
(select the two relevant atom groups, e.g. the charged side-chain atoms).

## Status

- **WT**: fully built and verified (tleap `Errors = 0`, disulfide geometry
  confirmed, restraints wired). 3 replicates submitted and queued.
- No other 1KS4 mutants defined yet.
