# EG1 domain-swap complexes — CHARMM36m MD on the WashU engineering cluster

Molecular dynamics for the EG1 domain-swap models and their docked ligands,
using **exactly the simulation protocol from `XYN_charmm_MD`** so the two sets of
runs can be compared directly, and running **production on a GPU**.

This document assumes you have never used this cluster. Every command is meant
to be copied as written, except where it says to fill something in.

---

## 0. What you get

One test system — **`m02_X6`**: the m02 domain-swap model (1ks5 backbone,
1xyn insert) with xylohexaose in its top docked pose. Per system, **two**
independent replicates of:

| stage | length | restraints | purpose |
|---|---|---|---|
| minimisation | ≤50 000 steps | heavy atoms | remove build clashes |
| NVT | 500 ps | heavy atoms | reach 300 K |
| NPT | 2 ns | backbone only | reach the right density; ligand free to settle |
| production | 100 ns | none | the data |

2 fs steps, frames every 10 ps → 10 001 frames per replicate. CHARMM36m with
CHARMM-modified TIP3P, 150 mM NaCl, dodecahedral box with 1.2 nm padding, 300 K.

**Replicates are independent experiments and are reported individually.** Nothing
here averages them, and neither should your figures: a mean across runs hides
where they disagree, which is usually the interesting part.

Set `N_REPS` in `config.sh` to run more.

---

## 1. Get onto the cluster

```bash
ssh <netID>@ssh-shell-3.engr.wustl.edu
```

Your password will not show as you type. That is normal.

**Where to put things.** Your home directory is capped at ~10 GB and is for
scripts and config only. All simulation data goes in the project space:

```bash
cd /engrfs/project/joshua.yuan/Haina
mkdir -p $USER && cd $USER
```

A 100 ns replicate is a couple of GB. Three replicates × several systems will not fit
in home, and filling home breaks your login in confusing ways.

**Never run simulations on the login node.** It is for editing, copying and
submitting. Work runs through SLURM.

---

## 2. One-time setup

### 2a. What has to be installed

Almost nothing. The pipeline deliberately depends only on what a login shell
already has:

| needed | notes |
|---|---|
| `bash` | the scripts are plain bash, no arrays-of-arrays or GNU-only flags |
| `python3` | **standard library only** — `os`, `sys`, `re`, `math`, `glob`. No numpy, scipy, rdkit, MDAnalysis. Nothing to `pip install`, no conda environment, no `module load python` unless `python3` is missing entirely |
| `apptainer` (or `singularity`) | to run the GROMACS container |
| SLURM | `sbatch`, `squeue`, `scontrol` |

Check in one line:

```bash
python3 -c "import os,sys,re,math,glob; print('python', sys.version.split()[0], 'ok')"
```

Python 3.6 or newer (the scripts use f-strings). Anything on a current cluster
is far newer than that.

### 2b. Copy this pipeline over

From your laptop:

```bash
git clone <your-repo-url>            # or: scp -r EG1_charmm_MD <netID>@ssh-shell-3.engr.wustl.edu:...
cd EG1_charmm_MD
```

### 2c. Check your account

`config.sh` is set to `engr-lab-joshua.yuan`, taken from the lab's working GPU
job (`Docking/run_GNINA_4HK8_validation.slurm`). If you submit under a different
PI, change it. Confirm with:

```bash
sacctmgr -n show assoc user=$USER format=account%30
```

### 2d. GPU requests here are SHARDS, not whole GPUs

This is the one piece of syntax that catches people. The usual
`--gres=gpu:1` is **not** what this cluster uses. The lab's working job asks for

```
#SBATCH --gres="shard:a6000:1"      # the lab's GNINA example
#SBATCH --gres="shard:a100-sxm4:1"  # what this pipeline asks for
```

A *shard* is a fraction of a physical GPU. A plain `--gres=gpu:1` may be
rejected outright, or sit in the queue forever, if the partition only advertises
shards. Confirm that `a100-sxm4` shards really are offered before submitting a
batch:

```bash
sinfo -p general-gpu -o "%P %G %l"
```

If the shard type is spelled differently there, change `GRES` in `config.sh` to
match it exactly. One shard is enough for one GROMACS replicate; ask for more
only if a run turns out to be GPU-bound.

### No walltime is requested

These jobs pass **no `--time`**, so each takes the partition default. The chain
does not assume what that is: at startup `run_md.batch` asks the scheduler what
this job was actually granted and sets mdrun `-maxh` half an hour below it.

That margin is the point. If mdrun is still running when SLURM kills the job,
the script dies with it, the self-resubmit never happens, and you come back to a
replicate that simply stopped with nothing in any log explaining why. Stopping
early instead means a clean checkpoint and a new job.

If the scheduler reports no limit at all, `MAXH` in `config.sh` (23.5 h) is the
fallback. Being wrong there costs at most 15 minutes: `-cpt 15` checkpoints that
often, and the next link resumes from the checkpoint.

### 2e. Check GROMACS runs

GROMACS on this cluster is **a container, not a module**. `gmx` alone is not on
the container's PATH — that is the error the lab wiki documents:

```
FATAL: "gmx": executable file not found in $PATH
```

The container holds one build per CPU instruction set, and you must name the
full path. Check which set your nodes have:

```bash
lscpu | grep -o 'avx512[a-z]*' | head -1        # avx512f -> use avx_512
```

Then confirm the container runs:

```bash
apptainer exec /engrfs/project/joshua.yuan/Haina/software/containers/gromacs.sif \
  /usr/local/gromacs/avx_512/bin/gmx --version | head -20
```

Look at two lines of the output:

- **GROMACS version** — record it in your methods.
- **GPU support** — if this says `disabled`, the container has no CUDA build and
  production will silently run on CPU. Stop and ask the cluster admins, or set
  `GPU_FLAGS=""` in `config.sh` and accept CPU speed.

You do not type this path anywhere else. `gmx.sh` builds it from `config.sh`
(`CONTAINER`, `GMX_ARCH`) and every script uses that.

### 2f. Unpack the force field

```bash
./01_get_charmm36m.sh
```

This puts `charmm36-jul2022.ff/` next to the scripts, where `pdb2gmx` will find
it via `GMXLIB`.

---

## 3. Clean the receptor PDB

`pdb2gmx` will not build a raw model file reliably, and the ways it fails are
quiet ones. Clean it first:

```bash
python3 prep_receptor.py inputs/receptors/<model>.pdb \
                         inputs/receptors/<model>_clean.pdb --chain A
```

What it removes, and why each matters:

| removed | why |
|---|---|
| `HETATM` records | waters, ions and the docked ligand. The ligand arrives separately, in CHARMM naming, via `02_build_glycan.py` |
| all but one chain | a second copy doubles the system and renumbers everything after it |
| alternate locations | `pdb2gmx` would build both copies of the side chain |
| hydrogens | `pdb2gmx` rebuilds them for the force field actually in use; keeping the model's own is how naming mismatches get in |

What it **refuses** to do rather than guess:

- a **non-standard residue** — nothing in CHARMM36 to build it with;
- a residue with **no CA** — a missing backbone atom fails much later and much
  less clearly;
- it warns on **numbering gaps**, because `pdb2gmx` joins them as if continuous.
  That is right only if the gap is renumbering, not missing residues.

Those checks are the point of the script. A model that is already clean passes
through unchanged — the m02 model here does:

```
residues  : 223  (1-223)
atoms kept: 1717   (dropped 0 HETATM, 0 hydrogens)
```

`new_system.sh` runs exactly this step for you, so the standalone call above is
only needed when you want the `_clean.pdb` as a separate file to inspect or
reuse.

---

## 4. Set up a system

**This is the only step where you name an input file.** Nothing downstream is
edited by hand.

The system in this repository was made with:

```bash
./new_system.sh m02_X6 \
   --receptor inputs/receptors/fold_m02_1ks5backbone_1xyninsert_model_0_clean.pdb \
   --ligand X6 \
   --pose inputs/poses/fold_m02_1ks5backbone_1xyninsert_model_0_X6_ad4_out.pdbqt \
   --pose-index 1
```

Apo (protein only) is the same without `--ligand`/`--pose`:

```bash
./new_system.sh m02_apo \
   --receptor inputs/receptors/fold_m02_1ks5backbone_1xyninsert_model_0_clean.pdb
```

Options: `--chain B`, `--pose-index 3`, `--force` to replace an existing system.

This writes `systems/<name>.conf` and `build/<name>/protein.pdb` (plus
`ligand_af3.pdb` for holo), and prints the residue count, histidine count and
net charge it derived from your structure. Those are **expectations**, and
`prep_system.sh` checks them against what `pdb2gmx` actually builds. If they
disagree, the structure is not what the file says it is — find out why before
continuing.

> `EXPECT_DISULFIDES` is written as `0`. **If your protein has disulfides, edit
> the conf now.** The build stops if `pdb2gmx` finds a different number.

To add a system later, just run `new_system.sh` again with a new name. There is
no list to update: "all systems" means every file in `systems/`.

---

## 5. Build the ligand topology (holo only)

CHARMM36's GROMACS port ships single sugars but **no glycosidic patch**, so
`pdb2gmx` physically cannot build an oligosaccharide. The topology comes from
CHARMM-GUI once per molecule, and is then reused for every pose.

- **X6 (xylohexaose)** — already done, in `ligands/X6_charmm/`, reused unchanged
  from the XYN project.
- **cellose (cellobiose)** — **not yet generated.** See `ligands/README.md`.
  Holo-cellose systems cannot be built until you make it.

Put the docked pose into that topology:

```bash
LIGAND=X6 python3 02_build_glycan.py eg1_m08_X6
```

The script re-derives every atom name from the connectivity, checks all 24
stereocentres, checks that every bond the topology declares really is a bond in
this pose, and refuses if any of that fails.

### If it says REFUSED

A common refusal is:

```
*** REFUSED: 1 placed hydrogens are within 1.20 A of a non-bonded atom
```

This is about **one pose**, not your docking run. AutoDock rotates the
glycosidic torsions without sugar-specific torsional preferences, so some poses
come back eclipsed, and two ring hydrogens on neighbouring sugars end up on top
of each other once hydrogens are built. See which poses are usable:

```bash
./check_poses.sh m02_X6 inputs/poses/fold_m02_1ks5backbone_1xyninsert_model_0_X6_ad4_out.pdbqt
```

For the m02 receptor here, **pose 1 (−9.477 kcal/mol) builds cleanly** and is
what `systems/m02_X6.conf` uses; poses 2, 3 and 7 are refused. It is not always
the top pose: on the 1ks5 wild-type docking, 5 of 8 were refused and the best
usable one was pose 2, 0.06 kcal/mol behind pose 1 — far below the method's
accuracy. Take the clean pose; do not lower the threshold to force a strained
geometry through.

---

## 6. Build the solvated system

```bash
./prep_system.sh m02_X6
```

Runs `pdb2gmx` → restraints → box → solvate → ions, and finishes with a real
`grompp` of `nvt.mdp` so a topology that cannot run is caught now rather than
after the job queues. Output lands in `build/<name>/system.{gro,top}`.

This is CPU work and takes a few minutes. It is small enough for an interactive
session:

```bash
srun -p general-cpu -c 4 -A <your-account> -J prep --pty /bin/bash
```

---

## 7. Run

```bash
./submit.sh m02_X6              # this system, 2 replicates
./submit.sh                     # every system in systems/
```

Each replicate gets `runs/<system>/rep<N>/` and its own GPU job. 100 ns does not
fit in one job and is not meant to: each job runs until its wall limit,
checkpoints, and **submits its own continuation** (up to `MAXCHAIN`). You submit
once.

`gen_seed = -1` gives each replicate different starting velocities, which is
what makes them independent.

### Watching it

```bash
squeue -u $USER
tail -f runs/m02_X6/rep1/gmx_*.log      # stdout+stderr, merged on purpose
grep -h "Performance" runs/*/rep1/gmx_*.log     # ns/day
grep -h "ALL STAGES COMPLETE" runs/*/rep*/gmx_*.log
```

Check the first job's log early for these two lines:

```
GPU support: enabled
1 GPU selected for this run
```

If mdrun reports no GPU, it is running on CPU at perhaps a tenth of the speed.
The usual cause is a missing `--nv` on `apptainer exec`, which `gmx.sh` passes —
so if it happens, check that the node actually has a GPU allocated (`nvidia-smi`
inside the job).

---

## 8. Analyse

```bash
./submit_analysis.sh m02_X6
```

Per replicate: backbone RMSD, per-residue backbone RMSF, and **per-residue amide
hydrogen RMSF** for comparison with HDX. Outputs are tagged
`m02_X6_rep1_c36m`, `m02_X6_rep2_c36m`. Details in the header of `run_analysis.sh`.

---

## 9. Changing the protocol

Everything is in `config.sh` — box, salt, temperature, stage lengths, replicate
count, walltime, GPU offload. The `mdp/` files read those numbers, so change
them in one place.

**If you change anything in PART 2 of `config.sh`, the runs are no longer
comparable with `XYN_charmm_MD`.** That is a decision, not a detail.

`GPU_UPDATE=1` moves the integrator onto the GPU as well. It is off by default
because some builds and ensembles refuse it. Turn it on, run one replicate, and
read what mdrun says at startup before trusting it.

---

## 10. Troubleshooting

| symptom | cause | fix |
|---|---|---|
| `"gmx": executable file not found in $PATH` | container has no bare `gmx` | already handled by `gmx.sh`; if you ran the command by hand, use the full `/usr/local/gromacs/<arch>/bin/gmx` path |
| `Illegal instruction` from gmx | wrong `GMX_ARCH` for the node | check `lscpu`, set `GMX_ARCH` in `config.sh` |
| mdrun runs but reports no GPU | no `--nv`, or no GPU allocated | confirm the shard reached the job: `scontrol show job <id> \| grep -i gres`, then `nvidia-smi` inside the job |
| job pending forever, or `Invalid generic resource` | asked for `gpu:1` where the partition offers shards | `sinfo -p general-gpu -o "%P %G"`, set `GRES` in `config.sh` |
| every job rejected on walltime | `WALL_H` exceeds the partition limit | `sinfo -p general-gpu -o "%P %l"`, lower `WALL_H` |
| `ACCOUNT is still the placeholder` | step 2b not done | set `ACCOUNT` in `config.sh` |
| `REFUSED: placed hydrogens ... within 1.20 A` | eclipsed glycosidic torsion in that pose | `./check_poses.sh`, use a clean pose |
| `no CHARMM topology for 'cellose'` | topology not generated yet | `ligands/README.md` |
| build stops on residue/charge/histidine mismatch | structure ≠ what the conf claims | do not "fix" the conf until you know why they differ |
| job dies instantly with a Spack/module error | not applicable here — this cluster uses the container, no modules | — |

---

## 11. What has and has not been tested

Tested on a laptop, without GROMACS:

- cleaning, pose extraction and registration of the real `m02_X6` system
  (223 residues, 55 heavy atoms, −9.477 kcal/mol matching
  `models_docking_summary.tsv`);
- the CHARMM glycan build for all 9 of its X6 poses (6 clean, 3 refused), with
  the top pose clean and in use;
- every script parses (`bash -n`, `ast.parse`);
- **the whole pipeline run from a copied-out directory with GROMACS and SLURM
  stubbed**, to catch anything that only breaks after transfer: `run_md.batch`
  resolved its root, system and replicate from the directory and ran
  min -> NVT -> NPT -> production; `run_analysis.sh` produced all three outputs;
  `submit.sh` staged both replicates and passed the right sbatch options
  (`--gres=shard:a100-sxm4:1`, the account, no `--time`);
- no absolute paths to any local machine, and no macOS-only shell idioms
  (`sed -i ''`, `readlink -f`, `stat -f`).

**Not tested**, because it needs the cluster:

- the container path, `GMX_ARCH`, and whether that build has GPU support;
- that GROMACS works under `--gres=shard:a100-sxm4:1`. The partition, account and
  shard syntax come from a job that really runs on this cluster, but that job was
  GNINA on an **a6000** shard, not GROMACS on an **a100-sxm4** one. Confirm the
  shard type exists (`sinfo -p general-gpu -o "%P %G"`) before a large batch;
- `pdb2gmx`, solvation, ions, and the `grompp` dry run in `prep_system.sh`;
- mdrun performance and the self-resubmitting chain.

Do step 2e and one `prep_system.sh` before submitting a large batch.

---

## Repository layout and what not to commit

This folder is self-contained: `inputs/` holds every structure and docked pose,
so a fresh clone can rebuild every system without reaching outside the repo.

`.gitignore` keeps simulation output out of git — `runs/`, trajectories,
checkpoints, and the solvated `system.gro`/`system.top`, all of which are
regenerable and far larger than the rest of the repository. What IS committed is
everything needed to reproduce them: inputs, `systems/*.conf`, the registered
`build/<system>/protein.pdb` and glycan files, and the scripts.

Only the m02 receptor and its X6 poses are included, since that is the system
being tested. Adding another is two file copies and one `new_system.sh` call —
see `inputs/README.md`.

---

## File map

```
config.sh            every setting: cluster (PART 1) + protocol (PART 2)
gmx.sh               resolves gmx inside the Apptainer container (sourced)
new_system.sh        register a system from any PDB (+ optional docked pose)
prep_receptor.py       clean a receptor into build/<sys>/protein.pdb
pose_from_docking.py   extract one docked pose's heavy atoms
check_poses.sh       which docked poses can be built in CHARMM
02_build_glycan.py   docked pose -> CHARMM glycan topology
merge_holo.py          merge glycan into the protein .gro/.top
prep_system.sh       pdb2gmx -> box -> solvate -> ions -> grompp test
submit.sh            submit replicates (GPU)
run_md.batch           the chain: min -> NVT -> NPT -> production
submit_analysis.sh   submit analysis
run_analysis.sh        RMSD + RMSF + amide-H RMSF per replicate
01_get_charmm36m.sh  unpack CHARMM36m
mdp/                 the four .mdp files (identical to XYN_charmm_MD)
inputs/              receptors/ (m02 model, raw + _clean), poses/ (X6 .pdbqt), summary tsv
ligands/             per-ligand CHARMM topologies (X6 ready; cellose not)
systems/             one .conf per system — written by new_system.sh
build/               prepared systems
runs/                <system>/rep<N>/ — the simulations
```
