# Inputs

Everything the test system is built from, so a fresh clone reproduces it.

| file | what it is |
|---|---|
| `receptors/fold_m02_1ks5backbone_1xyninsert_model_0.pdb` | AlphaFold3 model: 1ks5 (Eg1A) backbone with the 1xyn insert, chain A, 223 residues |
| `receptors/..._model_0_clean.pdb` | the same file after `prep_receptor.py` — this is what the build uses |
| `poses/..._model_0_X6_ad4_out.pdbqt` | AutoDock4 poses of xylohexaose in that model, 9 poses, best first |
| `models_docking_summary.tsv` | best affinity for all 24 (model, ligand) pairs from the docking run |
| `box_centers.tsv` | docking box centres |

Only this one receptor and ligand are here on purpose. To add more, copy the
model PDB into `receptors/`, its `*_ad4_out.pdbqt` into `poses/`, and run
`./new_system.sh` again — nothing else changes.

## The ligand

**X6** = xylohexaose, 55 heavy atoms. Its CHARMM topology is ready in
`../ligands/X6_charmm/`.

**cellose** (cellobiose) has **no CHARMM topology yet** — see
`../ligands/README.md` before trying to build a cellose system.

## Which pose to use

The top pose is not automatically usable. AutoDock rotates glycosidic torsions
without sugar-specific torsional preferences, so a pose can come back eclipsed
and fail to build in CHARMM.

For this receptor, **pose 1 (-9.477 kcal/mol) builds cleanly** and is what
`systems/m02_X6.conf` uses. Poses 2, 3 and 7 are refused. Check any new pose
file before committing to it:

```bash
./check_poses.sh <system> inputs/poses/<file>.pdbqt
```
