# Ligand topologies

CHARMM36's GROMACS port ships single sugars (`carb.rtp`) but **no 1->4 glycosidic
patch**, so `pdb2gmx` cannot build any oligosaccharide. Every ligand here needs a
topology generated once by CHARMM-GUI Glycan Reader & Modeler, then reused for
every pose of that ligand.

| ligand | what it is | topology | status |
|---|---|---|---|
| `X6` | xylohexaose, 105 atoms (C30H50O25) | `X6_charmm/X6.itp` | **ready** — reused from XYN_charmm_MD |
| `cellose` | cellobiose, 45 atoms (C12H22O11) | `cellose_charmm/cellose.itp` | **you must generate it** (see below) |

`X6.itp` is the CHARMM-GUI `CARB.itp` from the XYN project, unchanged. The
molecule is identical in every docked pose, so only the coordinates differ, and
`02_build_glycan.py` is what puts the topology into a given pose.

## Generating a missing ligand topology

1. https://charmm-gui.org -> Input Generator -> **Glycan Reader & Modeler**.
2. Build the glycan (cellobiose = two beta-D-glucose joined 1->4), solvate with
   defaults, and request **GROMACS** output.
3. Download the archive and take `gromacs/toppar/CARB.itp` from it.
4. Put it here as `<ligand>_charmm/<ligand>.itp`, and the glycan's own PDB as
   `<ligand>_charmm/<ligand>_charmm_reference.pdb`.
5. Check the atom count matches the docked ligand's heavy-atom count plus its
   hydrogens; `prep_system.sh` will refuse to build if they disagree.
