#!/bin/bash
# Converts the AMBER (ff14SB, via tleap_1KS4.in) topology into GROMACS format
# so the existing 1KS4_GROMACS/*.mdp + run_1KS4.batch pipeline can run it with
# gmx_mpi -- matches how the XYN WT/Q4 systems were actually built (ff14SB via
# tleap) and simulated (GROMACS engine), rather than AMBER's pmemd/sander.
#
# Usage: ./01_amber_to_gromacs.sh <CASE_DIR>
#   e.g. ./01_amber_to_gromacs.sh WT
# Run from 1KS4_GROMACS/ (same convention as submit_replicates.sh), AFTER
# system.prmtop/system.inpcrd already exist in <CASE_DIR> (from tleap).
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR=$1

if [ -z "${CASE_DIR}" ]; then
  echo "Usage: ./01_amber_to_gromacs.sh <CASE_DIR>"
  exit 1
fi

cd "${CASE_DIR}"

PRMTOP=system.prmtop
INPCRD=system.inpcrd

if [ ! -f ${PRMTOP} ] || [ ! -f ${INPCRD} ]; then
  echo "ERROR: ${CASE_DIR}/${PRMTOP}/${INPCRD} not found -- run tleap in ${CASE_DIR} first"
  exit 1
fi

# ParmEd (bundled with AMBER, via amber.python) converts AMBER prmtop/inpcrd
# -> GROMACS gro/top directly, preserving ff14SB parameters (bonded/nonbonded
# terms, including the CYX disulfide).
amber.python << 'PYEOF'
import parmed as pmd
amber = pmd.load_file('system.prmtop', 'system.inpcrd')
amber.save('system.top', overwrite=True)
amber.save('system.gro', overwrite=True)
print("SUCCESS")
PYEOF

if [ ! -s system.top ] || [ ! -s system.gro ]; then
  echo "ERROR: system.top/system.gro missing or empty after ParmEd conversion"
  exit 1
fi

# Unlike pdb2gmx, ParmEd's conversion does NOT generate a posre.itp, so
# nvt.mdp's `define = -DPOSRES` would silently restrain nothing without this.
# Building a heavy-atom restraint file the same way make_posre_backbone.sh
# builds the backbone-only one, selecting by group NAME (not index -- this
# converted topology's group numbering doesn't necessarily match pdb2gmx's).
gmx_mpi make_ndx -f system.gro -o index.ndx << EOF
q
EOF
gmx_mpi genrestr -f system.gro -n index.ndx -o posre.itp -fc 1000 1000 1000 << EOF
Protein
EOF

# Wire posre.itp into system.top automatically (finds the protein
# moleculetype by atom count rather than assuming it's named "Protein").
# Referenced via SCRIPT_DIR so this works regardless of how CASE_DIR is
# nested (this was a real bug in an earlier version of this script).
"${SCRIPT_DIR}/insert_posre_include.sh" system.top POSRES posre.itp

echo "Done. ${CASE_DIR}/system.gro and system.top ready (ff14SB, converted from AMBER)."
echo "Next: ./make_posre_backbone.sh ${CASE_DIR}, then ./submit_replicates.sh ${CASE_DIR} 3"
