#!/bin/bash
# Generates a backbone-only position-restraint file for the NPT equilibration
# step (equi.mdp uses -DPOSRES_BB, separate from the heavy-atom posre.itp
# used by -DPOSRES in nvt.mdp).
#
# Usage: ./make_posre_backbone.sh <CASE_DIR>
#   e.g. ./make_posre_backbone.sh WT
# Run from 1KS4_GROMACS/ (same convention as submit_replicates.sh), AFTER
# system.gro/system.top already exist in <CASE_DIR>.
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR=$1

if [ -z "${CASE_DIR}" ]; then
  echo "Usage: ./make_posre_backbone.sh <CASE_DIR>"
  exit 1
fi

cd "${CASE_DIR}"

if [ ! -f system.gro ] || [ ! -f system.top ]; then
  echo "ERROR: ${CASE_DIR}/system.gro or system.top not found -- build the system first"
  exit 1
fi

# index.ndx may already exist (e.g. from 01_amber_to_gromacs.sh) -- reuse it
# instead of regenerating, since it has no reason to differ.
if [ ! -f index.ndx ]; then
  gmx_mpi make_ndx -f system.gro -o index.ndx << EOF
q
EOF
fi

# Selecting by group NAME, not index -- group numbering isn't guaranteed to
# match between a pdb2gmx-built topology and an acpype-converted one.
gmx_mpi genrestr -f system.gro -n index.ndx -o posre_backbone.itp -fc 1000 1000 1000 << EOF
Backbone
EOF

# Wire posre_backbone.itp into system.top automatically, same mechanism as
# 01_amber_to_gromacs.sh uses for posre.itp -- finds the protein moleculetype
# by atom count rather than assuming a fixed name.
"${SCRIPT_DIR}/insert_posre_include.sh" system.top POSRES_BB posre_backbone.itp

echo "Done. ${CASE_DIR}/posre_backbone.itp written and wired into system.top."
echo "Next: ./submit_replicates.sh ${CASE_DIR} 3"
