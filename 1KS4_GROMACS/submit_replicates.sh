#!/bin/bash
# Sets up and submits N independent replicates for one case.
#
# Usage: ./submit_replicates.sh <CASE_DIR> [N_REPS]
#   e.g. ./submit_replicates.sh WT 3
#
# Expects <CASE_DIR>/system.gro and system.top to already exist (run
# 00_prepare_system.sh and make_posre_backbone.sh once inside CASE_DIR first).
# Creates CASE_DIR/rep1 .. rep<N_REPS>, copies the built system + restraint
# .itp files into each, and submits run_1KS4.batch from each rep directory so
# every replicate gets its own independent min->nvt->equi->prod run (nvt.mdp's
# gen_seed=-1 gives each one different starting velocities).
set -e

CASE_DIR=$1
NREPS=${2:-3}

if [ -z "${CASE_DIR}" ]; then
  echo "Usage: ./submit_replicates.sh <CASE_DIR> [N_REPS]"
  exit 1
fi

for f in system.gro system.top; do
  if [ ! -f "${CASE_DIR}/${f}" ]; then
    echo "ERROR: ${CASE_DIR}/${f} not found -- build the system in ${CASE_DIR} first"
    exit 1
  fi
done

for f in posre.itp posre_backbone.itp; do
  if [ ! -f "${CASE_DIR}/${f}" ]; then
    echo "ERROR: ${CASE_DIR}/${f} not found. Without it, nvt.mdp/equi.mdp's"
    echo "-DPOSRES/-DPOSRES_BB restraints silently apply to nothing (no error,"
    echo "just unrestrained equilibration). Run make_posre_backbone.sh (and, if"
    echo "using the AMBER route, confirm 01_amber_to_gromacs.sh's posre.itp step"
    echo "ran) in ${CASE_DIR} first."
    exit 1
  fi
done
echo "Also confirm system.top actually #includes both posre.itp and"
echo "posre_backbone.itp under the protein moleculetype -- this script can't"
echo "verify that part for you (see make_posre_backbone.sh's printed note)."

for i in $(seq 1 ${NREPS}); do
  REPDIR=${CASE_DIR}/rep${i}
  mkdir -p ${REPDIR}
  cp ${CASE_DIR}/system.gro ${CASE_DIR}/system.top ${REPDIR}/
  cp ${CASE_DIR}"/"*.itp ${REPDIR}/ 2>/dev/null || true
  echo "Submitting ${REPDIR}..."
  (cd ${REPDIR} && sbatch ../../run_1KS4.batch)
done
