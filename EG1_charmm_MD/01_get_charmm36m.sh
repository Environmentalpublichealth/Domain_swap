#!/bin/bash
# Installs the CHARMM36m force field into this pipeline directory.
#
# CHARMM36m is NOT bundled with GROMACS -- unlike the ff99/oplsaa sets that ship
# in $GMXDATA/top. It has to be downloaded from the MacKerell lab as a .ff
# directory that pdb2gmx then finds in the CURRENT WORKING DIRECTORY or in
# $GMXLIB.
#
# This pipeline needs the carbohydrate half of it as well as the protein half:
# carb.rtp supplies the beta-D-xylose parameters the xylohexaose is built from.
#
# Usage:
#   ./01_get_charmm36m.sh                       # try to download
#   ./01_get_charmm36m.sh /path/to/charmm36-jul2022.ff.tgz   # use a local copy
#
# If the download fails (the MacKerell site reorganises its URLs periodically,
# and compute nodes are often firewalled), fetch it by hand from
#   http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs
# choose the newest "CHARMM36m" GROMACS port, scp it to this directory, and
# re-run this script with the tarball path as the argument.
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

LOCAL_TGZ=$1

# Known-good as of writing; override with CHARMM_URL=... if the site moved.
CHARMM_URL=${CHARMM_URL:-"http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz"}

existing=$(ls -d charmm36*.ff 2>/dev/null | head -1 || true)
if [ -n "${existing}" ]; then
  echo "Already present: ${existing}"
  echo "Delete it first if you want to reinstall."
  exit 0
fi

if [ -n "${LOCAL_TGZ}" ]; then
  [ -f "${LOCAL_TGZ}" ] || { echo "ERROR: ${LOCAL_TGZ} not found"; exit 1; }
  echo "Unpacking local ${LOCAL_TGZ} ..."
  tar xzf "${LOCAL_TGZ}"
else
  echo "Downloading CHARMM36m from:"
  echo "  ${CHARMM_URL}"
  if ! curl -fsSL -o charmm36.ff.tgz "${CHARMM_URL}"; then
    echo ""
    echo "DOWNLOAD FAILED. This is expected on a firewalled compute node, and"
    echo "the MacKerell site also moves these URLs periodically."
    echo ""
    echo "Do this instead, from a machine with web access:"
    echo "  1. open http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs"
    echo "  2. download the newest CHARMM36m GROMACS port (charmm36-*.ff.tgz)"
    echo "  3. scp it to ${SCRIPT_DIR}/"
    echo "  4. ./01_get_charmm36m.sh charmm36-<date>.ff.tgz"
    exit 1
  fi
  tar xzf charmm36.ff.tgz
  rm -f charmm36.ff.tgz
fi

FF=$(ls -d charmm36*.ff 2>/dev/null | head -1 || true)
[ -n "${FF}" ] || { echo "ERROR: no charmm36*.ff directory after unpacking"; exit 1; }

# --- Verify it is actually usable, rather than just present.
fail=0
for f in forcefield.itp forcefield.doc aminoacids.rtp ions.itp tip3p.itp; do
  if [ -s "${FF}/${f}" ]; then
    echo "  ok      ${FF}/${f}"
  else
    echo "  MISSING ${FF}/${f}"
    fail=1
  fi
done

# CHARMM36m vs plain CHARMM36. The decisive evidence is the PARAMETER FILE the
# port was built from: par_all36m_prot.prm is the 36m protein set (modified
# backbone CMAP), par_all36_prot.prm is plain C36. Merely finding the string
# "CHARMM36m" is NOT sufficient -- a plain-C36 port still cites the Huang 2016
# CHARMM36m paper in its reference list, so that test passes either way.
if grep -q "par_all36m_prot" "${FF}/forcefield.doc" 2>/dev/null; then
  echo "  ok      built from par_all36m_prot.prm -- this IS CHARMM36m"
elif grep -q "par_all36_prot" "${FF}/forcefield.doc" 2>/dev/null; then
  echo "  *** WARNING: built from par_all36_prot.prm -- this is plain CHARMM36,"
  echo "  ***          NOT CHARMM36m. The ADHA document specifies CHARMM36m."
  echo "  ***          Download a newer port before running production."
else
  echo "  *** WARNING: could not determine C36 vs C36m from forcefield.doc."
  echo "  ***          Check by hand:  grep par_all36 ${FF}/forcefield.doc"
fi

# --- Ion residue names. CHARMM uses SOD/CLA, not Amber's NA/CL. genion fails
# --- with an unhelpful message if you hand it the wrong name, so surface the
# --- real contents of ions.itp here.
echo ""
echo "Ion moleculetypes available in ${FF}/ions.itp:"
awk '/^\[ *moleculetype *\]/{want=1;next} want && !/^;/ && NF {print "  " $1; want=0}' \
    "${FF}/ions.itp" 2>/dev/null | sort -u | paste -sd' ' - | fold -s -w 76

[ "${fail}" -eq 0 ] || { echo ""; echo "FATAL: force field incomplete."; exit 1; }

echo ""
echo "Done. Force field installed as ${SCRIPT_DIR}/${FF}"
echo "pdb2gmx will find it because prep_system.sh runs with -ff ${FF%.ff}"
echo "from a directory where this .ff is visible."
echo ""
echo "Next: sbatch prep_all.batch      (builds all six systems, ~1 h)"
echo "      tail -f prep_*.log"
