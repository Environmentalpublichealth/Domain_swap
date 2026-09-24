#!/bin/bash
#SBATCH --partition=general-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
# No --time: the analysis job takes the partition default, same as the MD jobs.
#SBATCH --output=analysis_out_%j.log
#SBATCH --error=analysis_err_%j.log
# ============================================================================
# Backbone RMSD, per-residue backbone RMSF and per-residue AMIDE HYDROGEN RMSF
# for ONE replicate of the CHARMM36m XYN runs.
#
#   cd runs/1xyn_apo/rep1 && sbatch ../../../run_analysis.sh
#
# SYSTEM and REP come from the directory path, not from arguments, so the
# output tag can never disagree with the directory the numbers came from.
# Submit every replicate at once with ../../submit_analysis.sh instead.
#
# Options (environment variables):
#   B_NS=50      discard the first 50 ns from RMSD/RMSF (default 0: keep all)
#   EXTRAS=1     also compute radius of gyration and DSSP secondary structure
#
# Outputs, all tagged <system>_<rep>_c36m:
#   rmsd_<tag>.xvg          backbone RMSD vs time, ns / nm
#   rmsf_<tag>.xvg          per-residue backbone RMSF, nm
#   rmsf_amideH_<tag>.xvg   per-residue amide-H RMSF, nm  <- HDX comparison
#   amideH_atoms_<tag>.pdb  one record per selected amide H, to audit by eye
#   frame0_<tag>.pdb        first frame, to map MD numbering onto HDX peptides
#
# Nothing here averages replicates. Each replicate is analysed on its own and
# reported on its own; a mean across them would hide where they disagree.
# ============================================================================
set -u

# ---------------------------------------------------------------- where am I
REP=$(basename "$PWD")
SYSTEM=$(basename "$(dirname "$PWD")")
PIPELINE_ROOT=$(cd ../../.. 2>/dev/null && pwd)

case "${REP}" in rep[0-9]*) ;; *)
  echo "FATAL: '${REP}' is not a rep directory."
  echo "       Run this from runs/<system>/rep<N>/, e.g."
  echo "         cd runs/1xyn_apo/rep1 && sbatch ../../../run_analysis.sh"
  exit 1 ;;
esac
if [ ! -f "${PIPELINE_ROOT}/systems/${SYSTEM}.conf" ]; then
  echo "FATAL: no systems/${SYSTEM}.conf under ${PIPELINE_ROOT}"
  echo "       Expected layout: <pipeline>/runs/<system>/rep<N>/"
  exit 1
fi
TAG="${SYSTEM}_${REP}_c36m"

# shellcheck source=/dev/null
[ -r "${PIPELINE_ROOT}/config.sh" ] && source "${PIPELINE_ROOT}/config.sh"
# the per-system conf too: N_RESIDUES is what the CA count is checked against,
# and HAS_LIGAND says whether "Protein" leaves a ligand out of the fit group
# shellcheck source=/dev/null
source "${PIPELINE_ROOT}/systems/${SYSTEM}.conf"
# GROMACS from the Apptainer container (see gmx.sh); no module on this cluster.
# shellcheck source=/dev/null
source "${PIPELINE_ROOT}/gmx.sh"

gmx_check || exit 1

# Discard equilibration from RMSD/RMSF if asked. Empty by default: the whole
# production run is analysed, and any truncation is a visible, deliberate choice.
B_NS=${B_NS:-0}
TRJ_OPTS=""
[ "${B_NS}" != "0" ] && TRJ_OPTS="-b $(awk -v n="${B_NS}" 'BEGIN{print n*1000}')"

echo "=========================================================="
echo "system     : ${SYSTEM}"
echo "replicate  : ${REP}"
echo "tag        : ${TAG}"
echo "pipeline   : ${PIPELINE_ROOT}"
echo "analysing  : $( [ "${B_NS}" = "0" ] && echo "the whole trajectory" \
                      || echo "from ${B_NS} ns onward (B_NS=${B_NS})" )"
echo "=========================================================="

[ -s prod_1.tpr ] || { echo "FATAL: no prod_1.tpr here."; exit 1; }
[ -s prod_1.xtc ] || { echo "FATAL: no prod_1.xtc here."; exit 1; }
if [ ! -s prod_1.gro ]; then
  echo "*** WARNING: prod_1.gro absent -- production has not reached ${PROD_NS:-100} ns."
  echo "*** Analysing the partial trajectory anyway; note the length in the RMSD plot."
fi

# ============================================================ periodic images
# A. put the protein back in one piece and centre it
echo ""; echo "=== A. centering ==="
echo "Protein System" | ${GMX} trjconv -s prod_1.tpr -f prod_1.xtc \
  -o prod_center.xtc -pbc mol -center \
  || { echo "FATAL: trjconv -pbc mol failed."; exit 1; }

# B. remove overall rotation/translation, fitting on Backbone. Every quantity
#    below is computed from this ONE fitted trajectory, so RMSD, backbone RMSF
#    and amide-H RMSF all share the same reference frame.
echo ""; echo "=== B. rot+trans fit ==="
echo "Backbone System" | ${GMX} trjconv -s prod_1.tpr -f prod_center.xtc \
  -o prod_final.xtc -fit rot+trans \
  || { echo "FATAL: trjconv -fit rot+trans failed."; exit 1; }

# ==================================================================== 1. RMSD
# The tool is `gmx rms`, not `gmx rmsd`.
echo ""; echo "=== 1. backbone RMSD ==="
echo "Backbone Backbone" | ${GMX} rms -s prod_1.tpr -f prod_final.xtc ${TRJ_OPTS} \
  -o rmsd_${TAG}.xvg -tu ns \
  || { echo "FATAL: gmx rms failed."; exit 1; }
echo "wrote rmsd_${TAG}.xvg"

# ==================================================================== 2. RMSF
echo ""; echo "=== 2. per-residue backbone RMSF ==="
echo "Backbone" | ${GMX} rmsf -s prod_1.tpr -f prod_final.xtc ${TRJ_OPTS} \
  -o rmsf_${TAG}.xvg -res \
  || { echo "FATAL: gmx rmsf failed."; exit 1; }
echo "wrote rmsf_${TAG}.xvg"

# ========================================================== 3. amide-H RMSF
# The backbone amide hydrogen is in none of GROMACS's stock groups, so it needs
# its own index.
#
# Naming: CHARMM36 calls it HN, Amber/ff14SB calls it H. Matching both means
# this script also runs unchanged on the ff14SB trajectories.
#
# Two exclusions, both physically required:
#   * PROLINE has no amide hydrogen (its N is tertiary). It cannot exchange and
#     must not appear in the output.
#   * The N-TERMINUS is an -NH3+ group (HT1-3 in CHARMM, H1-3 in Amber), not a
#     peptide amide. Matching the exact names H/HN already excludes it.
echo ""; echo "=== 3. amide-hydrogen RMSF ==="
${GMX} select -s prod_1.tpr -on amideH_${TAG}.ndx \
  -select 'group "Protein" and name H HN and not resname PRO' \
  || { echo "FATAL: amide-H selection failed."; exit 1; }

# Check the selection rather than trusting it: a silently wrong selection makes
# a plausible-looking plot of the wrong atoms, the worst failure mode there is.
${GMX} select -s prod_1.tpr -on _ca_${TAG}.ndx \
  -select 'group "Protein" and name CA' >/dev/null 2>&1
${GMX} select -s prod_1.tpr -on _pro_${TAG}.ndx \
  -select 'group "Protein" and name CA and resname PRO' >/dev/null 2>&1
count_ndx () { grep -v '^\[' "$1" 2>/dev/null | wc -w | tr -d ' '; }
nH=$(count_ndx amideH_${TAG}.ndx)
nCA=$(count_ndx _ca_${TAG}.ndx)
nPRO=$(count_ndx _pro_${TAG}.ndx)
expect=$((nCA - nPRO - 1))          # -1 for the N-terminal NH3+
echo "residues (CA)          : ${nCA}   (systems/${SYSTEM}.conf says ${N_RESIDUES:-?})"
if [ -n "${N_RESIDUES:-}" ] && [ "${nCA}" -ne "${N_RESIDUES}" ]; then
  echo "*** WARNING: the trajectory has ${nCA} residues, the conf expects ${N_RESIDUES}."
  echo "***          Check you are analysing the system you think you are."
fi
echo "prolines               : ${nPRO}"
echo "amide H atoms selected : ${nH}   (expected ${expect} = residues - prolines - N-term)"
if [ "${nH}" -eq 0 ]; then
  echo "*** FATAL: no amide hydrogens selected. Check the atom naming with:"
  echo "***   gmx_mpi dump -s prod_1.tpr | grep -m20 'name=\"H'"
  exit 1
fi
if [ "${nH}" -ne "${expect}" ]; then
  echo "*** WARNING: amide-H count is off by $((nH - expect)). Benign causes: extra"
  echo "***          chains (one N-term each) or a capped terminus. Audit the pdb below."
fi
rm -f _ca_${TAG}.ndx _pro_${TAG}.ndx

# -res averages per residue. With one amide H per residue that is the atom's own
# RMSF, but it puts RESIDUE NUMBERS on the x-axis, which is what mapping onto HDX
# peptides needs. Prolines are simply absent, and that gap is meaningful.
#
# No -fit: prod_final.xtc is already rot+trans fitted, exactly like the backbone
# RMSF above. Fitting twice would make the two incomparable.
echo 0 | ${GMX} rmsf -s prod_1.tpr -f prod_final.xtc ${TRJ_OPTS} \
  -n amideH_${TAG}.ndx -o rmsf_amideH_${TAG}.xvg -res \
  || { echo "FATAL: amide-H RMSF failed."; exit 1; }
echo "wrote rmsf_amideH_${TAG}.xvg"

# An explicit atom -> residue map, so the selection can be checked by eye.
echo 0 | ${GMX} trjconv -s prod_1.tpr -f prod_final.xtc \
  -n amideH_${TAG}.ndx -o amideH_atoms_${TAG}.pdb -dump 0 >/dev/null 2>&1 \
  && echo "wrote amideH_atoms_${TAG}.pdb  (one record per selected H)"

# ================================================== 4. numbering map for HDX
echo ""; echo "=== 4. first frame (MD numbering -> HDX peptides) ==="
echo "Protein" | ${GMX} trjconv -s prod_1.tpr -f prod_final.xtc \
  -o frame0_${TAG}.pdb -dump 0 >/dev/null 2>&1 \
  && echo "wrote frame0_${TAG}.pdb"

# ======================================================== 5. optional extras
if [ "${EXTRAS:-0}" = "1" ]; then
  echo ""; echo "=== 5. radius of gyration + DSSP ==="
  echo "Protein" | ${GMX} gyrate -s prod_1.tpr -f prod_final.xtc \
    -o rg_${TAG}.xvg >/dev/null 2>&1 && echo "wrote rg_${TAG}.xvg"
  echo "Protein" | ${GMX} do_dssp -s prod_1.tpr -f prod_final.xtc \
    -o ss_${TAG}.xpm -sc scount_${TAG}.xvg >/dev/null 2>&1 \
    && echo "wrote ss_${TAG}.xpm + scount_${TAG}.xvg" \
    || echo "*** do_dssp failed -- load a dssp module or set \$DSSP, then rerun."
fi

# prod_final.xtc is KEPT: it is the fitted trajectory every number above came
# from, and HDXer reads it directly. Only the intermediate goes.
rm -f prod_center.xtc

echo ""
echo "=========================================================="
echo "DONE  ${TAG}"
ls -1 rmsd_${TAG}.xvg rmsf_${TAG}.xvg rmsf_amideH_${TAG}.xvg 2>/dev/null | sed 's/^/  /'
echo "=========================================================="
