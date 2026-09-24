#!/bin/bash
# Build one solvated, ionised, restraint-wired CHARMM36m system.
#
#   ./prep_system.sh q4_holo
#   ./prep_system.sh                 # every system in config.sh SYSTEMS
#
# Same script for apo and holo -- the only difference is whether the glycan gets
# merged in, and that is decided by HAS_LIGAND in systems/<name>.conf, not by a
# separate code path. One protocol, six systems.
#
# Every tool here is SERIAL. Run it as a batch job (./prep_all.batch) or in a
# 1-task interactive allocation:
#     srun -p general-cpu -A <account> -c 4 --mem=8G --pty /bin/bash
#     module load gromacs/2022.5_plumed_gcc_9.5.0_openmpi_4.1.5
#     ./prep_system.sh
#
# A NOTE ON gmx_mpi AND YOUR SHELL: a GROMACS fatal error inside an MPI build
# calls MPI_ABORT, which kills every process in the MPI job -- including an
# interactive shell started with `srun --pty`. That is why pdb2gmx's output goes
# to a FILE here and is printed afterwards, and why prep_all.batch exists. If a
# build fails interactively and you land back on the login node, the log is
# still on disk.
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${SCRIPT_DIR}/config.sh"

# config.sh leaves SYSTEMS empty on purpose: the systems/ directory IS the list,
# so ./new_system.sh is the only place a system gets declared. Setting SYSTEMS in
# config.sh narrows "all" to that subset.
if [ -z "${SYSTEMS:-}" ]; then
  SYSTEMS="$(ls "${SCRIPT_DIR}/systems" 2>/dev/null | sed 's/\.conf$//' | tr '\n' ' ')"
fi

# GROMACS comes from the Apptainer container on this cluster; gmx.sh resolves it
# (and honours a GMX already set in the environment, for other clusters).
PIPELINE_ROOT="${SCRIPT_DIR}"
# shellcheck source=/dev/null
source "${SCRIPT_DIR}/gmx.sh"

export GMX_MAXBACKUP=-1      # overwrite cleanly instead of piling up #file.N#

# ---------------------------------------------------------------- dispatch
if [ "$#" -gt 1 ]; then
  echo "Usage: ./prep_system.sh [system]"; exit 1
fi
if [ "$#" -eq 0 ]; then
  rc=0
  for s in ${SYSTEMS}; do
    echo ""
    echo "##########################################################################"
    echo "### ${s}"
    echo "##########################################################################"
    "${BASH_SOURCE[0]}" "${s}" || { rc=1; echo "*** ${s} FAILED (continuing)"; }
  done
  echo ""
  if [ "${rc}" -ne 0 ]; then
    echo "One or more systems failed. Nothing was submitted."; exit 1
  fi
  echo "All systems built. Next:  ./submit.sh"
  exit 0
fi

SYS=$1
CONF="${SCRIPT_DIR}/systems/${SYS}.conf"
[ -f "${CONF}" ] || {
  echo "ERROR: no systems/${SYS}.conf"
  echo "Known systems: ${SYSTEMS}"; exit 1; }
# shellcheck source=/dev/null
source "${CONF}"

BUILD="${SCRIPT_DIR}/build/${SYS}"
[ -f "${BUILD}/protein.pdb" ] || {
  echo "ERROR: ${BUILD}/protein.pdb not found."
  echo "       Run:  python3 00_prep_inputs.py ${SYS}"; exit 1; }

gmx_check || exit 1

if [ -n "${SLURM_NTASKS}" ] && [ "${SLURM_NTASKS}" -gt 1 ]; then
  echo "*** WARNING: SLURM_NTASKS=${SLURM_NTASKS}, but every tool here is serial."
  echo "***          Build in a 1-task allocation (see the header)."
fi

# ---------------------------------------------------------------- force field
FF_DIR=$(cd "${SCRIPT_DIR}" && ls -d charmm36*.ff 2>/dev/null | head -1 || true)
if [ -z "${FF_DIR}" ] || [ ! -d "${SCRIPT_DIR}/${FF_DIR}" ]; then
  echo "ERROR: no charmm36*.ff in ${SCRIPT_DIR}"
  echo "       Run:  ./01_get_charmm36m.sh charmm36-jul2022.ff.tgz"; exit 1
fi
export GMXLIB="${SCRIPT_DIR}${GMXDATA:+:${GMXDATA}/top}"
FF_NAME="${FF_DIR%.ff}"

# Read one setting out of an mdp. Comments are stripped FIRST: a naive filter
# applied to `ref_t = 300 300 ; config.sh TEMP_K` keeps the dot out of
# "config.sh" and yields "300 300 .", which then fails to compare against 300.
mdp_get () {   # mdp_get <file> <key>
  awk -F';' '{print $1}' "$1" \
    | awk -v k="$2" -F'=' '
        { key=$1; gsub(/^[ \t]+|[ \t]+$/, "", key)
          if (tolower(key) == tolower(k)) {
            v=$2; gsub(/^[ \t]+|[ \t]+$/, "", v); gsub(/[ \t]+/, " ", v)
            print v; exit } }'
}

echo "Force field : ${FF_DIR}   (pdb2gmx -ff ${FF_NAME})"
echo "Water model : ${WATER_MODEL}  (CHARMM-modified TIP3P)"
echo "Ligand      : $([ "${HAS_LIGAND}" = "1" ] && echo "xylohexaose (merged from build/${SYS}/glycan.itp)" || echo "none (apo)")"

# ------------------------------------------------ GUARD: mdp agrees with config
MDP_BAD=0
for m in nvt npt prod; do
  f="${SCRIPT_DIR}/mdp/${m}.mdp"
  for v in $(mdp_get "${f}" ref_t); do
    [ "${v%%.*}" = "${TEMP_K%%.*}" ] || {
      echo "ERROR: mdp/${m}.mdp has ref_t = ${v} but config.sh TEMP_K = ${TEMP_K}"; MDP_BAD=1; }
  done
done
g=$(mdp_get "${SCRIPT_DIR}/mdp/nvt.mdp" gen_temp)
if [ -n "${g}" ] && [ "${g%%.*}" != "${TEMP_K%%.*}" ]; then
  echo "ERROR: mdp/nvt.mdp has gen_temp = ${g} but config.sh TEMP_K = ${TEMP_K}"; MDP_BAD=1
fi
[ "${MDP_BAD}" -eq 0 ] || { echo "Fix mdp/ or config.sh so they agree."; exit 1; }
echo "Guard ok    : ref_t / gen_temp all match TEMP_K = ${TEMP_K} K"

# ------------------------------------------------ GUARD: ion names in this port
for ion in "${SALT_CATION}" "${SALT_ANION}"; do
  grep -qE "^[[:space:]]*${ion}[[:space:]]" "${SCRIPT_DIR}/${FF_DIR}/ions.itp" || {
    echo "ERROR: '${ion}' is not a moleculetype in ${FF_DIR}/ions.itp."
    echo "       CHARMM36 names ions SOD/CLA, not Amber's NA/CL. Available:"
    awk '/^\[ *moleculetype *\]/{want=1;next} want && !/^;/ && NF {print "         " $1; want=0}' \
        "${SCRIPT_DIR}/${FF_DIR}/ions.itp" | sort -u
    exit 1; }
done
echo "Guard ok    : ion names ${SALT_CATION}/${SALT_ANION} exist in ${FF_DIR}"

# ------------------------------------------------ GUARD: the glycan, for holo
if [ "${HAS_LIGAND}" = "1" ]; then
  for f in glycan.pdb glycan.itp posre_lig.itp; do
    [ -s "${BUILD}/${f}" ] || {
      echo ""
      echo "ERROR: ${BUILD}/${f} not found."
      echo "       The CHARMM36 GROMACS port has beta-D-xylose as a MONOSACCHARIDE"
      echo "       only -- no 1->4 glycosidic patch anywhere in it -- so pdb2gmx"
      echo "       cannot build xylohexaose. The topology comes from CHARMM-GUI."
      echo "       See README step 2, then:  python3 02_build_glycan.py ${SYS}"
      exit 1; }
  done
  echo "Guard ok    : glycan topology and restraints present"
fi

cd "${BUILD}"
rm -f processed.gro complex.gro boxed.gro solvated.gro system.gro \
      system.top system.top.bak posre.itp posre_backbone.itp \
      pdb2gmx.log ions.tpr test_nvt.tpr mdout.mdp 2>/dev/null || true

# -------------------------------------------------------------------- pdb2gmx
echo ""
echo "=== pdb2gmx ==="
# -ignh          drop AlphaFold3's (absent) hydrogens and let CHARMM36m
#                templates build a consistent set
# -ter + indices the terminus menu order is port-specific; see the long note in
#                config.sh. Asserted by NAME below, so a wrong guess fails loudly
#                instead of silently building a different molecule.
# no -ss         disulfides are auto-detected from specbond.dat, which is what
#                we want for q4; the result is counted below rather than trusted
#
# Output goes to a FILE, not through tee: a fatal error here calls MPI_ABORT,
# which would take an interactive shell down with it and could lose a buffered
# pipe. The log is printed straight afterwards either way.
set +e
printf '%s\n%s\n' "${NTER_INDEX}" "${CTER_INDEX}" \
  | ${GMX} pdb2gmx -f protein.pdb -o processed.gro -p system.top -i posre.itp \
      -ff "${FF_NAME}" -water "${WATER_MODEL}" -ignh -ter > pdb2gmx.log 2>&1
pdb2gmx_rc=$?
set -e
if [ "${pdb2gmx_rc}" -ne 0 ] || [ ! -s system.top ]; then
  echo "--- pdb2gmx.log ---------------------------------------------------"
  cat pdb2gmx.log
  echo "-------------------------------------------------------------------"
  echo ""
  echo "FATAL: pdb2gmx exited ${pdb2gmx_rc}."
  echo ""
  echo "If the error mentions a terminus or a building block, the menu indices"
  echo "in config.sh (NTER_INDEX=${NTER_INDEX}, CTER_INDEX=${CTER_INDEX}) are"
  echo "wrong for this port. The menus it offered are above -- read the numbers"
  echo "next to NH3+ and COO- off them and put those in config.sh:"
  awk '/Select (start|end) terminus/ {show=1; print "    " $0; next}
       show && /^[[:space:]]*[0-9]+:/ {print "    " $0; next}
       show {show=0}' pdb2gmx.log 2>/dev/null
  exit 1
fi

got_n=$(awk '/Start terminus/{print $NF; exit}' pdb2gmx.log)
got_c=$(awk '/End terminus/{print $NF; exit}'   pdb2gmx.log)
echo "Termini     : requested ${NTER_INDEX}/${CTER_INDEX}, applied ${got_n:-?} / ${got_c:-?}"
if [ "${got_n}" != "NH3+" ] || [ "${got_c}" != "COO-" ]; then
  echo "*** FATAL: expected NH3+ / COO-. A wrong terminus is not cosmetic -- it"
  echo "*** changes the net charge. Menus this port offered:"
  awk '/Select (start|end) terminus/ {show=1; print "    " $0; next}
       show && /^[[:space:]]*[0-9]+:/ {print "    " $0; next}
       show {show=0}' pdb2gmx.log
  echo "*** Put the right numbers in config.sh NTER_INDEX / CTER_INDEX."
  exit 1
fi

# ---------------------------------------------------------------- topology checks
echo ""
echo "=== topology checks ==="

# THE BOND ITSELF is the test, not a residue name.
#
# pdb2gmx sets the CYS2 *rtp block* for a disulfide-bonded cysteine but leaves
# the residue NAME in [ atoms ] as CYS. An earlier version of this check counted
# CYS2 residues, found zero, and declared a perfectly good q4 build broken --
# while pdb2gmx's own log said
#     Linking CYS-97 SG-737 and CYS-141 SG-1078...
# So ask the topology whether two SG atoms are bonded TO EACH OTHER. That is the
# thing that actually matters and no naming convention can fool it.
ndisulf=$(awk '
  /^\[ *atoms *\]/  {sec="atoms"; next}
  /^\[ *bonds *\]/  {sec="bonds"; next}
  /^\[/              {sec="";      next}
  sec=="atoms" && $1 ~ /^[0-9]+$/ && $5=="SG"                        {sg[$1]=1}
  sec=="bonds" && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ &&
                 ($1 in sg) && ($2 in sg)                            {n++}
  END {print n+0}' system.top)

# Cross-check against what pdb2gmx said it did. The two must agree: a "Linking"
# line with no bond in the topology (or the reverse) means something is wrong
# that neither test alone would catch.
nlink=$(grep -c '^Linking CYS' pdb2gmx.log 2>/dev/null || true)
nlink=${nlink:-0}
if [ "${ndisulf}" -ne "${nlink}" ]; then
  echo "   *** INCONSISTENT: ${ndisulf} SG-SG bond(s) in system.top but pdb2gmx"
  echo "   *** reported ${nlink} 'Linking CYS' line(s). Stop and look at pdb2gmx.log."
  exit 1
fi
echo "Disulfide bridges : ${ndisulf}   (${SYS}.conf expects ${EXPECT_DISULFIDES})"
if [ "${ndisulf}" -ne "${EXPECT_DISULFIDES}" ]; then
  echo "   *** MISMATCH. ${ndisulf} SG-SG bond(s) found in the topology;"
  echo "   *** pdb2gmx.log reports ${nlink} 'Linking CYS' line(s)."
  if [ "${EXPECT_DISULFIDES}" -gt 0 ]; then
    echo "   *** The engineered Cys97-Cys141 bond was NOT made. Those two would"
    echo "   *** run as free thiols -- a different molecule from the one designed,"
    echo "   *** and the entire point of the q4 variant would be missing."
    echo "   *** Look for these lines in pdb2gmx.log -- they say what it did:"
    echo "   ***   grep -A6 'Special Atom Distance matrix' ${BUILD}/pdb2gmx.log"
    echo "   ***   grep 'Linking CYS' ${BUILD}/pdb2gmx.log"
    echo "   *** 00_prep_inputs.py measured the SG-SG distance at ~0.2 nm, which is"
    echo "   *** exactly specbond.dat's ideal, so the geometry is not the problem."
  fi
  exit 1
fi

# Accept every spelling this port can produce. aminoacids.rtp defines
# HSD/HSE/HSP/HSPM while aminoacids.r2b maps HISD/HISE/HISH/HIS1 onto them, and
# pdb2gmx writes the ALIAS into the topology -- a check that accepted only
# HSD/HSE/HSP once scored a perfectly good build as zero histidines.
HIS_RE='^(HIS|HISD|HISE|HISH|HIS1|HSD|HSE|HSP|HSPM|HID|HIE|HIP)$'
echo "Histidine states:"
awk '/^\[ *atoms *\]/,/^\[ *bonds *\]/' system.top \
  | awk -v re="${HIS_RE}" '$4 ~ re {print $4, $3}' | sort -u \
  | awk '{print $1}' | sort | uniq -c | awk '{printf "   %6d  %s\n", $1, $2}'
nhis=$(awk '/^\[ *atoms *\]/,/^\[ *bonds *\]/' system.top \
       | awk -v re="${HIS_RE}" '$4 ~ re {print $3}' | sort -u | wc -l | tr -d ' ')
echo "Total histidines  : ${nhis}   (config.sh expects ${EXPECT_HIS})"
if [ "${nhis}" -ne "${EXPECT_HIS}" ]; then
  echo "   *** MISMATCH. Residue names actually in the topology:"
  awk '/^\[ *atoms *\]/,/^\[ *bonds *\]/' system.top \
    | awk 'NF>6 && $1 ~ /^[0-9]+$/ {print $4}' | sort -u | tr '\n' ' ' \
    | fold -s -w 66 | sed 's/^/   ***   /'
  echo ""; exit 1
fi
if awk '/^\[ *atoms *\]/,/^\[ *bonds *\]/' system.top \
     | awk '$4=="HSP"||$4=="HSPM"||$4=="HISH"||$4=="HIP"' | grep -q .; then
  echo "   *** A PROTONATED histidine (+1) is present. config.sh builds all three"
  echo "   *** neutral, which is correct at pH 7. Stop and decide deliberately."
  exit 1
fi

CHG=$(grep -iE "total charge" pdb2gmx.log | tail -1 \
      | grep -oE '[-+]?[0-9]+\.[0-9]+' | tail -1)
if [ -n "${CHG}" ]; then
  CHG_INT=$(awk -v c="${CHG}" 'BEGIN{printf "%d", (c<0? c-0.5 : c+0.5)}')
  echo "pdb2gmx net charge: ${CHG}   (config.sh expects ${EXPECT_CHARGE})"
  [ "${CHG_INT}" -eq "${EXPECT_CHARGE}" ] || {
    echo "   *** MISMATCH. Check the histidine tautomers and the termini before"
    echo "   *** going further: all six systems must reach ${EXPECT_CHARGE}."
    exit 1; }
else
  echo "*** WARNING: no total-charge line in pdb2gmx.log. Expected ${EXPECT_CHARGE}."
fi

# ----------------------------------------------------------------- restraints
echo ""
echo "=== restraint files ==="
# Generated against processed.gro (PROTEIN ONLY). genrestr writes indices
# relative to the structure it is given, and this .itp is #included inside the
# PROTEIN moleculetype, where indices must be protein-local. Running it on the
# merged or solvated file happens to work only because the protein comes first,
# which is a coincidence of ordering, not a guarantee.
echo "Backbone" | ${GMX} genrestr -f processed.gro -o posre_backbone.itp \
     -fc 1000 1000 1000 > genrestr.log 2>&1
[ -s posre_backbone.itp ] || {
  echo "FATAL: genrestr produced no posre_backbone.itp"; cat genrestr.log; exit 1; }

NPROT=$(awk '
  /^\[ *atoms *\]/      {inatoms=1; n=0; next}
  /^\[/                 {if (inatoms && n>max) max=n; inatoms=0; next}
  inatoms && $1 ~ /^[0-9]+$/ {n++}
  END                   {if (inatoms && n>max) max=n; print max+0}' system.top)
echo "protein moleculetype: ${NPROT} atoms"
for f in posre.itp posre_backbone.itp; do
  hi=$(grep -E '^[[:space:]]*[0-9]+' "${f}" | awk '{if($1>m) m=$1} END{print m+0}')
  lo=$(grep -E '^[[:space:]]*[0-9]+' "${f}" | awk 'NR==1{m=$1} $1<m{m=$1} END{print m+0}')
  echo "  ${f}: indices ${lo}..${hi}"
  if [ "${hi}" -gt "${NPROT}" ] || [ "${lo}" -lt 1 ]; then
    echo "*** FATAL: ${f} indexes atom ${hi}, outside the protein moleculetype"
    echo "*** (1..${NPROT}). The restraints would land on the wrong atoms."; exit 1
  fi
done
nheavy=$(grep -cE '^[[:space:]]*[0-9]+' posre.itp || true)
nbb=$(grep -cE '^[[:space:]]*[0-9]+' posre_backbone.itp || true)
echo "posre.itp          (-DPOSRES,    protein heavy atoms) : ${nheavy}"
echo "posre_backbone.itp (-DPOSRES_BB, protein backbone)    : ${nbb}"
echo "   (${N_RESIDUES} residues -> expect ~$((N_RESIDUES * 3)) backbone atoms)"
[ "${nbb}" -lt "${nheavy}" ] || echo "*** WARNING: backbone count >= heavy count -- check genrestr's selection."

if grep -q "#ifdef POSRES_BB" system.top; then
  echo "system.top already wires POSRES_BB -- not inserting again."
else
  "${SCRIPT_DIR}/insert_posre_include.sh" system.top POSRES_BB posre_backbone.itp \
    | sed 's/^/   /'
fi

# ------------------------------------------------------------ merge the glycan
START_GRO=processed.gro
if [ "${HAS_LIGAND}" = "1" ]; then
  echo ""
  echo "=== merging the xylohexaose ==="
  python3 "${SCRIPT_DIR}/merge_holo.py" processed.gro glycan.pdb glycan.itp \
          system.top complex.gro
  START_GRO=complex.gro
fi

# ------------------------------------------------------------- box + solvate
echo ""
echo "=== box (${BOX_TYPE}, ${BOX_PAD} nm padding) + solvation ==="
${GMX} editconf -f "${START_GRO}" -o boxed.gro -c -d "${BOX_PAD}" -bt "${BOX_TYPE}" \
    > editconf.log 2>&1
${GMX} solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p system.top \
    > solvate.log 2>&1
[ -s solvated.gro ] || { echo "FATAL: solvate produced nothing"; cat solvate.log; exit 1; }

# ------------------------------------------------------------------------ ions
echo ""
echo "=== ions (neutralise ${EXPECT_CHARGE} + ${SALT_MM} mM ${SALT_CATION}${SALT_ANION}) ==="
# The bulk-salt count is computed from the WATER COUNT, not from genion's -conc.
# -conc sizes salt from the BOX VOLUME (C x V x N_A), which counts the volume the
# protein occupies as if it were solvent; the water-count recipe does not. The
# two disagree by ~10% here, and with six systems of slightly different size that
# would put a systematic ionic-strength difference between them.
NSOL=$(awk '/^\[ *molecules *\]/,0' system.top | awk '$1=="SOL"{n=$2} END{print n+0}')
[ "${NSOL:-0}" -gt 0 ] || { echo "FATAL: no SOL count in [ molecules ]"; exit 1; }
NSALT=$(awk -v w="${NSOL}" -v c="${SALT_MM}" 'BEGIN{printf "%d", (c/1000.0)*w/55.5 + 0.5}')
echo "waters      : ${NSOL}"
echo "bulk salt   : ${NSALT} ${SALT_CATION} / ${NSALT} ${SALT_ANION} pairs for ${SALT_MM} mM"

${GMX} grompp -f "${SCRIPT_DIR}/mdp/min.mdp" -c solvated.gro -p system.top \
    -o ions.tpr -maxwarn 2 > grompp_ions.log 2>&1 || {
  echo "FATAL: grompp before genion failed:"; tail -40 grompp_ions.log; exit 1; }
echo "SOL" | ${GMX} genion -s ions.tpr -o system.gro -p system.top \
    -pname "${SALT_CATION}" -nname "${SALT_ANION}" \
    -np "${NSALT}" -nn "${NSALT}" -neutral > genion.log 2>&1
[ -s system.gro ] || { echo "FATAL: genion produced nothing"; cat genion.log; exit 1; }

WANT_CAT=$(awk -v n="${NSALT}" -v q="${EXPECT_CHARGE}" 'BEGIN{print n + (q<0 ? -q : 0)}')
WANT_ANI=$(awk -v n="${NSALT}" -v q="${EXPECT_CHARGE}" 'BEGIN{print n + (q>0 ?  q : 0)}')
GOT_CAT=$(awk '/^\[ *molecules *\]/,0' system.top | awk -v a="${SALT_CATION}" '$1==a{n+=$2} END{print n+0}')
GOT_ANI=$(awk '/^\[ *molecules *\]/,0' system.top | awk -v b="${SALT_ANION}"  '$1==b{n+=$2} END{print n+0}')
printf '   %-4s %s  (expected %s = %s bulk + neutralising)\n' "${SALT_CATION}" "${GOT_CAT}" "${WANT_CAT}" "${NSALT}"
printf '   %-4s %s  (expected %s = %s bulk + neutralising)\n' "${SALT_ANION}" "${GOT_ANI}" "${WANT_ANI}" "${NSALT}"
if [ "${GOT_CAT}" -ne "${WANT_CAT}" ] || [ "${GOT_ANI}" -ne "${WANT_ANI}" ]; then
  echo "*** MISMATCH between ions requested and ions placed."; exit 1
fi

# --------------------------------------------------- the check that matters most
# A real grompp of the FIRST MD STAGE, against the finished system. This is the
# only test that exercises everything at once: that -DPOSRES resolves, that the
# restraint indices are in range for BOTH molecules, that `Protein` and
# `Non-Protein` exist as coupling groups now that a ligand is in the box, and
# that every #include resolves. Finding any of that out here costs seconds;
# finding it out after the job queues costs hours.
echo ""
echo "=== dry-run grompp of nvt.mdp (the real test) ==="
if ${GMX} grompp -f "${SCRIPT_DIR}/mdp/nvt.mdp" -c system.gro -r system.gro \
      -p system.top -o test_nvt.tpr -maxwarn 1 > grompp_test.log 2>&1; then
  nrestr=$(grep -ciE 'position restraint' grompp_test.log || true)
  echo "  ok  nvt.mdp grompps cleanly against the finished system"
  awk '/Number of degrees of freedom/{print "  " $0}' grompp_test.log | head -3
  rm -f test_nvt.tpr
else
  echo "*** FATAL: grompp of nvt.mdp failed. The system is NOT ready."
  echo "--- grompp_test.log (tail) ----------------------------------------"
  tail -40 grompp_test.log
  echo "-------------------------------------------------------------------"
  exit 1
fi

for d in POSRES POSRES_BB; do
  grep -q "#ifdef ${d}" system.top || { echo "FATAL: no #ifdef ${d} in system.top"; exit 1; }
done

# ---------------------------------------------------------------- final report
echo ""
echo "=== ${SYS}: final system ==="
natoms=$(sed -n '2p' system.gro | tr -d ' ')
printf "   total atoms      %s\n" "${natoms}"
printf "   box vectors      %s\n" "$(tail -1 system.gro)"
echo "[ molecules ]:"
awk '/^\[ *molecules *\]/,0' system.top | grep -v '^\[' | grep -v '^;' \
  | awk 'NF{printf "   %s\n", $0}'
echo ""
echo "Built. Next:  ./submit.sh ${SYS}"
