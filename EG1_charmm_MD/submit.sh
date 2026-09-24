#!/bin/bash
# Set up and submit the replicates for one system, or for all six.
#
#   ./submit.sh                 # all systems in config.sh SYSTEMS
#   ./submit.sh q4_holo         # one system
#   ./submit.sh q4_holo 5       # one system, 5 replicates instead of N_REPS
#
# Each replicate gets its own copy of the built system and its own full
# min -> nvt -> npt -> prod chain, so nvt.mdp's gen_seed = -1 gives each one
# genuinely independent starting velocities. They are INDEPENDENT EXPERIMENTS
# and get reported individually -- never collapsed into a mean with error bars.
#
# Each chain self-resubmits across wall boundaries until its 100 ns is done.
# You submit once.
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

cd "${SCRIPT_DIR}"

TARGET=$1
NREPS=${2:-${N_REPS:-3}}

if [ -n "${TARGET}" ]; then
  case " ${SYSTEMS} " in
    *" ${TARGET} "*) TARGETS="${TARGET}" ;;
    *) echo "Usage: ./submit.sh [system] [N_REPS]"
       echo "  '${TARGET}' is not one of: ${SYSTEMS}"
       exit 1 ;;
  esac
else
  TARGETS="${SYSTEMS}"
fi

# ---------------------------------------------------------------------------
# Guard the mdp set against config.sh. There is ONE mdp directory shared by all
# six systems -- that is what "a standard protocol for all" means here -- so the
# only way the protocol can drift is if the mdp files and config.sh stop
# agreeing. Values are READ from the mdp, never assumed: a guard that hardcodes
# dt = 2 fs while the mdp says otherwise is worse than no guard at all.
mdpval () {   # mdpval <file> <key>
  awk -F';' '{print $1}' "$1" \
    | awk -v k="$2" -F'=' '
        { key=$1; gsub(/^[ \t]+|[ \t]+$/, "", key)
          if (tolower(key) == tolower(k)) {
            v=$2; gsub(/^[ \t]+|[ \t]+$/, "", v); gsub(/[ \t]+/, " ", v)
            print v; exit } }'
}

bad=0
for stage in min nvt npt prod; do
  [ -s "mdp/${stage}.mdp" ] || { echo "ERROR: mdp/${stage}.mdp missing."; bad=1; }
done
[ "${bad}" -eq 0 ] || exit 1

DT=$(mdpval mdp/prod.mdp dt)
[ -n "${DT}" ] || { echo "ERROR: no dt in mdp/prod.mdp."; exit 1; }

check_steps () {   # check_steps <file> <wanted_steps> <label>
  local got; got=$(mdpval "$1" nsteps)
  if [ "${got}" != "$2" ]; then
    echo "ERROR: $1 nsteps=${got}, but config.sh says $3 -> ${2} steps at dt=${DT}."
    bad=1
  fi
}
check_steps mdp/nvt.mdp  "$(awk -v ps="${NVT_PS}"  -v dt="${DT}" 'BEGIN{printf "%d", ps/dt}')"      "NVT_PS=${NVT_PS} ps"
check_steps mdp/npt.mdp  "$(awk -v ns="${NPT_NS}"  -v dt="${DT}" 'BEGIN{printf "%d", ns*1000/dt}')" "NPT_NS=${NPT_NS} ns"
check_steps mdp/prod.mdp "$(awk -v ns="${PROD_NS}" -v dt="${DT}" 'BEGIN{printf "%d", ns*1000/dt}')" "PROD_NS=${PROD_NS} ns"

want_out=$(awk -v ps="${FRAME_PS}" -v dt="${DT}" 'BEGIN{printf "%d", ps/dt}')
for stage in nvt npt prod; do
  got=$(mdpval "mdp/${stage}.mdp" nstxout-compressed)
  [ "${got}" = "${want_out}" ] || {
    echo "ERROR: mdp/${stage}.mdp nstxout-compressed=${got}, but config.sh"
    echo "       FRAME_PS=${FRAME_PS} ps means ${want_out} steps at dt=${DT}."
    echo "       Mismatched frame spacing between systems makes them impossible"
    echo "       to put through one analysis, which is the point of running them."
    bad=1; }
done

for stage in nvt npt prod; do
  for v in $(mdpval "mdp/${stage}.mdp" ref_t); do
    [ "${v%%.*}" = "${TEMP_K%%.*}" ] || {
      echo "ERROR: mdp/${stage}.mdp ref_t=${v} but config.sh TEMP_K=${TEMP_K}."; bad=1; }
  done
done

# Production must be unrestrained. A stray define here restrains the whole 100 ns
# and the RMSF comes out near zero -- silent, expensive, plausible-looking.
d=$(mdpval mdp/prod.mdp define)
[ -z "${d}" ] || { echo "ERROR: mdp/prod.mdp sets define='${d}'. Production must be UNRESTRAINED."; bad=1; }
# ...and the two equilibration stages must actually be restrained.
mdpval mdp/nvt.mdp define | grep -q '\-DPOSRES\b'    || { echo "ERROR: mdp/nvt.mdp does not define -DPOSRES."; bad=1; }
mdpval mdp/npt.mdp define | grep -q '\-DPOSRES_BB\b' || { echo "ERROR: mdp/npt.mdp does not define -DPOSRES_BB."; bad=1; }

if [ "${bad}" -ne 0 ]; then
  echo ""
  echo "Nothing submitted. Six systems that do not share one protocol cannot"
  echo "answer the question this directory exists to ask."
  exit 1
fi
echo "Checked: one shared mdp set, ${PROD_NS} ns production at ${FRAME_PS} ps/frame,"
echo "         ${NVT_PS} ps NVT + ${NPT_NS} ns NPT, ${TEMP_K} K, dt = ${DT} ps."
echo ""

# ---------------------------------------------------------------------------
# Is a job for this replicate already in the queue? The run-output guard below
# only sees replicates that have STARTED. A replicate that is queued but has not
# begun has no .cpt, no .gro and no job log -- so without this check, re-running
# submit.sh (after a partial rejection, say) would queue a SECOND job for every
# pending replicate, and two mdruns in one directory corrupt each other's
# trajectory and checkpoints.
queued () {   # queued <job-name>
  command -v squeue >/dev/null 2>&1 || return 1
  [ -n "$(squeue -h -u "${USER:-$(id -un)}" -n "$1" -o '%i' 2>/dev/null)" ]
}

submitted=0
skipped=0
failed=""
for SYS in ${TARGETS}; do
  # shellcheck source=/dev/null
  source "${SCRIPT_DIR}/systems/${SYS}.conf"
  BUILD="build/${SYS}"

  NEED="system.gro system.top posre.itp posre_backbone.itp"
  [ "${HAS_LIGAND}" = "1" ] && NEED="${NEED} glycan.itp posre_lig.itp"
  missing=0
  for f in ${NEED}; do
    [ -f "${BUILD}/${f}" ] || { echo "  MISSING ${BUILD}/${f}"; missing=1; }
  done
  if [ "${missing}" -ne 0 ]; then
    echo "SKIPPING ${SYS}: not built. Run  ./prep_system.sh ${SYS}"
    echo ""
    skipped=$((skipped + 1))
    continue
  fi
  for dd in POSRES POSRES_BB; do
    grep -q "#ifdef ${dd}" "${BUILD}/system.top" || {
      echo "SKIPPING ${SYS}: system.top has no '#ifdef ${dd}' block -- the"
      echo "  restraint file exists but is not #included, so it would silently"
      echo "  do nothing."
      missing=1; }
  done
  if [ "${HAS_LIGAND}" = "1" ]; then
    grep -qE '^[[:space:]]*#include[[:space:]]+"glycan\.itp"' "${BUILD}/system.top" || {
      echo "SKIPPING ${SYS}: holo system whose system.top does not include"
      echo "  glycan.itp. It would run as apo under a holo name."
      missing=1; }
  fi
  if [ "${missing}" -ne 0 ]; then echo ""; skipped=$((skipped + 1)); continue; fi

  echo "${SYS}: submitting ${NREPS} replicate(s), ~$(awk -v a="${NVT_PS}" -v b="${NPT_NS}" \
        -v c="${PROD_NS}" 'BEGIN{printf "%.1f", a/1000+b+c}') ns each"

  for i in $(seq 1 "${NREPS}"); do
    REPDIR="runs/${SYS}/rep${i}"

    # Do NOT overwrite a replicate that has already started. Re-running this over
    # a live run would copy system.gro/system.top out from under a running mdrun
    # AND queue a second job writing the same files -- two mdruns in one
    # directory corrupt each other's trajectory and checkpoints.
    # Look for RUN artifacts specifically: a prepared but never started replicate
    # already holds a copied system.gro. And `ls a b` is unusable here -- it exits
    # non-zero whenever ANY operand is missing, so it would report "clean" for a
    # directory holding a checkpoint but no .gro, which is exactly the mid-run
    # state this guard exists to catch.
    if queued "${SYS}_rep${i}"; then
      echo "  SKIP rep${i}: a job named ${SYS}_rep${i} is already in the queue."
      echo "       Submitting again would put two mdruns in one directory."
      if [ "${FORCE:-0}" != "1" ]; then skipped=$((skipped + 1)); continue; fi
      echo "       FORCE=1 set -- submitting anyway. This is almost certainly wrong."
    fi

    started=$(find "${REPDIR}" -maxdepth 1 \( -name '*.cpt' -o -name 'min.gro' \
                -o -name 'nvt.gro' -o -name 'npt.gro' -o -name 'prod_1.gro' \
                -o -name 'gmx_*.log' \) -print -quit 2>/dev/null || true)
    if [ -n "${started}" ]; then
      echo "  SKIP rep${i}: already contains run output."
      echo "       continue it: (cd ${REPDIR} && sbatch --job-name=${SYS}_rep${i} ${SCRIPT_DIR}/run_md.batch)"
      echo "       start over : rm -rf ${REPDIR}   then re-run this script"
      if [ "${FORCE:-0}" != "1" ]; then skipped=$((skipped + 1)); continue; fi
      echo "       FORCE=1 set -- overwriting."
    fi

    mkdir -p "${REPDIR}"
    cp "${BUILD}/system.gro" "${BUILD}/system.top" "${REPDIR}/"
    cp "${BUILD}/"*.itp "${REPDIR}/" 2>/dev/null || true
    # Command-line options override the #SBATCH directives inside run_md.batch,
    # so config.sh really is the single source of truth for the SLURM request.
    # run_md.batch re-applies the same three on every self-resubmit.
    if (cd "${REPDIR}" && sbatch --job-name="${SYS}_rep${i}" \
         --account="${ACCOUNT}" \
         --partition="${PARTITION_GPU:-general-gpu}" \
         --gres="${GRES:-shard:a100-sxm4:1}" \
         --nodes=1 --ntasks=1 --cpus-per-task="${GPU_CPUS:-8}" \
         --mem="${MEM:-16G}" \
         --export=ALL,EG1_CHAIN=1 "${SCRIPT_DIR}/run_md.batch"); then
      submitted=$((submitted + 1))
    else
      # A rejected submission must NOT pass quietly. With 18 chains going out at
      # once, a queue or QOS limit can refuse some of them, and the error scrolls
      # past between the ones that succeeded. Silently ending up with two
      # replicates instead of three -- and not noticing for days -- is exactly
      # the failure this whole directory is built to avoid.
      failed="${failed}${SYS}/rep${i} "
    fi
  done
  echo ""
done

echo "=========================================================="
echo "submitted ${submitted} replicate chain(s); skipped ${skipped}."
if [ -n "${failed}" ]; then
  echo ""
  echo "**********************************************************************"
  echo "*** $(echo ${failed} | wc -w | tr -d ' ') SUBMISSION(S) WERE REJECTED BY THE SCHEDULER."
  echo "***"
  for r in ${failed}; do echo "***   ${r}"; done
  echo "***"
  echo "*** Those replicates are NOT queued and will never run. The usual cause"
  echo "*** is a per-user job or QOS limit -- check with:"
  echo "***   sacctmgr show assoc user=\${USER} format=user,maxjobs,maxsubmit"
  echo "*** Re-run this script to queue just the missing ones -- it skips any"
  echo "*** replicate that already has run output OR an entry in the queue, so"
  echo "*** it will not double-submit the ones that went through."
  echo "**********************************************************************"
fi
echo ""
echo "Each chain re-queues itself until its ${PROD_NS} ns is finished -- you do"
echo "not resubmit anything by hand."
echo ""
echo "Monitor:"
echo "  squeue -u \$USER"
echo "  grep -h 'next job resumes here\|ALL STAGES COMPLETE' runs/*/rep*/gmx_*.log | tail -30"
echo ""
echo "Check on the FIRST job that starts (both matter):"
echo "  grep 'wall budget' runs/*/rep1/gmx_*.log     # SLURM's real limit, reconciled"
echo "  grep -A2 'Performance' runs/*/rep1/gmx_*.log # ns/day -> is NTASKS=${NTASKS} right?"
echo "  grep -iE 'gen.seed' runs/*/rep*/nvt.log      # a DIFFERENT number per replicate"
