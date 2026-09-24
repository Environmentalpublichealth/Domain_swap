#!/bin/bash
# ============================================================================
# Submit run_analysis.sh for every finished replicate.
#
#   ./submit_analysis.sh                  all systems in config.sh SYSTEMS
#   ./submit_analysis.sh 1xyn_apo q4_holo just these
#   B_NS=50 ./submit_analysis.sh          discard the first 50 ns from RMSD/RMSF
#   EXTRAS=1 ./submit_analysis.sh         also Rg + DSSP
#   FORCE=1 ./submit_analysis.sh          re-analyse replicates already done
#
# One job per replicate. They are independent analyses of independent
# simulations and are reported separately, never averaged together.
# ============================================================================
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
[ -r "${SCRIPT_DIR}/config.sh" ] && source "${SCRIPT_DIR}/config.sh"

# config.sh leaves SYSTEMS empty on purpose: the systems/ directory IS the list,
# so ./new_system.sh is the only place a system gets declared. Setting SYSTEMS in
# config.sh narrows "all" to that subset.
if [ -z "${SYSTEMS:-}" ]; then
  SYSTEMS="$(ls "${SCRIPT_DIR}/systems" 2>/dev/null | sed 's/\.conf$//' | tr '\n' ' ')"
fi

command -v sbatch >/dev/null 2>&1 || { echo "ERROR: no sbatch on PATH."; exit 1; }

LIST=${*:-${SYSTEMS:-}}
[ -n "${LIST}" ] || { echo "ERROR: no systems given and config.sh sets no SYSTEMS."; exit 1; }

sub=0; skip=0; miss=0
for SYS in ${LIST}; do
  [ -f "${SCRIPT_DIR}/systems/${SYS}.conf" ] || { echo "  ${SYS}: no systems/${SYS}.conf -- skipping"; continue; }
  for i in $(seq 1 "${N_REPS:-3}"); do
    REPDIR="${SCRIPT_DIR}/runs/${SYS}/rep${i}"
    TAG="${SYS}_rep${i}_c36m"

    if [ ! -s "${REPDIR}/prod_1.xtc" ]; then
      echo "  ${TAG}  no prod_1.xtc -- not started, skipping"; miss=$((miss+1)); continue
    fi
    # An unfinished run is still worth analysing; say so rather than skipping it.
    [ -s "${REPDIR}/prod_1.gro" ] || echo "  ${TAG}  NOTE: production unfinished, analysing what exists"
    if [ -s "${REPDIR}/rmsf_amideH_${TAG}.xvg" ] && [ "${FORCE:-0}" != "1" ]; then
      echo "  ${TAG}  already analysed -- skipping (FORCE=1 to redo)"; skip=$((skip+1)); continue
    fi

    jid=$(cd "${REPDIR}" && sbatch --parsable \
            --export=ALL,B_NS="${B_NS:-0}",EXTRAS="${EXTRAS:-0}" \
            --account="${ACCOUNT}" \
            --partition="${PARTITION_CPU:-general-cpu}" \
            --cpus-per-task="${PREP_CPUS:-4}" \
            --job-name="an_${SYS}_rep${i}" \
            "${SCRIPT_DIR}/run_analysis.sh")
    echo "  ${TAG}  submitted job ${jid}"
    sub=$((sub+1))
  done
done

echo ""
echo "${sub} submitted, ${skip} already done, ${miss} not started."
echo ""
echo "watch:    squeue -u ${USER}"
echo "check:    grep -l 'DONE' runs/*/rep*/analysis_out_*.log | wc -l"
echo "failures: grep -h 'FATAL\|WARNING' runs/*/rep*/analysis_out_*.log | sort | uniq -c"
