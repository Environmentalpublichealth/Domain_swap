#!/usr/bin/env bash
# Submit RepairPDB (2 tasks) -> BuildModel (10 tasks) -> summary on Hellbender.
# Run from anywhere:  bash foldx_swaps/hpc/submit.sh
# Extra sbatch options (e.g. --account=xxx) are passed through:  bash hpc/submit.sh --account=xxx
set -euo pipefail
cd "$(dirname "$0")/.."          # foldx_swaps; all jobs run from here
[ -s mutants.tsv ] || { echo "mutants.tsv missing -- run: python3 make_foldx_inputs.py"; exit 1; }
n=$(( $(wc -l < mutants.tsv) - 1 ))
mkdir -p hpc/logs
rep=$(sbatch --parsable "$@" hpc/repair.sbatch)
bld=$(sbatch --parsable "$@" --array=0-$((n - 1)) --dependency=afterok:$rep hpc/buildmodel.sbatch)
sum=$(sbatch --parsable "$@" --dependency=afterany:$bld hpc/summary.sbatch)
echo "repair=$rep (2 tasks)  buildmodel=$bld ($n tasks)  summary=$sum"
echo "monitor: squeue -u \$USER ; results: foldx_swaps/foldx_summary.tsv"
