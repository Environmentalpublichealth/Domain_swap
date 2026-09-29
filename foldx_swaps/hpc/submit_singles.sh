#!/usr/bin/env bash
# Single-substitution scan: BuildModel on singles/<mutant>/ (10 tasks, one per swap;
# each task builds that swap's substitutions one at a time). Reuses repair/*_Repair.pdb
# from the first run, so no repair job is submitted.
# Run from anywhere:  bash foldx_swaps/hpc/submit_singles.sh [--account=xxx]
set -euo pipefail
cd "$(dirname "$0")/.."          # foldx_swaps; all jobs run from here
[ -d singles ] || { echo "singles/ missing -- run: python3 make_foldx_inputs.py"; exit 1; }
for p in 1XYN 1KS5; do
    [ -s repair/${p}_Repair.pdb ] || { echo "repair/${p}_Repair.pdb missing -- run hpc/submit.sh first"; exit 1; }
done
n=$(( $(wc -l < mutants.tsv) - 1 ))
mkdir -p hpc/logs
bld=$(sbatch --parsable "$@" --job-name=foldx_singles --array=0-$((n - 1)) \
      --export=ALL,SET=singles hpc/buildmodel.sbatch)
echo "singles=$bld ($n tasks, $(cat singles/*/individual_list.txt | wc -l | tr -d ' ') substitutions)"
echo "monitor: squeue -u \$USER ; when done copy foldx_swaps/singles/ back and run: python3 compare_singles.py"
