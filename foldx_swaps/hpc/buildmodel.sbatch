#!/usr/bin/env bash
#SBATCH --job-name=foldx_build
#SBATCH --partition=requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --array=0-9
#SBATCH --output=hpc/logs/%x_%A_%a.out
# One task per mutant: task i -> data row i of mutants.tsv. FoldX is single-threaded.
set -euo pipefail
source hpc/env.sh

mapfile -t ROWS < <(tail -n +2 mutants.tsv)
IFS=$'\t' read -r MUT BB NSUBS SPAN MUTS <<< "${ROWS[$SLURM_ARRAY_TASK_ID]}"
REP=repair/${BB}_Repair.pdb
[ -s "$REP" ] || { echo "FATAL: $REP missing -- repair job failed?"; exit 1; }

cd "mutants/$MUT"
if [ -s "Dif_${BB}_Repair.fxout" ]; then echo "exists, skipping: mutants/$MUT"; exit 0; fi
ln -sf "../../$REP" "${BB}_Repair.pdb"
link_rotabase

echo "$(date '+%F %T')  BuildModel $MUT on $BB ($NSUBS substitutions, $NRUNS runs) on $(hostname)"
echo "  $MUTS"
"$FOLDX" --command=BuildModel --pdb="${BB}_Repair.pdb" --mutant-file=individual_list.txt \
         --numberOfRuns="$NRUNS" > buildmodel.log 2>&1
[ -s "Dif_${BB}_Repair.fxout" ] || { echo "FATAL: no Dif file -- see mutants/$MUT/buildmodel.log"; tail -20 buildmodel.log; exit 1; }
echo "$(date '+%F %T')  done"
