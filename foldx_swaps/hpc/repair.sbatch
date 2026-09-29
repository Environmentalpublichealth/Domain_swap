#!/usr/bin/env bash
#SBATCH --job-name=foldx_repair
#SBATCH --partition=requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --array=0-1
#SBATCH --output=hpc/logs/%x_%A_%a.out
# RepairPDB on each cleaned crystal structure: task 0 -> 1XYN, task 1 -> 1KS5.
set -euo pipefail
source hpc/env.sh

PDBS=(1XYN 1KS5)
P=${PDBS[$SLURM_ARRAY_TASK_ID]}

cd repair
[ -s "${P}.pdb" ] || { echo "FATAL: repair/${P}.pdb missing -- run make_foldx_inputs.py first"; exit 1; }
if [ -s "${P}_Repair.pdb" ]; then echo "exists, skipping: repair/${P}_Repair.pdb"; exit 0; fi
link_rotabase

echo "$(date '+%F %T')  RepairPDB $P on $(hostname)"
"$FOLDX" --command=RepairPDB --pdb="${P}.pdb" > "${P}_repair.log" 2>&1
[ -s "${P}_Repair.pdb" ] || { echo "FATAL: no ${P}_Repair.pdb -- see repair/${P}_repair.log"; tail -20 "${P}_repair.log"; exit 1; }
echo "$(date '+%F %T')  wrote repair/${P}_Repair.pdb"
