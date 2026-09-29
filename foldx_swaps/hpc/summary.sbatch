#!/usr/bin/env bash
#SBATCH --job-name=foldx_summary
#SBATCH --partition=requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --output=hpc/logs/%x_%j.out
# Collect the total-energy ddG of every BuildModel run into foldx_summary.tsv.
# Each run is reported on its own row (no averaging), so run-to-run spread stays visible.
set -euo pipefail

out=foldx_summary.tsv
printf "mutant\tbackbone\tn_subs\trun\tddG_kcal_mol\n" > "$out"
tail -n +2 mutants.tsv | while IFS=$'\t' read -r MUT BB NSUBS SPAN MUTS; do
    dif=mutants/$MUT/Dif_${BB}_Repair.fxout
    if [ ! -s "$dif" ]; then
        printf "%s\t%s\t%s\tNA\tMISSING\n" "$MUT" "$BB" "$NSUBS" >> "$out"; continue
    fi
    # Rows look like "<BB>_Repair_1_<run>.pdb<TAB><total ddG><TAB>...", after the "Pdb" header line
    awk -F'\t' -v m="$MUT" -v b="$BB" -v n="$NSUBS" '
        $1 == "Pdb" { go = 1; next }
        go && $1 ~ /\.pdb$/ { r = $1; sub(/\.pdb$/, "", r); sub(/.*_/, "", r)
                              printf "%s\t%s\t%s\t%d\t%.2f\n", m, b, n, r + 1, $2 }' "$dif" >> "$out"
done

column -t -s $'\t' "$out"
echo
echo "ddG = mutant - WT total energy (kcal/mol); positive = destabilizing."
