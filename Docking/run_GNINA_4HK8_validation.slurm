#!/bin/bash
#SBATCH -p general-gpu
#SBATCH -A engr-lab-joshua.yuan
#SBATCH --gres="shard:a6000:1"
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH -J GNINA_4HK8
#SBATCH --output=GNINA_4HK8_%j.out
#SBATCH --error=GNINA_4HK8_%j.err

set -euo pipefail

BASE=/engrfs/project/joshua.yuan/tanyuxuan/domain_swap
SIF=$BASE/containers/gnina.sif
WORK=$BASE/docking/validation/4HK8
OUT=$WORK/results/exh32/reproducibility

mkdir -p "$OUT"

# 4HK8-X6 redocking validation
# GNINA v1.3.3
# Experimental X6 defines the docking box
# exhaustiveness = 32
# 4 independent random seeds
# 9 output poses per seed

for SEED in 0 1 2 3
do
    echo "=== 4HK8-X6 | exhaustiveness 32 | seed $SEED ==="

    apptainer exec --nv \
        --bind "$BASE:$BASE" \
        "$SIF" gnina \
        -r "$WORK/4HK8_receptor_clean.pdb" \
        -l "$WORK/X6.sdf" \
        --autobox_ligand "$WORK/X6.sdf" \
        --seed "$SEED" \
        --cpu 1 \
        --exhaustiveness 32 \
        --num_modes 9 \
        -o "$OUT/4HK8_X6_exh32_seed${SEED}.sdf.gz"
done
