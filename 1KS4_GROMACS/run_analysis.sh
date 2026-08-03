#!/bin/bash
#SBATCH --job-name=1KS4_analysis
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=analysis_out_%j.log
#SBATCH --error=analysis_err_%j.log

# Usage: sbatch run_analysis.sh <CASE> <REP>
#   e.g. sbatch run_analysis.sh WT 1
# Run from inside that case/replicate's own working directory (where prod_1.xtc
# / prod_1.tpr live). Output filenames embed enzyme+case+replicate so results
# from different cases/replicates/proteins never collide -- see the AdhA/XYN
# rmsf_WT.xvg mix-up from earlier for why this matters.

CASE=$1
REP=$2
if [ -z "$CASE" ] || [ -z "$REP" ]; then
  echo "Usage: sbatch run_analysis.sh <CASE> <REP>"
  exit 1
fi

TAG="1KS4_${CASE}_rep${REP}"

module purge
module load gromacs/2022.5_plumed_gcc_9.5.0_openmpi_4.1.5

echo "Starting analysis for ${TAG}..."

# Using group NAMES (not index numbers) below -- group ordering isn't
# guaranteed to match between a pdb2gmx-built topology and an acpype-converted
# one, so a hardcoded index could silently select the wrong atoms.

# A. Center the protein in the box
echo "Protein System" | gmx_mpi trjconv -s prod_1.tpr -f prod_1.xtc -o prod_center.xtc -pbc mol -center

# B. Remove rotation/translation (fit on Backbone)
echo "Backbone System" | gmx_mpi trjconv -s prod_1.tpr -f prod_center.xtc -o prod_final.xtc -fit rot+trans

# C. RMSD (backbone)
echo "Backbone Backbone" | gmx_mpi rms -s prod_1.tpr -f prod_final.xtc -o rmsd_${TAG}.xvg -tu ns

# D. Per-residue RMSF (backbone)
echo "Backbone" | gmx_mpi rmsf -s prod_1.tpr -f prod_final.xtc -o rmsf_${TAG}.xvg -res

# E. Radius of gyration
echo "Protein" | gmx_mpi gyrate -s prod_1.tpr -f prod_final.xtc -o rg_${TAG}.xvg

# F. Secondary structure occupancy (needed to confirm helix/loop loss directly,
#    not just infer it from RMSF -- see conversation notes on why RMSF alone
#    was ambiguous for the AdhA mutation-site comparison).
#    gmx 2022.5 uses do_dssp + an external `mkdssp`/`dssp` binary on PATH
#    (module load a dssp module, or point DSSP env var at the binary).
#    If your GROMACS module is >=2023, replace with: gmx dssp -s ... -f ... -o ...
echo "Protein" | gmx_mpi do_dssp -s prod_1.tpr -f prod_final.xtc -o ss_${TAG}.xpm -sc scount_${TAG}.xvg

echo "Cleaning up intermediate centered file..."
rm -f prod_center.xtc

echo "Pipeline complete: rmsd_${TAG}.xvg, rmsf_${TAG}.xvg, rg_${TAG}.xvg, ss_${TAG}.xpm, scount_${TAG}.xvg"
