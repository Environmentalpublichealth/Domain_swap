#!/bin/bash
#SBATCH --job-name=fetch_uniref30
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --output=fetch_uniref30_%j.out
#SBATCH --error=fetch_uniref30_%j.err

echo "Starting UniRef30 download"

# 3. Download the database
# The -c flag ensures that if the connection drops, wget will resume where it left off
wget -c http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/UniRef30_2020_06_hhsuite.tar.gz

echo "Download completed"

