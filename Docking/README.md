# Run docking with GNINA
## 1. Prepare receptor files
Download AF3 prediction from the alphafold server. The download folder includes 5 folds ranking from 0 to 4, we choose the best one (rank 0). All the coordinate files are in CIF format, pymol can open it but we need pdb as inputs so we convert CIF into PDB format using pymol export. 

The AF3 structures are mostly clean, but it is missing hydrogen. For docking, we need to add hydrogens into the pdb files. The program 'reduce' can do it.
```bash
# enter the vina conda environment
conda activate vina # you should see (vina)
conda install bioconda::reduce
# add H to the structure
cd /path/to/the/folder/store/all/pdb/
reduce -FLIP Eg1A_WT_AF3_model_0.pdb > Eg1A_WT_AF3_model_0_H.pdb # let's take the EG structure as an example, change the file name for other proteins
```
There should be a output file `Eg1A_WT_AF3_model_0_H.pdb` produced in the same directory.

## 2. Create box site for docking
Each structure should have their own box site to dock, we create it by running the script in here.
```
# download the make_site_box.sh script
wget 
make_site_box.sh
```
This will output a folder named `autobox`, and it should have one _site.pdb file for each protein structure. 

The script will loop over all the pdb files in the working directory. 

## 3. Run GNINA
This script will loop over all the pdb files in the working directory. 
```bash
BASE=/path/to/yout/docking/folder # change this to your own path
SIF=$BASE/Docking/gnina_build/gnina.sif
WORK=$BASE/Docking
OUT=$WORK/results

mkdir -p "$OUT"
for PDB in "$WORK"/*.pdb; do
    M=$(basename "$PDB" .pdb) # define the pdb structure
    BOX="$WORK/autobox/${M}_site.pdb" # define the box files
    apptainer exec --nv \
        --bind "$BASE:$BASE" \
        "$SIF" gnina \
        -r "$WORK/${M}_H.pdb" \
        -l "$WORK/X6.sdf" \
        --autobox_ligand "$BOX" \
        --seed 42 \
        --cpu 4 \
        --exhaustiveness 32 \
        --num_modes 9 \
        -o "$OUT/${M}_X6_exh32.sdf.gz"
done
```
Make sure all the file paths are correct. 
