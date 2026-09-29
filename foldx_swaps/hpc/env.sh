# Shared settings for the FoldX jobs on Hellbender.
FOLDX_DIR=/home/yjx9t/data/EG1/foldx
FOLDX=$FOLDX_DIR/foldx_20261231
NRUNS=5   # BuildModel repeats per mutant (side-chain packing is stochastic); report each run

[ -x "$FOLDX" ] || { echo "FATAL: $FOLDX not found or not executable (chmod +x?)"; exit 1; }

# FoldX looks for rotabase.txt in the working directory. Link it in when the
# FoldX folder ships one (newer builds embed it and don't need the file).
link_rotabase() {
    [ -f "$FOLDX_DIR/rotabase.txt" ] && ln -sf "$FOLDX_DIR/rotabase.txt" rotabase.txt
    return 0
}
