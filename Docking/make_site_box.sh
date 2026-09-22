#!/usr/bin/env bash
# Build a pocket PDB for GNINA --autobox_ligand when no reference ligand exists.
# Runs on every *.pdb in this directory except reduce output (*_H.pdb). For each model it finds the catalytic
# Glu pair automatically (the only two Glu whose carboxyl centers are within
# PAIR_MAX A of each other), then selects every residue with any atom within
# CUTOFF A of the midpoint of that pair -- in the model's OWN coordinate frame
# (AF3 models are not superimposed, so one shared box would miss the pocket).
# A model with zero or several candidate pairs is flagged and skipped, never
# guessed.
# Usage: bash make_site_box.sh [CUTOFF]   (default 10 A)
# Output: autobox/<model>_site.pdb ; pass it to gnina --autobox_ligand

set -euo pipefail
cd "$(dirname "$0")"
CUTOFF=${1:-10}
PAIR_MAX=7.5   # retaining GH11/GH12 pairs here are 6.3-6.9 A apart
mkdir -p autobox

shopt -s nullglob extglob
# skip reduce output (*_H.pdb): same coordinates as the heavy-atom model
pdbs=(!(*_H).pdb)
[ ${#pdbs[@]} -gt 0 ] || { echo "ERROR: no .pdb files in $(pwd)"; exit 1; }

failed=0
for pdb in "${pdbs[@]}"; do
  model=${pdb%.pdb}

  # candidate catalytic pairs: "chain resA resB dist", one line per pair
  pairs=$(awk -v max="$PAIR_MAX" '
    /^ATOM/ && substr($0,18,3)=="GLU" && (substr($0,13,4)==" OE1" || substr($0,13,4)==" OE2") {
      k = substr($0,22,1) SUBSEP (substr($0,23,4)+0)
      n[k]++; x[k]+=substr($0,31,8); y[k]+=substr($0,39,8); z[k]+=substr($0,47,8)
    }
    END {
      for (a in n) for (b in n) {
        split(a, A, SUBSEP); split(b, B, SUBSEP)
        if (A[1] != B[1] || A[2]+0 >= B[2]+0 || n[a] != 2 || n[b] != 2) continue
        d = sqrt((x[a]/2-x[b]/2)^2 + (y[a]/2-y[b]/2)^2 + (z[a]/2-z[b]/2)^2)
        if (d <= max) printf "%s %d %d %.1f\n", A[1], A[2], B[2], d
      }
    }' "$pdb")

  npairs=$(printf "%s" "$pairs" | grep -c . || true)
  if [ "$npairs" -ne 1 ]; then
    echo "SKIP  $model: found $npairs Glu pairs within $PAIR_MAX A (need exactly 1)${pairs:+:}"
    [ -n "$pairs" ] && printf "%s\n" "$pairs" | sed 's/^/        chain /'
    failed=1
    continue
  fi
  read -r chain e1 e2 dist <<< "$pairs"

  awk -v c="$chain" -v a="$e1" -v b="$e2" -v cut="$CUTOFF" '
    # pass 1: midpoint of the four OE1/OE2 atoms of the catalytic pair
    NR==FNR {
      if (/^ATOM/ && substr($0,18,3)=="GLU" && substr($0,13,3)==" OE" && substr($0,22,1)==c) {
        r = substr($0,23,4)+0
        if (r==a || r==b) { cx+=substr($0,31,8); cy+=substr($0,39,8); cz+=substr($0,47,8); n++ }
      }
      next
    }
    FNR==1 { cx/=n; cy/=n; cz/=n }
    # pass 2: flag residues with any atom within cutoff
    /^ATOM/ {
      r = substr($0,22,1) substr($0,23,4); line[++m] = $0; res[m] = r
      d = sqrt((substr($0,31,8)-cx)^2 + (substr($0,39,8)-cy)^2 + (substr($0,47,8)-cz)^2)
      if (d <= cut) keep[r] = 1
    }
    END {
      for (i=1; i<=m; i++) if (res[i] in keep) print line[i]
      print "END"
    }' "$pdb" "$pdb" > "autobox/${model}_site.pdb"

  # report the box GNINA will build (extent of the site atoms + default 4 A padding)
  awk -v model="$model" -v pair="GLU${e1}/GLU${e2} (${dist} A)" '
    /^ATOM/ { x=substr($0,31,8)+0; y=substr($0,39,8)+0; z=substr($0,47,8)+0
      if (!n++) { x0=x1=x; y0=y1=y; z0=z1=z }
      if (x<x0) x0=x; if (x>x1) x1=x; if (y<y0) y0=y; if (y>y1) y1=y; if (z<z0) z0=z; if (z>z1) z1=z
      r[substr($0,22,5)] = 1 }
    END { nr=0; for (k in r) nr++
      printf "OK    %-24s %-22s %3d res  center %6.1f %6.1f %6.1f  box %4.1f x %4.1f x %4.1f A\n",
        model, pair, nr, (x0+x1)/2, (y0+y1)/2, (z0+z1)/2, x1-x0+8, y1-y0+8, z1-z0+8 }' \
    "autobox/${model}_site.pdb"
done

[ "$failed" -eq 0 ] || { echo "Some models were skipped -- see SKIP lines above."; exit 1; }
