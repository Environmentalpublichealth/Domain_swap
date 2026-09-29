#!/usr/bin/env python3
"""Build FoldX BuildModel inputs for the 1XYN <-> 1KS5 region-swap mutants.

Every swap replaces a segment with one of the same length, so each mutant is
a set of point substitutions on the crystal backbone. Mutations are derived by
comparing each FASTA sequence with the sequence of its crystal structure; the
region ranges written in the FASTA headers are only reported, never used (the
M03 headers disagree with the M03 sequences).

Standard library only. Run from anywhere:
    python3 make_foldx_inputs.py [path/to/swaps.fasta]

Writes, next to this script:
    repair/<PDB>.pdb                 cleaned crystal chain A (no water/ions/altlocs)
    mutants/<name>/individual_list.txt
    mutants.tsv                      one row per mutant (row order = SLURM array index)
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FASTA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "1KS5_1XYN_region_swap_10_sequences.fasta")
BACKBONES = ["1XYN", "1KS5"]
CHAIN = "A"

THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
         "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
         "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def clean_pdb(src, dst):
    """Keep chain A protein ATOM records only; first altloc; renumber nothing.
    Returns {resnum: one-letter}."""
    seq, out = {}, []
    for line in open(src):
        if not line.startswith("ATOM") or line[21] != CHAIN:
            continue
        if line[16] not in (" ", "A"):
            continue
        if line[17:20] not in THREE:
            sys.exit(f"ERROR: non-standard residue {line[17:20]} in {src}")
        if line[26] != " ":
            sys.exit(f"ERROR: insertion code at residue {line[22:27]} in {src}")
        out.append(line[:16] + " " + line[17:])
        if line[12:16] == " CA ":
            seq[int(line[22:26])] = THREE[line[17:20]]
    open(dst, "w").write("".join(out) + "TER\nEND\n")
    nums = sorted(seq)
    if nums != list(range(nums[0], nums[-1] + 1)):
        sys.exit(f"ERROR: {src} has missing residues; FoldX mutations would be misnumbered")
    return seq


def read_fasta(path):
    recs, name = [], None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            recs.append([line[1:], ""])
        elif line:
            recs[-1][1] += line
    return recs


os.makedirs(os.path.join(HERE, "repair"), exist_ok=True)
wt = {}
for bb in BACKBONES:
    wt[bb] = clean_pdb(os.path.join(HERE, "pdb", f"{bb}.pdb"), os.path.join(HERE, "repair", f"{bb}.pdb"))
    print(f"{bb}: chain {CHAIN}, residues {min(wt[bb])}-{max(wt[bb])} -> repair/{bb}.pdb")

rows, warnings = [], []
mutseq = {}
for header, seq in read_fasta(FASTA):
    name = header.split("|")[0].strip()
    bb = name.split("_")[0]
    if bb not in wt:
        sys.exit(f"ERROR: {name}: backbone {bb} not in {BACKBONES}")
    ref = wt[bb]
    first = min(ref)
    if len(seq) != len(ref):
        sys.exit(f"ERROR: {name}: length {len(seq)} != {bb} length {len(ref)}; "
                 "not a same-length swap, BuildModel cannot model it")
    muts = []
    for i, aa in enumerate(seq):
        pos = first + i
        if aa not in THREE.values():
            sys.exit(f"ERROR: {name}: invalid residue '{aa}' at {pos}")
        if aa != ref[pos]:
            muts.append((pos, ref[pos], aa))
    if not muts:
        sys.exit(f"ERROR: {name}: identical to {bb}")
    span = (muts[0][0], muts[-1][0])
    mutseq[name] = (bb, seq, span)

    hdr = re.search(r"%s:(\d+)-(\d+)" % bb, header)
    hdr_span = (int(hdr.group(1)), int(hdr.group(2))) if hdr else None
    if hdr_span and not (hdr_span[0] <= span[0] and span[1] <= hdr_span[1]):
        warnings.append(f"{name}: header says {bb}:{hdr_span[0]}-{hdr_span[1]}, "
                        f"but the sequence changes {bb}:{span[0]}-{span[1]} (using the sequence)")

    d = os.path.join(HERE, "mutants", name)
    os.makedirs(d, exist_ok=True)
    mut_str = ",".join(f"{w}{CHAIN}{p}{m}" for p, w, m in muts)
    open(os.path.join(d, "individual_list.txt"), "w").write(mut_str + ";\n")
    rows.append((name, bb, len(muts), f"{span[0]}-{span[1]}", mut_str))

# Reciprocal check: 1XYN_Mxx's insert must be 1KS5 wild type over 1KS5_Mxx's span, and vice versa.
for name, (bb, seq, span) in mutseq.items():
    tag = name.split("_", 1)[1]
    other = [n for n in mutseq if n != name and n.split("_", 1)[1] == tag]
    if len(other) != 1:
        warnings.append(f"{name}: no reciprocal partner found for {tag}")
        continue
    obb, _, ospan = mutseq[other[0]]
    first = min(wt[bb])
    insert = seq[span[0] - first: span[1] + 1 - first]
    donor = "".join(wt[obb][p] for p in range(ospan[0], ospan[1] + 1))
    if insert != donor:
        warnings.append(f"{name}: insert {insert} != {obb}:{ospan[0]}-{ospan[1]} {donor}")

with open(os.path.join(HERE, "mutants.tsv"), "w") as f:
    f.write("mutant\tbackbone\tn_subs\tchanged_span\tmutations\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")

print(f"\n{len(rows)} mutants -> mutants.tsv, mutants/<name>/individual_list.txt")
for name, bb, n, span, _ in rows:
    print(f"  {name:10s} {bb}  {n:2d} substitutions  {span}")
if warnings:
    print("\nWARNINGS:")
    for w in warnings:
        print("  " + w)
