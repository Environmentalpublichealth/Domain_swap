#!/usr/bin/env python3
"""Per-position analysis of the single-substitution FoldX scan.

Reads singles/<mutant>/Dif_<BB>_Repair.fxout (every substitution of each swap
built on its own) and mutants/<mutant>/Dif_<BB>_Repair.fxout (whole swap).
Runs are reported as min-max ranges, never averaged.

Per substitution:
  total / steric / packing   as in compare_mutants.py
  phi                        backbone phi in the repaired host (Gly with phi > 0 flagged)
  class                      HOTSPOT   every run > 2 kcal/mol
                             tolerated every run < 1 kcal/mol
                             stabilizing every run < -0.5 kcal/mol
                             moderate  otherwise
Per swap:
  sum of singles vs whole swap, run by run (run r of each single summed; runs
  are interchangeable replicates, so the pairing by index is arbitrary but
  gives five independent sums). combined - sum > 0 means the substitutions
  hurt more together than alone (they interact); < 0 means they buffer.
  Partial-swap suggestion: the swap with every HOTSPOT reverted to the host.

Standard library only:  python3 compare_singles.py
Writes foldx_singles_runs.tsv (one row per substitution per run).
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# "backbone clash" is reported by FoldX but not included in "total energy", so it is left out
STERIC = ("Van der Waals clashes", "torsional clash")
HOT, TOL, STAB = 2.0, 1.0, -0.5


def read_dif(path):
    """{(line_index, run_index): {term: value}} from a BuildModel Dif file."""
    lines = open(path).read().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("Pdb\t"))
    head = lines[i].split("\t")
    out = {}
    for l in lines[i + 1:]:
        if not l.strip():
            continue
        v = l.split("\t")
        m = re.search(r"_(\d+)_(\d+)\.pdb$", v[0])
        out[(int(m.group(1)), int(m.group(2)))] = {head[j]: float(v[j]) for j in range(1, len(head)) if v[j] != ""}
    return out


def terms(d):
    steric = sum(d.get(k, 0.0) for k in STERIC)
    return d["total energy"], steric, d["total energy"] - steric - d.get("disulfide", 0.0), d.get("disulfide", 0.0)


def phis(pdb):
    bb = {}
    for l in open(pdb):
        if l.startswith("ATOM") and l[12:16].strip() in ("N", "CA", "C"):
            bb.setdefault(int(l[22:26]), {})[l[12:16].strip()] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]

    def dih(p0, p1, p2, p3):
        b0 = [p0[i] - p1[i] for i in range(3)]
        b1 = [p2[i] - p1[i] for i in range(3)]
        b2 = [p3[i] - p2[i] for i in range(3)]
        n = math.sqrt(sum(x * x for x in b1))
        b1 = [x / n for x in b1]
        v = [b0[i] - sum(b0[j] * b1[j] for j in range(3)) * b1[i] for i in range(3)]
        w = [b2[i] - sum(b2[j] * b1[j] for j in range(3)) * b1[i] for i in range(3)]
        c = [b1[1] * v[2] - b1[2] * v[1], b1[2] * v[0] - b1[0] * v[2], b1[0] * v[1] - b1[1] * v[0]]
        return math.degrees(math.atan2(sum(c[i] * w[i] for i in range(3)), sum(v[i] * w[i] for i in range(3))))

    return {r: dih(bb[r - 1]["C"], bb[r]["N"], bb[r]["CA"], bb[r]["C"]) for r in bb if r - 1 in bb}


def rng(v):
    return f"{min(v):6.1f} to {max(v):5.1f}"


flags = {}
if os.path.exists(os.path.join(HERE, "af3_segment_check.tsv")):
    for l in open(os.path.join(HERE, "af3_segment_check.tsv")).read().splitlines()[1:]:
        flags[l.split("\t")[0]] = l.split("\t")[-1]

rows = [l.split("\t") for l in open(os.path.join(HERE, "mutants.tsv")).read().splitlines()[1:]]
phi = {bb: phis(os.path.join(HERE, "repair", f"{bb}_Repair.pdb")) for bb in ("1XYN", "1KS5")}

missing = [r[0] for r in rows if not os.path.exists(os.path.join(HERE, "singles", r[0], f"Dif_{r[1]}_Repair.fxout"))]
if missing:
    sys.exit("ERROR: no singles Dif file for " + ", ".join(missing) + " -- scan not finished or not copied back")

tsv = open(os.path.join(HERE, "foldx_singles_runs.tsv"), "w")
tsv.write("mutant\tbackbone\tsubstitution\tposition\tphi\trun\ttotal\tsteric\tpacking\tdisulfide\tclass\n")

for name, bb, nsubs, span, muts in rows:
    subs = muts.split(",")
    single = read_dif(os.path.join(HERE, "singles", name, f"Dif_{bb}_Repair.fxout"))
    whole = read_dif(os.path.join(HERE, "mutants", name, f"Dif_{bb}_Repair.fxout"))
    runs = sorted({r for _, r in single})
    if len(single) != len(subs) * len(runs):
        sys.exit(f"ERROR: {name}: expected {len(subs)} x {len(runs)} rows in singles Dif, found {len(single)}")

    flag = flags.get(name, "NA")
    print(f"\n=== {name}  ({bb} {span}, {nsubs} substitutions)  backbone flag: {flag}"
          + ("  <- AF3 predicts the segment refolds; fixed-backbone ddG not reliable" if flag == "AF3_refold" else ""))
    print(f"  {'subst':8s} {'phi':>5s}  {'total':>14s}  {'steric':>14s}  {'packing':>14s}  class")
    hot = []
    for i, s in enumerate(subs, 1):
        pos = int(s[2:-1])
        T = [terms(single[(i, r)]) for r in runs]
        tot = [t[0] for t in T]
        cls = ("HOTSPOT" if min(tot) > HOT else "stabilizing" if max(tot) < STAB
               else "tolerated" if max(tot) < TOL else "moderate")
        note = ""
        if s[0] == "G" and phi[bb].get(pos, -1) > 0:
            note = "  Gly with phi>0"
        if s[0] == "C":
            note = "  removes Cys"
        if cls == "HOTSPOT":
            hot.append(s)
        print(f"  {s:8s} {phi[bb].get(pos, float('nan')):5.0f}  {rng(tot)}  {rng([t[1] for t in T])}  "
              f"{rng([t[2] for t in T])}  {cls}{note}")
        for r, t in zip(runs, T):
            tsv.write(f"{name}\t{bb}\t{s}\t{pos}\t{phi[bb].get(pos, float('nan')):.0f}\t{r + 1}\t"
                      f"{t[0]:.2f}\t{t[1]:.2f}\t{t[2]:.2f}\t{t[3]:.2f}\t{cls}\n")

    sums = [sum(terms(single[(i, r)])[0] for i in range(1, len(subs) + 1)) for r in runs]
    comb = [terms(whole[(1, r)])[0] for r in sorted({r for _, r in whole})]
    diff = [c - s for c, s in zip(comb, sums)]
    print(f"  sum of singles {rng(sums)}   whole swap {rng(comb)}   whole - sum {rng(diff)}"
          + ("   -> substitutions interact" if min(abs(d) for d in diff) > 2 else "   -> roughly additive"))
    if hot:
        keep = [s for s in subs if s not in hot]
        print(f"  partial swap without hotspots ({len(keep)}/{len(subs)} kept; sum of singles "
              f"{rng([sum(terms(single[(subs.index(s) + 1, r)])[0] for s in keep) for r in runs])}): "
              + (",".join(keep) if keep else "nothing left"))
tsv.close()
print("\nper-run table: foldx_singles_runs.tsv")
print("Partial-swap sums assume additivity; build the partial swap itself before trusting it.")
