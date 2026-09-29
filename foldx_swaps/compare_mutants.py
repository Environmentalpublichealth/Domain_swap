#!/usr/bin/env python3
"""Compare FoldX BuildModel results across the swap mutants.

For every mutant and every run (runs are reported individually, never averaged):
  total      FoldX ddG (mutant - WT built in the same run), kcal/mol
  steric     Van der Waals clashes + torsional clash:
             penalty from forcing the segment onto the fixed host backbone;
             would partly relax if the backbone could move
  disulfide  disulfide term (non-zero only when a Cys is removed)
  packing    total - steric - disulfide: H-bonds, solvation, vdW, entropy;
             packing/polarity mismatch that backbone relaxation would not fix
  per_sub    total / number of substitutions (size-adjusted cost)

Backbone flag (from af3_segment_check.tsv, if present):
  native_backbone     AF3 keeps the segment in the host conformation; ddG usable
  low_AF3_confidence  AF3 keeps it but with segment pLDDT < 70; ddG usable with care
  AF3_refold          AF3 predicts the segment refolds; the fixed-backbone ddG
                      only says "does not fit the native backbone" -- do not
                      trust its value, and it is left out of rankings/comparisons

Standard library only:  python3 compare_mutants.py
Writes foldx_compare_runs.tsv (one row per run) and prints the per-mutant
ranges, the reciprocal-pair table, and which mutants differ beyond run noise.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
# "backbone clash" is reported by FoldX but not included in "total energy", so it is left out
STERIC = ("Van der Waals clashes", "torsional clash")


def read_dif(path):
    lines = open(path).read().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("Pdb\t"))
    head = lines[i].split("\t")
    runs = []
    for l in lines[i + 1:]:
        if not l.strip():
            continue
        v = l.split("\t")
        runs.append({head[j]: float(v[j]) for j in range(1, len(head)) if v[j] != ""})
    return runs


flags = {}
fp = os.path.join(HERE, "af3_segment_check.tsv")
if os.path.exists(fp):
    for l in open(fp).read().splitlines()[1:]:
        v = l.split("\t")
        flags[v[0]] = v[-1]
else:
    print("NOTE: af3_segment_check.tsv not found -- backbone flags not shown\n")

mutants = []
for line in open(os.path.join(HERE, "mutants.tsv")).read().splitlines()[1:]:
    name, bb, nsubs, span, muts = line.split("\t")
    runs = []
    for r, d in enumerate(read_dif(os.path.join(HERE, "mutants", name, f"Dif_{bb}_Repair.fxout")), 1):
        total = d["total energy"]
        steric = sum(d.get(k, 0.0) for k in STERIC)
        ss = d.get("disulfide", 0.0)
        runs.append(dict(run=r, total=total, steric=steric, disulfide=ss,
                         packing=total - steric - ss, per_sub=total / int(nsubs)))
    mutants.append(dict(name=name, bb=bb, swap=name.split("_", 1)[1], nsubs=int(nsubs),
                        span=span, runs=runs, flag=flags.get(name, "NA")))

cols = ("total", "steric", "disulfide", "packing", "per_sub")
with open(os.path.join(HERE, "foldx_compare_runs.tsv"), "w") as f:
    f.write("mutant\tbackbone\tswap\tn_subs\tspan\tbackbone_flag\trun\t" + "\t".join(cols) + "\n")
    for m in mutants:
        for r in m["runs"]:
            f.write(f"{m['name']}\t{m['bb']}\t{m['swap']}\t{m['nsubs']}\t{m['span']}\t{m['flag']}\t{r['run']}\t"
                    + "\t".join(f"{r[c]:.2f}" for c in cols) + "\n")


def rng(m, key):
    v = [r[key] for r in m["runs"]]
    return min(v), max(v)


def fmt(lo, hi):
    return f"{lo:6.1f} to {hi:5.1f}"


def tier(hi):
    return ("tolerated" if hi < 2 else "moderate" if hi < 7 else
            "strong" if hi < 15 else "incompatible")


def driver(m):
    st, pk, ss = (max(r[k] for r in m["runs"]) for k in ("steric", "packing", "disulfide"))
    parts = []
    if st >= 3:
        parts.append("steric")
    if pk >= 3:
        parts.append("packing")
    if ss >= 1:
        parts.append("disulfide")
    return "+".join(parts) or "-"


print("Per-mutant ranges over runs (kcal/mol; min to max, not averaged)\n")
print(f"{'mutant':9s} {'subs':>4s}  {'total':>14s}  {'steric':>14s}  {'packing':>14s}  "
      f"{'per_sub':>14s}  {'SS':>4s}  tier          {'driver':22s}  backbone flag")
for m in sorted(mutants, key=lambda m: rng(m, "total")[1]):
    print(f"{m['name']:9s} {m['nsubs']:4d}  {fmt(*rng(m, 'total'))}  {fmt(*rng(m, 'steric'))}  "
          f"{fmt(*rng(m, 'packing'))}  {fmt(*rng(m, 'per_sub'))}  {rng(m, 'disulfide')[1]:4.1f}  "
          f"{tier(rng(m, 'total')[1]):12s}  {driver(m):22s}  {m['flag']}"
          + ("  <- ddG not reliable" if m["flag"] == "AF3_refold" else ""))

print("\nReciprocal pairs: same segment exchanged in both directions (total ddG ranges)\n")
print(f"{'swap':5s} {'1KS5 segment -> 1XYN':>20s}   {'1XYN segment -> 1KS5':>20s}   better host")
for swap in sorted({m["swap"] for m in mutants}):
    by = {m["bb"]: m for m in mutants if m["swap"] == swap}
    if set(by) != {"1XYN", "1KS5"}:
        continue
    x, k = rng(by["1XYN"], "total"), rng(by["1KS5"], "total")
    better = "1XYN" if x[1] < k[0] else "1KS5" if k[1] < x[0] else "overlap"
    if "AF3_refold" in (by["1XYN"]["flag"], by["1KS5"]["flag"]):
        better = "not comparable (AF3 refold)"
    print(f"{swap:5s} {fmt(*x):>20s}   {fmt(*k):>20s}   {better}")

# Two mutants are called different only if their run ranges are separated by
# more than 1 kcal/mol (roughly FoldX's own error on a single mutation).
GAP = 1.0
print(f"\nWithin each backbone, pairs whose run ranges are NOT separated by >{GAP} kcal/mol"
      " (treat as indistinguishable; AF3_refold mutants excluded):\n")
for bb in ("1KS5", "1XYN"):
    ms = sorted((m for m in mutants if m["bb"] == bb and m["flag"] != "AF3_refold"), key=lambda m: rng(m, "total")[0])
    ties = []
    for i, a in enumerate(ms):
        for b in ms[i + 1:]:
            if rng(b, "total")[0] - rng(a, "total")[1] <= GAP:
                ties.append(f"{a['swap']}~{b['swap']}")
    print(f"  {bb}: {', '.join(ties) if ties else 'none (all separable)'}")
print("\nper-run table: foldx_compare_runs.tsv")
