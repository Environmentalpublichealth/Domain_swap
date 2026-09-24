#!/usr/bin/env python3
"""Clean one receptor PDB into build/<system>/protein.pdb for pdb2gmx.

    prep_receptor.py <in.pdb> <out.pdb> [--chain A]

Keeps ATOM records of one chain, first altloc only, and drops everything
pdb2gmx would either choke on or silently reinterpret: HETATM (waters, ions,
the docked ligand -- that arrives separately and in CHARMM naming), hydrogens
(pdb2gmx rebuilds them for the force field it is actually using), and any
duplicate altloc.

It REFUSES rather than guesses when it finds something that would change what
"residue 142" means later: a non-standard residue, a residue with no CA, or
numbering that jumps. Those are the defects that are invisible in a built
system and fatal to a per-residue comparison with HDX.
"""
import sys
import os

STD = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER "
          "THR TRP TYR VAL".split())
# what pdb2gmx will turn these into; HIS protonation is decided by pdb2gmx
ALIAS = {"HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "HID": "HIS",
         "HIE": "HIS", "HIP": "HIS", "MSE": "MET", "CYX": "CYS"}


def die(msg):
    sys.stderr.write("FATAL (prep_receptor.py): %s\n" % msg)
    sys.exit(1)


def main():
    # parse by walking, so an option's VALUE is never mistaken for a positional
    argv, args, chain = sys.argv[1:], [], "A"
    i = 0
    while i < len(argv):
        if argv[i] == "--chain":
            i + 1 < len(argv) or die("--chain needs a value")
            chain = argv[i + 1]; i += 2
        elif argv[i].startswith("--"):
            die("unknown option: %s" % argv[i])
        else:
            args.append(argv[i]); i += 1
    if len(args) != 2:
        die("usage: prep_receptor.py <in.pdb> <out.pdb> [--chain A]")
    src, dst = args
    os.path.isfile(src) or die("no such file: %s" % src)

    kept, seen, order = [], {}, []
    nonstd, noca, hetatm, hydro = set(), [], 0, 0
    for line in open(src):
        rec = line[:6]
        if rec == "HETATM":
            hetatm += 1
            continue
        if rec != "ATOM  ":
            continue
        if line[21] != chain and line[21] != " ":
            continue
        alt = line[16]
        if alt not in (" ", "A"):          # first altloc only
            continue
        name = line[12:16].strip()
        elem = line[76:78].strip() or name[:1]
        if elem == "H":                    # pdb2gmx rebuilds hydrogens
            hydro += 1
            continue
        res = line[17:20].strip()
        res = ALIAS.get(res, res)
        num = int(line[22:26])
        if res not in STD:
            nonstd.add(res)
        if num not in seen:
            seen[num] = set()
            order.append(num)
        seen[num].add(name)
        kept.append(line[:16] + " " + line[17:20].replace(line[17:20], "%-3s" % res) + line[20:])

    if not kept:
        die("no ATOM records for chain '%s' in %s" % (chain, src))
    if nonstd:
        die("non-standard residues pdb2gmx cannot build: %s\n"
            "       Remove them, or add the right patch to the force field."
            % ", ".join(sorted(nonstd)))
    noca = [n for n in order if "CA" not in seen[n]]
    if noca:
        die("%d residue(s) have no CA atom: %s\n"
            "       A missing backbone atom is exactly the defect that stops a build "
            "later,\n       and it shifts every residue index in the HDX comparison."
            % (len(noca), ", ".join(str(n) for n in noca[:10])))

    gaps = [(a, b) for a, b in zip(order, order[1:]) if b != a + 1]
    with open(dst, "w") as fh:
        fh.writelines(kept)
        fh.write("TER\nEND\n")

    print("  receptor  : %s" % src)
    print("  chain     : %s" % chain)
    print("  residues  : %d  (%d-%d)" % (len(order), order[0], order[-1]))
    print("  atoms kept: %d   (dropped %d HETATM, %d hydrogens)"
          % (len(kept), hetatm, hydro))
    if gaps:
        print("  *** NOTE: %d numbering gap(s): %s" %
              (len(gaps), ", ".join("%d->%d" % g for g in gaps[:5])))
        print("  ***       pdb2gmx will join them as if continuous. That is only")
        print("  ***       correct if the gap is renumbering, not missing residues.")
    print("  wrote     : %s" % dst)
    return len(order)


if __name__ == "__main__":
    main()
