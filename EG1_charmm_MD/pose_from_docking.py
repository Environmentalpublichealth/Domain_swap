#!/usr/bin/env python3
"""Extract one docked pose as the ligand coordinates the CHARMM build needs.

    pose_from_docking.py <docked.pdbqt|.pdb|.sdf> <out.pdb> [--pose 1]

AutoDock writes every pose into one .pdbqt as numbered MODELs, best first, with
AutoDock atom types and polar hydrogens. What 02_build_glycan.py wants is the
HEAVY ATOMS of a single pose as one residue: it re-derives every atom name from
the connectivity itself, so names and hydrogens here are irrelevant -- only the
heavy-atom positions matter.

--pose selects the MODEL (1 = best affinity, the default). The affinity of the
chosen pose is printed, because picking pose 3 and forgetting is a silent way to
simulate something you did not mean to.
"""
import sys
import os

def die(msg):
    sys.stderr.write("FATAL (pose_from_docking.py): %s\n" % msg)
    sys.exit(1)


def read_pdbqt(path, want):
    """Return (atom lines of the requested MODEL, its reported affinity)."""
    models, cur, aff, curaff = [], None, [], None
    for line in open(path):
        if line.startswith("MODEL"):
            cur, curaff = [], None
        elif line.startswith("ENDMDL"):
            if cur is not None:
                models.append(cur); aff.append(curaff)
            cur = None
        elif line.startswith("REMARK") and "VINA RESULT" in line:
            try:
                curaff = float(line.split()[3])
            except (IndexError, ValueError):
                pass
        elif line[:6] in ("ATOM  ", "HETATM") and cur is not None:
            cur.append(line)
        elif line[:6] in ("ATOM  ", "HETATM") and cur is None and not models:
            cur = [line]                      # single-pose file with no MODEL
    if cur:
        models.append(cur); aff.append(curaff)
    if not models:
        die("no atom records found in %s" % path)
    if want < 1 or want > len(models):
        die("--pose %d requested but the file has %d pose(s)" % (want, len(models)))
    return models[want - 1], aff[want - 1], len(models)


def read_sdf(path):
    lines = open(path).read().splitlines()
    try:
        natoms = int(lines[3][:3])
    except (IndexError, ValueError):
        die("cannot read the atom count from %s" % path)
    out = []
    for ln in lines[4:4 + natoms]:
        p = ln.split()
        out.append((float(p[0]), float(p[1]), float(p[2]), p[3]))
    return out


def main():
    # parse by walking: "--pose 3" must not leave "3" looking like a filename
    argv, args, want = sys.argv[1:], [], 1
    i = 0
    while i < len(argv):
        if argv[i] == "--pose":
            i + 1 < len(argv) or die("--pose needs a value")
            try:
                want = int(argv[i + 1])
            except ValueError:
                die("--pose needs a number, got '%s'" % argv[i + 1])
            i += 2
        elif argv[i].startswith("--"):
            die("unknown option: %s" % argv[i])
        else:
            args.append(argv[i]); i += 1
    if len(args) != 2:
        die("usage: pose_from_docking.py <docked.pdbqt|.pdb|.sdf> <out.pdb> [--pose N]")
    src, dst = args
    os.path.isfile(src) or die("no such file: %s" % src)

    heavy = []
    if src.endswith(".sdf"):
        for x, y, z, el in read_sdf(src):
            if el != "H":
                heavy.append((x, y, z, el))
        npose, aff = 1, None
    else:
        lines, aff, npose = read_pdbqt(src, want)
        for ln in lines:
            # In .pdbqt the AutoDock type sits in the last column; H and HD are
            # the hydrogens. Fall back to the element column for plain .pdb.
            tail = ln[77:].split()
            adtype = tail[-1] if tail else (ln[76:78].strip() or ln[12:16].strip()[:1])
            if adtype in ("H", "HD"):
                continue
            el = "C" if adtype in ("C", "A") else adtype[0]
            heavy.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54]), el))

    if not heavy:
        die("no heavy atoms left after dropping hydrogens")

    with open(dst, "w") as fh:
        for i, (x, y, z, el) in enumerate(heavy, 1):
            fh.write("HETATM%5d %-4s LIG A   1    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
                     % (i, ("%s%d" % (el, i))[:4], x, y, z, el))
        fh.write("END\n")

    print("  pose file : %s" % src)
    print("  pose      : %d of %d%s" % (want, npose,
          "  (affinity %.3f kcal/mol)" % aff if aff is not None else ""))
    print("  heavy atoms: %d" % len(heavy))
    print("  wrote     : %s" % dst)


if __name__ == "__main__":
    main()
