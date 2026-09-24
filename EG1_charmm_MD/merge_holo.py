#!/usr/bin/env python3
"""Merge the built glycan into the protein's coordinates and topology.

    merge_holo.py <processed.gro> <glycan.pdb> <glycan.itp> <system.top> <out.gro>

Called by prep_system.sh for holo systems only. Two jobs:

COORDINATES  Append the glycan's 105 atoms to the protein .gro. Order matters
             absolutely: GROMACS matches coordinates to [ molecules ] by
             POSITION, not by name, so the glycan must sit directly after the
             protein and before the solvent that gets added next.

TOPOLOGY     Insert `#include "glycan.itp"` after the protein moleculetype and
             before [ system ] (a moleculetype has to be defined before it is
             used), then add the glycan to [ molecules ] immediately after the
             protein line -- again so the topology order matches the coordinate
             order.

The .top is rewritten in place after a .bak copy, and the script refuses to run
twice on the same file.
"""
import sys
import os
import re


def die(msg):
    sys.stderr.write("FATAL (merge_holo.py): %s\n" % msg)
    sys.exit(1)


def read_gro(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    if len(lines) < 3:
        die("%s is too short to be a .gro" % path)
    try:
        n = int(lines[1].strip())
    except ValueError:
        die("%s line 2 is not an atom count" % path)
    if len(lines) < n + 3:
        die("%s claims %d atoms but has %d lines" % (path, n, len(lines) - 3))
    return lines[0], n, lines[2:2 + n], lines[2 + n]


def read_pdb(path):
    out = []
    with open(path) as fh:
        for line in fh:
            if line[:6] in ('ATOM  ', 'HETATM'):
                out.append({
                    'resid': int(line[22:26]),
                    # cols 18-21: CHARMM residue names are four characters.
                    'resname': line[17:21].strip(),
                    'name': line[12:16].strip(),
                    # PDB is A, .gro is nm
                    'xyz': (float(line[30:38]) / 10.0,
                            float(line[38:46]) / 10.0,
                            float(line[46:54]) / 10.0),
                })
    return out


def molecule_name(itp):
    """The [ moleculetype ] name declared in the glycan .itp."""
    section = None
    with open(itp) as fh:
        for raw in fh:
            line = raw.split(';', 1)[0].strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('[') and line.endswith(']'):
                section = line[1:-1].strip().lower()
                continue
            if section == 'moleculetype':
                return line.split()[0]
    die("no [ moleculetype ] name in %s" % itp)


def main():
    if len(sys.argv) != 6:
        die("usage: merge_holo.py <processed.gro> <glycan.pdb> <glycan.itp> "
            "<system.top> <out.gro>")
    pro_gro, gly_pdb, gly_itp, top, out_gro = sys.argv[1:]

    for f in (pro_gro, gly_pdb, gly_itp, top):
        if not os.path.exists(f):
            die("%s not found" % f)

    molname = molecule_name(gly_itp)
    title, npro, pro_atoms, box = read_gro(pro_gro)
    gly = read_pdb(gly_pdb)
    if not gly:
        die("%s contains no atoms" % gly_pdb)

    # ---- coordinates -------------------------------------------------------
    # Residue numbers continue from the protein so the merged .gro reads
    # sensibly in a viewer. GROMACS itself does not use them.
    last_resid = 0
    if pro_atoms:
        try:
            last_resid = int(pro_atoms[-1][0:5])
        except ValueError:
            last_resid = 0
    resid_map, nxt = {}, last_resid
    for a in gly:
        if a['resid'] not in resid_map:
            nxt += 1
            resid_map[a['resid']] = nxt

    lines = list(pro_atoms)
    for i, a in enumerate(gly):
        idx = npro + i + 1
        lines.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f"
                     % (resid_map[a['resid']] % 100000, a['resname'][:5],
                        a['name'][:5].rjust(5), idx % 100000,
                        a['xyz'][0], a['xyz'][1], a['xyz'][2]))

    with open(out_gro, 'w') as fh:
        fh.write("%s + %s (%d atoms)\n" % (title.strip(), molname, len(gly)))
        fh.write("%5d\n" % (npro + len(gly)))
        fh.write("\n".join(lines) + "\n")
        fh.write(box + "\n")

    # ---- topology ----------------------------------------------------------
    with open(top) as fh:
        toplines = fh.read().splitlines()

    inc = '#include "glycan.itp"'
    if any(l.strip() == inc for l in toplines):
        die("%s already includes glycan.itp. prep_system.sh should have started "
            "from a clean pdb2gmx topology; delete the build directory and "
            "re-run rather than merging twice." % top)

    sys_idx = next((i for i, l in enumerate(toplines)
                    if re.match(r'^\s*\[\s*system\s*\]', l)), None)
    if sys_idx is None:
        die("no [ system ] section in %s" % top)

    mol_idx = next((i for i, l in enumerate(toplines)
                    if re.match(r'^\s*\[\s*molecules\s*\]', l)), None)
    if mol_idx is None:
        die("no [ molecules ] section in %s" % top)

    # First real entry under [ molecules ] is the protein; the glycan goes
    # directly after it, matching the coordinate order written above.
    ins_mol = None
    for i in range(mol_idx + 1, len(toplines)):
        s = toplines[i].split(';', 1)[0].strip()
        if s:
            ins_mol = i + 1
            break
    if ins_mol is None:
        die("[ molecules ] in %s has no entries" % top)

    width = max(len(molname), 20)
    new = (toplines[:sys_idx]
           + ["; ---- xylohexaose ligand, built by 02_build_glycan.py ----------------",
              "; Defined here, AFTER the protein moleculetype and BEFORE [ system ],",
              "; because a moleculetype must exist before [ molecules ] refers to it.",
              inc,
              ""]
           + toplines[sys_idx:ins_mol]
           + ["%-*s 1" % (width, molname)]
           + toplines[ins_mol:])

    os.replace(top, top + '.bak')
    with open(top, 'w') as fh:
        fh.write("\n".join(new) + "\n")

    print("  merged glycan '%s': %d protein + %d ligand = %d atoms"
          % (molname, npro, len(gly), npro + len(gly)))
    print("  topology: added %s before [ system ], and '%s 1' as the second"
          % (inc, molname))
    print("            entry in [ molecules ]  (backup: %s)" % os.path.basename(top + '.bak'))


if __name__ == '__main__':
    main()
