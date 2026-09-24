#!/usr/bin/env python3
"""
Turn the CHARMM-GUI xylohexaose topology into a ready-to-use ligand for each
holo system, in that system's own AlphaFold3 binding pose.

    python3 02_build_glycan.py            # all holo systems
    python3 02_build_glycan.py q4_holo    # just one

WHAT THIS EXISTS FOR
--------------------
The CHARMM36 GROMACS port ships beta-D-xylose as a MONOSACCHARIDE (carb.rtp
[ BXYL ]) and nothing that joins two of them: carb.n.tdb and carb.c.tdb contain
only [ None ], and there is no 1->4 glycosidic patch anywhere in the port. So
pdb2gmx can build the protein but physically cannot build xylohexaose. The
linkage parameters have to come from outside, and CHARMM-GUI's Glycan Reader &
Modeler is where they come from.

You only run CHARMM-GUI ONCE. The molecule is the same in all three holo
systems -- only its pose differs -- so one topology is reused for all three, and
this script is what puts it into each pose.

WHAT IT NEEDS FROM YOU
----------------------
One .itp file under ligands/<LIGAND>_charmm/, holding a moleculetype with exactly
105 atoms (C30 H50 O25 -- six xyloses minus five bridging waters). Coordinates
are OPTIONAL: if a matching structure file is found the hydroxyl hydrogens
inherit their torsions from it, and if not they are built anti, which costs
nothing because the restrained NVT stage leaves hydrogens free to relax.

HOW THE POSE TRANSFER WORKS
---------------------------
AlphaFold3 writes the ligand as one residue called LIG with invented names
(C1..C30, O1..O25) and no hydrogens. CHARMM needs six residues with carbohydrate
names. Rather than guess a mapping, both molecules are labelled CANONICALLY from
their own connectivity, which for a linear xylo-oligosaccharide is unambiguous:

  * a pyranose ring is the unique 6-cycle carrying one oxygen
  * within a ring, C1 is the anomeric carbon -- the only ring carbon bonded to
    two oxygens -- and O5 is the ring oxygen next to it; walking the ring away
    from O5 gives C2, C3, C4, C5
  * the reducing-end residue is the one whose anomeric oxygen is terminal;
    residue i's O4 bridges to residue i+1's C1, which fixes the chain direction
  * CHARMM's own convention is then followed exactly: the bridging oxygen is
    named O4 of the ACCEPTOR residue, residue 1 keeps O1/HO1, residue 6 keeps
    O4/HO4, and the five donors lose theirs -- which is precisely why the whole
    molecule comes to 105 atoms and not 120

Because that labelling is derived from connectivity on BOTH sides, the map
between them is exact rather than fitted, and it is verified afterwards: every
name the .itp declares must be produced, every bond in the .itp must exist in
the AlphaFold3 pose, and the handedness at all 24 stereocentres must agree
across all six residues.

Hydrogens are then placed from local geometry -- not copied by superposition --
so a ring pucker that differs from CHARMM-GUI's idealised build cannot smuggle a
distorted hydrogen in.
"""
import os
import sys
import math
import glob

HERE = os.path.dirname(os.path.abspath(__file__))

BOND_CUT = 1.75          # heavy-atom covalent cutoff, A
CH_BOND = 1.09           # C-H, A
OH_BOND = 0.96           # O-H, A
TETRA = math.radians(109.5)
COH_ANGLE = math.radians(108.0)

EXPECT_ATOMS = 105       # C30 H50 O25
EXPECT_C, EXPECT_H, EXPECT_O = 30, 50, 25
N_RES = 6


class Problem(Exception):
    pass


# ------------------------------------------------------------------ vector ops
def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    n = norm(a)
    if n < 1e-9:
        raise Problem("degenerate vector")
    return scale(a, 1.0 / n)


def dist(a, b):
    return norm(sub(a, b))


def dihedral(p0, p1, p2, p3):
    """Signed dihedral p0-p1-p2-p3 in radians."""
    b0 = sub(p0, p1)
    b1 = unit(sub(p2, p1))
    b2 = sub(p3, p2)
    v = sub(b0, scale(b1, dot(b0, b1)))
    w = sub(b2, scale(b1, dot(b2, b1)))
    return math.atan2(dot(cross(b1, v), w), dot(v, w))


def place_from_dihedral(a, b, c, bond, angle, tors):
    """Position of atom d given a-b-c, |c-d|, angle b-c-d and dihedral a-b-c-d."""
    bc = unit(sub(c, b))
    ab = sub(b, a)
    n = cross(ab, bc)
    if norm(n) < 1e-6:                       # a, b, c collinear
        seed = (1.0, 0.0, 0.0) if abs(bc[0]) < 0.9 else (0.0, 1.0, 0.0)
        n = cross(bc, seed)
    n = unit(n)
    m = cross(n, bc)
    d2 = (-bond * math.cos(angle),
          bond * math.sin(angle) * math.cos(tors),
          bond * math.sin(angle) * math.sin(tors))
    return (c[0] + d2[0] * bc[0] + d2[1] * m[0] + d2[2] * n[0],
            c[1] + d2[0] * bc[1] + d2[1] * m[1] + d2[2] * n[1],
            c[2] + d2[0] * bc[2] + d2[1] * m[2] + d2[2] * n[2])


# ---------------------------------------------------------------- file parsing
def parse_structure(path):
    """Read a .pdb or .gro. Returns [{'name','resname','resid','xyz'}] in A."""
    atoms = []
    if path.endswith('.gro'):
        with open(path) as fh:
            lines = fh.read().splitlines()
        try:
            n = int(lines[1].strip())
        except (IndexError, ValueError):
            return []
        for line in lines[2:2 + n]:
            if len(line) < 44:
                continue
            atoms.append({
                'resid': int(line[0:5]),
                'resname': line[5:10].strip(),
                'name': line[10:15].strip(),
                # .gro is nm; everything in this script works in A
                'xyz': (float(line[20:28]) * 10.0,
                        float(line[28:36]) * 10.0,
                        float(line[36:44]) * 10.0),
            })
    else:
        with open(path) as fh:
            for line in fh:
                if line[:6] not in ('ATOM  ', 'HETATM'):
                    continue
                if line[16] not in (' ', 'A'):
                    continue
                atoms.append({
                    'resid': int(line[22:26]),
                    # cols 18-21, not 18-20: CHARMM residue names run to four
                    # characters (BXYL), and CHARMM-GUI writes them that way.
                    # Column 21 is blank in a standard 3-character PDB, so
                    # .strip() reads both correctly.
                    'resname': line[17:21].strip(),
                    'name': line[12:16].strip(),
                    'xyz': (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                })
    return atoms


def element_of(name, atomtype=''):
    """Element from a CHARMM atom name. H first: HO1, H51 etc. all start H."""
    n = name.strip().upper()
    if n.startswith('H'):
        return 'H'
    if n.startswith('C'):
        return 'C'
    if n.startswith('O'):
        return 'O'
    t = atomtype.strip().upper()
    for e in ('H', 'C', 'O'):
        if t.startswith(e) or t.startswith(e + 'C'):
            return e
    raise Problem("cannot tell the element of atom '%s' (type '%s')" % (name, atomtype))


def parse_itp_moleculetypes(path):
    """Split an .itp into moleculetypes with their atoms and bonds."""
    mols, cur, section = [], None, None
    with open(path) as fh:
        for raw in fh:
            line = raw.split(';', 1)[0].rstrip()
            if not line.strip():
                continue
            if line.lstrip().startswith('#'):
                continue
            s = line.strip()
            if s.startswith('[') and s.endswith(']'):
                section = s[1:-1].strip().lower()
                if section == 'moleculetype':
                    cur = {'name': None, 'atoms': [], 'bonds': [], 'path': path}
                    mols.append(cur)
                continue
            if cur is None:
                continue
            f = s.split()
            if section == 'moleculetype' and cur['name'] is None:
                cur['name'] = f[0]
            elif section == 'atoms' and len(f) >= 5 and f[0].isdigit():
                cur['atoms'].append({
                    'nr': int(f[0]), 'type': f[1], 'resnr': int(f[2]),
                    'resname': f[3], 'name': f[4],
                    'charge': float(f[6]) if len(f) > 6 else 0.0,
                })
            elif section == 'bonds' and len(f) >= 2 and f[0].isdigit() and f[1].isdigit():
                cur['bonds'].append((int(f[0]), int(f[1])))
    return mols


def find_glycan_itp(root):
    """Locate the 105-atom moleculetype under charmm_gui/. Reports what it saw."""
    itps = sorted(glob.glob(os.path.join(root, '**', '*.itp'), recursive=True))
    itps += sorted(glob.glob(os.path.join(root, '**', '*.top'), recursive=True))
    if not itps:
        raise Problem(
            "no .itp or .top files under %s.\n"
            "       Unpack the CHARMM-GUI download there:\n"
            "         see ligands/README.md for how to generate it" % root)

    seen, hits = [], []
    for path in itps:
        try:
            for mol in parse_itp_moleculetypes(path):
                if not mol['atoms']:
                    continue
                rel = os.path.relpath(path, root)
                seen.append((rel, mol['name'], len(mol['atoms'])))
                if len(mol['atoms']) == EXPECT_ATOMS:
                    hits.append(mol)
        except (ValueError, IndexError):
            continue

    if not hits:
        msg = ["no moleculetype with %d atoms (C30 H50 O25 = xylohexaose) found "
               "under %s." % (EXPECT_ATOMS, root),
               "       Moleculetypes that ARE there:"]
        for rel, nm, na in seen[:40]:
            msg.append("         %-46s %-12s %4d atoms" % (rel, nm, na))
        if not seen:
            msg.append("         (none -- the .itp files parsed to nothing)")
        msg.append("")
        msg.append("       A 105-atom entry is the fingerprint of the glycan. If the")
        msg.append("       closest thing above is 120 atoms, the six xyloses were built")
        msg.append("       UNLINKED; if it is ~20, only one xylose was built.")
        raise Problem("\n".join(msg))

    if len(hits) > 1:
        # The full CHARMM-GUI GROMACS output writes the same moleculetype into
        # topol.top AND into toppar/<name>.itp, so finding it twice is normal and
        # is not an ambiguity. Only genuinely DIFFERENT molecules are a problem.
        def signature(m):
            # Residue name, atom type and charge all belong here, not just the
            # atom names: two files can agree on every atom NAME and still be
            # different molecules (BXYL vs AXYL is the same 105 names with the
            # opposite anomeric configuration, and would be silently accepted by
            # a name-only signature).
            return (m['name'],
                    tuple((a['resnr'], a['resname'], a['name'], a['type'],
                           round(a['charge'], 4)) for a in m['atoms']),
                    tuple(sorted(m['bonds'])))
        sigs = {signature(m) for m in hits}
        if len(sigs) == 1:
            where = ', '.join(os.path.relpath(m['path'], root) for m in hits)
            print("  note                     the same molecule appears in %d files"
                  % len(hits))
            print("                           (%s);" % where)
            print("                           they are identical, using the first.")
        else:
            names = ', '.join("%s in %s" % (m['name'], os.path.relpath(m['path'], root))
                              for m in hits)
            raise Problem("found %d DIFFERENT 105-atom moleculetypes (%s).\n"
                          "       Leave only one CHARMM-GUI job unpacked under %s."
                          % (len(hits), names, root))
    return hits[0]


def find_template_coords(root, mol):
    """Optional: coordinates for the glycan, used only for hydroxyl torsions."""
    want = sorted((a['resnr'], a['name']) for a in mol['atoms'])
    cands = []
    for pat in ('*.pdb', '*.gro'):
        cands += glob.glob(os.path.join(root, '**', pat), recursive=True)
    best = None
    for path in sorted(cands):
        try:
            atoms = parse_structure(path)
        except (ValueError, IndexError):
            continue
        if not atoms:
            continue
        by_name = {}
        for a in atoms:
            by_name.setdefault((a['resid'], a['name']), a['xyz'])
        resids = sorted({r for r, _ in by_name})
        # Try every offset that could align the file's residue numbering with
        # the itp's (CHARMM-GUI often starts glycan residues at 1, but not
        # always).
        for off in {r - want[0][0] for r in resids}:
            shifted = [((r + off), n) for r, n in want]
            if all(k in by_name for k in shifted):
                coords = {want[i]: by_name[shifted[i]] for i in range(len(want))}
                if best is None or len(atoms) < best[1]:
                    best = (coords, len(atoms), path)
    return best


# ------------------------------------------------------- canonical sugar labels
def bonds_by_distance(coords):
    adj = {i: set() for i in range(len(coords))}
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            if dist(coords[i], coords[j]) <= BOND_CUT:
                adj[i].add(j)
                adj[j].add(i)
    return adj


def find_rings(adj, size=6):
    found = set()
    for start in adj:
        stack = [(start, [start])]
        while stack:
            node, path = stack.pop()
            if len(path) == size:
                if start in adj[node]:
                    found.add(frozenset(path))
                continue
            for nxt in adj[node]:
                if nxt > start and nxt not in path:
                    stack.append((nxt, path + [nxt]))
    return found


def canonical_labels(coords, elements):
    """Label a heavy-atom xylo-oligosaccharide the way CHARMM names it.

    Returns {(resid, atomname): atom_index}, resid counting 1..6 from the
    reducing end. The bridging oxygen between residue i and i+1 is named O4 of
    residue i -- CHARMM's convention -- and appears once, not twice.
    """
    adj = bonds_by_distance(coords)
    rings = find_rings(adj, 6)
    if len(rings) != N_RES:
        raise Problem("found %d six-membered rings, expected %d" % (len(rings), N_RES))

    ring_info = []
    for ring in rings:
        ring_o = [i for i in ring if elements[i] == 'O']
        if len(ring_o) != 1:
            raise Problem("a ring has %d oxygens; a pyranose has one" % len(ring_o))
        o5 = ring_o[0]
        ring_c = [i for i in ring if elements[i] == 'C']
        anomeric = [c for c in ring_c
                    if o5 in adj[c] and sum(1 for n in adj[c] if elements[n] == 'O') == 2]
        if len(anomeric) != 1:
            raise Problem("could not identify a unique anomeric carbon in a ring "
                          "(%d candidates)" % len(anomeric))
        c1 = anomeric[0]
        # Walk the ring away from O5: C1 -> C2 -> C3 -> C4 -> C5 -> O5.
        chain, prev, cur = [c1], o5, c1
        while len(chain) < 5:
            nxt = [n for n in adj[cur] if n in ring and n != prev]
            if len(nxt) != 1:
                raise Problem("ring walk is ambiguous at atom %d" % cur)
            prev, cur = cur, nxt[0]
            chain.append(cur)
        if o5 not in adj[chain[4]]:
            raise Problem("ring walk did not close back onto the ring oxygen")
        c1, c2, c3, c4, c5 = chain

        def exo_o(c, inring):
            got = [n for n in adj[c] if elements[n] == 'O' and n not in inring]
            if len(got) != 1:
                raise Problem("carbon %d carries %d exocyclic oxygens, expected 1"
                              % (c, len(got)))
            return got[0]

        ring_info.append({
            'ring': set(ring), 'C1': c1, 'C2': c2, 'C3': c3, 'C4': c4, 'C5': c5,
            'O5': o5,
            'Oano': exo_o(c1, ring), 'O2': exo_o(c2, ring),
            'O3': exo_o(c3, ring), 'O4': exo_o(c4, ring),
        })

    # Reducing end: its anomeric oxygen is terminal (bonded only to C1).
    reducing = [k for k, r in enumerate(ring_info) if len(adj[r['Oano']]) == 1]
    if len(reducing) != 1:
        raise Problem("expected exactly one free anomeric oxygen (the reducing "
                      "end); found %d. This is not a linear oligosaccharide with "
                      "a free reducing end." % len(reducing))

    # Order residues by following O4 -> next residue's C1.
    c1_owner = {r['C1']: k for k, r in enumerate(ring_info)}
    order, cur = [reducing[0]], reducing[0]
    while len(order) < N_RES:
        o4 = ring_info[cur]['O4']
        nxt = [c1_owner[n] for n in adj[o4] if n in c1_owner and c1_owner[n] != cur]
        if len(nxt) != 1:
            raise Problem("residue %d's O4 does not bridge to exactly one next "
                          "residue (found %d)" % (len(order), len(nxt)))
        cur = nxt[0]
        if cur in order:
            raise Problem("the residue chain loops back on itself")
        order.append(cur)

    labels = {}
    for pos, k in enumerate(order, start=1):
        r = ring_info[k]
        for nm in ('C1', 'C2', 'C3', 'C4', 'C5', 'O5', 'O2', 'O3', 'O4'):
            labels[(pos, nm)] = r[nm]
        if pos == 1:
            labels[(pos, 'O1')] = r['Oano']
    return labels, adj, ring_info, order


def chirality(coords, centre, a, b, c):
    """Handedness at a tetrahedral centre from three of its substituents."""
    v1, v2, v3 = (sub(coords[a], coords[centre]),
                  sub(coords[b], coords[centre]),
                  sub(coords[c], coords[centre]))
    t = dot(v1, cross(v2, v3))
    return 1 if t > 0 else -1


# --------------------------------------------------------------- H placement
def place_hydrogens(coords, labels, adj, elements, torsions):
    """Build every hydrogen the CHARMM topology needs, from local geometry."""
    h = {}

    def neighbours(idx):
        return sorted(adj[idx])

    for res in range(1, N_RES + 1):
        # --- methine hydrogens: three heavy neighbours fix the fourth vertex
        for cname, hname in (('C1', 'H1'), ('C2', 'H2'), ('C3', 'H3'), ('C4', 'H4')):
            c = labels[(res, cname)]
            nb = neighbours(c)
            if len(nb) != 3:
                raise Problem("res %d %s has %d heavy neighbours, expected 3"
                              % (res, cname, len(nb)))
            direction = (0.0, 0.0, 0.0)
            for n in nb:
                direction = add(direction, unit(sub(coords[n], coords[c])))
            h[(res, hname)] = add(coords[c], scale(unit(scale(direction, -1.0)), CH_BOND))

        # --- the C5 methylene: two heavy neighbours, two hydrogens
        c5 = labels[(res, 'C5')]
        nb = neighbours(c5)
        if len(nb) != 2:
            raise Problem("res %d C5 has %d heavy neighbours, expected 2" % (res, len(nb)))
        u1 = unit(sub(coords[nb[0]], coords[c5]))
        u2 = unit(sub(coords[nb[1]], coords[c5]))
        bisec = unit(scale(add(u1, u2), -1.0))
        perp = unit(cross(u1, u2))
        half = TETRA / 2.0
        for hname, sgn in (('H51', 1.0), ('H52', -1.0)):
            d = add(scale(bisec, math.cos(half)), scale(perp, sgn * math.sin(half)))
            h[(res, hname)] = add(coords[c5], scale(unit(d), CH_BOND))

        # --- hydroxyl hydrogens. Only the free ones exist: residue 1 keeps
        # --- HO1, residue N keeps HO4, and every residue keeps HO2/HO3.
        todo = [('O2', 'HO2', 'C2', 'C1'), ('O3', 'HO3', 'C3', 'C2')]
        if res == 1:
            todo.append(('O1', 'HO1', 'C1', 'O5'))
        if res == N_RES:
            todo.append(('O4', 'HO4', 'C4', 'C3'))
        for oname, hname, cname, refname in todo:
            o = labels[(res, oname)]
            c = labels[(res, cname)]
            ref = labels[(res, refname)]
            if len(adj[o]) != 1:
                raise Problem("res %d %s should be a free hydroxyl but has %d "
                              "heavy neighbours" % (res, oname, len(adj[o])))
            tors = torsions.get((res, hname), math.pi)   # anti by default
            h[(res, hname)] = place_from_dihedral(coords[ref], coords[c], coords[o],
                                                  OH_BOND, COH_ANGLE, tors)
    return h


def relax_hydroxyls(positions, bonded, report):
    """Rotate clashing hydroxyl hydrogens to the least-clashing torsion.

    Hydroxyl H positions come from the CHARMM-GUI reference conformer, or anti
    when there is none. Neither knows anything about THIS pose, so two -OH
    hydrogens can land on top of each other -- which the guard below then
    (correctly) refuses. A hydroxyl torsion is free to rotate at ~0 cost, so the
    honest fix is to turn it rather than to lower the threshold: scan the torsion
    and keep the angle that maximises the closest non-bonded contact.

    Only H atoms named HO* move, and only their torsion changes: every bond
    length and angle is untouched, so the topology still matches the geometry.
    """
    moved = []
    for key in sorted(positions):
        resnr, name = key
        if not name.startswith('HO'):
            continue
        oname = 'O' + name[2:]                      # HO2 -> O2
        cname = 'C' + name[2:]                      # HO2 -> C2
        ref = {'O1': 'O5', 'O2': 'C1', 'O3': 'C2', 'O4': 'C3'}.get(oname)
        for k in ((resnr, oname), (resnr, cname), (resnr, ref)):
            if k not in positions:
                ref = None
                break
        if ref is None:
            continue
        others = [k for k in positions
                  if k != key and (key, k) not in bonded and k != (resnr, oname)]

        def closest_for(pos):
            return min(dist(pos, positions[k]) for k in others)

        cur = closest_for(positions[key])
        if cur >= 1.20:
            continue                                # already fine, leave it alone
        best, best_d = positions[key], cur
        for step in range(36):                      # 10 degree scan
            cand = place_from_dihedral(positions[(resnr, ref)],
                                       positions[(resnr, cname)],
                                       positions[(resnr, oname)],
                                       OH_BOND, COH_ANGLE,
                                       math.radians(step * 10.0))
            d = closest_for(cand)
            if d > best_d:
                best, best_d = cand, d
        if best_d > cur:
            positions[key] = best
            moved.append((name, resnr, cur, best_d))
    if moved:
        report.append("hydroxyl relaxation      %d rotated (worst contact %.2f -> %.2f A)"
                      % (len(moved), min(m[2] for m in moved), min(m[3] for m in moved)))
    return positions


def template_torsions(tcoords, tlabels):
    """Read hydroxyl H-O-C-X torsions out of the CHARMM-GUI coordinates."""
    out = {}
    for res in range(1, N_RES + 1):
        todo = [('O2', 'HO2', 'C2', 'C1'), ('O3', 'HO3', 'C3', 'C2')]
        if res == 1:
            todo.append(('O1', 'HO1', 'C1', 'O5'))
        if res == N_RES:
            todo.append(('O4', 'HO4', 'C4', 'C3'))
        for oname, hname, cname, refname in todo:
            try:
                out[(res, hname)] = dihedral(tcoords[(res, refname)],
                                             tcoords[(res, cname)],
                                             tcoords[(res, oname)],
                                             tcoords[(res, hname)])
            except (KeyError, Problem):
                continue
    return out


# --------------------------------------------------------------------- output
def write_pdb(path, mol, positions, resnames, resids):
    with open(path, 'w') as fh:
        fh.write("REMARK   xylohexaose, CHARMM36 carbohydrate topology from "
                 "CHARMM-GUI,\n")
        fh.write("REMARK   placed in the docked pose by 02_build_glycan.py\n")
        for i, a in enumerate(mol['atoms'], start=1):
            key = (a['resnr'], a['name'])
            x, y, z = positions[key]
            elem = element_of(a['name'], a['type'])
            fh.write("ATOM  %5d %-4s %-4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f"
                     "          %2s\n"
                     % (i, (' ' + a['name'])[:4] if len(a['name']) < 4 else a['name'],
                        resnames[a['resnr']][:4], 'L', resids[a['resnr']],
                        x, y, z, 1.00, 0.00, elem))
        fh.write("TER\nEND\n")


def write_itp(path, src_path, molname, natoms):
    """Copy the CHARMM-GUI .itp, swapping its POSRES block for our own.

    CHARMM-GUI writes its own '#ifdef POSRES -> posre_<x>.itp'. Left in place it
    would collide with the protein's posre.itp, and -DPOSRES would then restrain
    either the wrong atoms or none. It is replaced with a block pointing at
    posre_lig.itp, whose indices are ligand-local.
    """
    keep, skipping = [], False
    with open(src_path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('#ifdef') and 'POSRES' in s.upper():
                skipping = True
                continue
            if skipping:
                if s.startswith('#endif'):
                    skipping = False
                continue
            keep.append(line)
    with open(path, 'w') as fh:
        fh.write("; xylohexaose -- CHARMM36 carbohydrate parameters.\n")
        fh.write("; Copied verbatim from the CHARMM-GUI output (%s), except that\n"
                 % os.path.basename(src_path))
        fh.write("; CHARMM-GUI's own POSRES block was removed and replaced below:\n")
        fh.write("; its posre file would have collided with the protein's.\n;\n")
        fh.write("; moleculetype %s, %d atoms.\n;\n" % (molname, natoms))
        fh.writelines(keep)
        fh.write("\n")
        fh.write("; Ligand heavy atoms, restrained during NVT only. npt.mdp uses\n")
        fh.write("; -DPOSRES_BB, which this file deliberately does NOT define, so the\n")
        fh.write("; ligand is free from the NPT stage onwards.\n")
        fh.write("#ifdef POSRES\n#include \"posre_lig.itp\"\n#endif\n")


def write_posre(path, mol, fc=1000):
    with open(path, 'w') as fh:
        fh.write("; Position restraints on the xylohexaose HEAVY atoms.\n")
        fh.write("; Indices are LIGAND-LOCAL: this file is #included inside the\n")
        fh.write("; ligand moleculetype, not the protein's.\n")
        fh.write("[ position_restraints ]\n")
        fh.write("; atom  type      fx      fy      fz\n")
        n = 0
        for a in mol['atoms']:
            if element_of(a['name'], a['type']) == 'H':
                continue
            fh.write("%6d     1  %6d  %6d  %6d\n" % (a['nr'], fc, fc, fc))
            n += 1
    return n


# ----------------------------------------------------------------------- main
def read_conf(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.split('#', 1)[0].strip()
            if '=' in line:
                k, v = line.split('=', 1)
                out[k.strip()] = v.strip().strip('"').strip("'")
    return out


def build_one(name, mol, tcoords, report):
    # A repaired pose wins over AlphaFold3's own, when 03_fix_ligand_stereo.py
    # has written one. Keeping both files means the substitution is always
    # visible on disk rather than buried in a log.
    lig_pdb = os.path.join(HERE, 'build', name, 'ligand_pose.pdb')
    if os.path.exists(lig_pdb):
        report.append("pose source              ligand_pose.pdb (repaired -- see "
                      "03_fix_ligand_stereo.py)")
    else:
        lig_pdb = os.path.join(HERE, 'build', name, 'ligand_af3.pdb')
        report.append("pose source              ligand_af3.pdb (docked pose, unmodified heavy atoms)")
    if not os.path.exists(lig_pdb):
        raise Problem("build/%s/ligand_af3.pdb not found. Run "
                      "`python3 00_prep_inputs.py %s` first." % (name, name))
    atoms = parse_structure(lig_pdb)
    coords = [a['xyz'] for a in atoms]
    elements = [element_of(a['name']) for a in atoms]
    if len(atoms) != EXPECT_C + EXPECT_O:
        raise Problem("ligand_af3.pdb has %d atoms, expected %d heavy atoms"
                      % (len(atoms), EXPECT_C + EXPECT_O))

    labels, adj, _ring_info, _order = canonical_labels(coords, elements)

    # Every name the topology declares must have been produced by the labeller.
    itp_names = {(a['resnr'] - min(x['resnr'] for x in mol['atoms']) + 1, a['name'])
                 for a in mol['atoms']}
    heavy_itp = {k for k in itp_names if not k[1].upper().startswith('H')}
    missing = sorted(heavy_itp - set(labels))
    extra = sorted(set(labels) - heavy_itp)
    if missing or extra:
        raise Problem("the canonical labelling does not match the topology.\n"
                      "       in the .itp but not built: %s\n"
                      "       built but not in the .itp: %s\n"
                      "       If the .itp numbers its residues from the "
                      "NON-reducing end, that is the cause -- tell me and I will "
                      "flip the direction." % (missing or 'none', extra or 'none'))

    # Chirality must be identical in all six residues: it is a homopolymer, so a
    # single inverted centre means one residue is not beta-D-xylose.
    signs = {}
    for res in range(1, N_RES + 1):
        signs[res] = (
            chirality(coords, labels[(res, 'C1')], labels[(res, 'O5')],
                      labels[(res, 'C2')], labels[(res, 'O1')] if res == 1
                      else labels[(res - 1, 'O4')] if res > 1 else labels[(res, 'O5')]),
            chirality(coords, labels[(res, 'C2')], labels[(res, 'C1')],
                      labels[(res, 'C3')], labels[(res, 'O2')]),
            chirality(coords, labels[(res, 'C3')], labels[(res, 'C2')],
                      labels[(res, 'C4')], labels[(res, 'O3')]),
            chirality(coords, labels[(res, 'C4')], labels[(res, 'C3')],
                      labels[(res, 'C5')], labels[(res, 'O4')]),
        )
    distinct = set(signs.values())
    if len(distinct) != 1:
        detail = '; '.join("res%d %s" % (r, s) for r, s in sorted(signs.items()))
        raise Problem("the six residues do not share one handedness (%s).\n"
                      "       beta-1,4-xylohexaose is a homopolymer, so at least one "
                      "unit here is\n       a different sugar -- AlphaFold3 does not "
                      "enforce ligand stereochemistry.\n"
                      "       Simulating it would put xylose parameters on an "
                      "arabinose geometry\n       for the whole run, and no log would "
                      "ever mention it.\n"
                      "       Diagnose and repair:  python3 03_fix_ligand_stereo.py"
                      % detail)
    report.append("stereocentres            24 centres, all six residues identical "
                  "(C1,C2,C3,C4 = %s)" % ','.join('%+d' % s for s in distinct.pop()))

    torsions = template_torsions(tcoords, labels) if tcoords else {}
    hpos = place_hydrogens(coords, labels, adj, elements, torsions)

    positions = {}
    base = min(a['resnr'] for a in mol['atoms']) - 1
    for a in mol['atoms']:
        key = (a['resnr'] - base, a['name'])
        if key in labels:
            positions[(a['resnr'], a['name'])] = coords[labels[key]]
        elif key in hpos:
            positions[(a['resnr'], a['name'])] = hpos[key]
        else:
            raise Problem("no position produced for %s of residue %d"
                          % (a['name'], a['resnr']))

    # Every bond the topology declares must be a real bond in this pose.
    idx = {a['nr']: (a['resnr'], a['name']) for a in mol['atoms']}
    worst, worst_pair = 0.0, None
    for i, j in mol['bonds']:
        d = dist(positions[idx[i]], positions[idx[j]])
        ei = element_of(idx[i][1])
        ej = element_of(idx[j][1])
        lim = 1.30 if 'H' in (ei, ej) else 1.90
        if d > lim:
            raise Problem("bond %s%d-%s%d in the topology is %.2f A in this "
                          "pose -- not a bond. The mapping is wrong."
                          % (idx[i][1], idx[i][0], idx[j][1], idx[j][0], d))
        if d > worst:
            worst, worst_pair = d, (idx[i], idx[j])
    report.append("topology bonds           all %d present in the pose (longest "
                  "%.2f A, %s%d-%s%d)"
                  % (len(mol['bonds']), worst, worst_pair[0][1], worst_pair[0][0],
                     worst_pair[1][1], worst_pair[1][0]))

    # Placed hydrogens must not be sitting on top of something else.
    bonded = set()
    for i, j in mol['bonds']:
        bonded.add((idx[i], idx[j]))
        bonded.add((idx[j], idx[i]))
    # Turn any clashing hydroxyl before judging the result: the torsion is free,
    # so a collision here is a placement artefact, not a bad pose.
    positions = relax_hydroxyls(positions, bonded, report)

    keys = list(positions)
    clashes = 0
    closest, closest_pair = 9.9, None
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = keys[i], keys[j]
            if (a, b) in bonded:
                continue
            if not (a[1].upper().startswith('H') or b[1].upper().startswith('H')):
                continue
            d = dist(positions[a], positions[b])
            if d < closest:
                closest, closest_pair = d, (a, b)
            if d < 1.20:
                clashes += 1
    if clashes:
        raise Problem("%d placed hydrogens are within 1.20 A of a non-bonded "
                      "atom (closest %.2f A, %s%d-%s%d)."
                      % (clashes, closest, closest_pair[0][1], closest_pair[0][0],
                         closest_pair[1][1], closest_pair[1][0]))
    report.append("hydrogen clashes         none below 1.20 A (closest non-bonded "
                  "contact %.2f A)" % closest)

    resnames = {a['resnr']: a['resname'] for a in mol['atoms']}
    resids = {r: r - base for r in resnames}
    outdir = os.path.join(HERE, 'build', name)
    write_pdb(os.path.join(outdir, 'glycan.pdb'), mol, positions, resnames, resids)
    write_itp(os.path.join(outdir, 'glycan.itp'), mol['path'], mol['name'],
              len(mol['atoms']))
    nrestr = write_posre(os.path.join(outdir, 'posre_lig.itp'), mol)
    report.append("wrote                    build/%s/{glycan.pdb, glycan.itp, "
                  "posre_lig.itp}" % name)
    report.append("                         %d atoms, %d heavy atoms restrained "
                  "under -DPOSRES" % (len(mol['atoms']), nrestr))
    return mol['name']


def main():
    shared = read_conf(os.path.join(HERE, 'config.sh'))
    # config.sh SYSTEMS is empty by default: the systems/ directory is the list.
    allsys = shared.get('SYSTEMS', '').split()
    if not allsys:
        allsys = sorted(os.path.splitext(f)[0]
                        for f in os.listdir(os.path.join(HERE, 'systems'))
                        if f.endswith('.conf'))
    holo = []
    for s in allsys:
        c = read_conf(os.path.join(HERE, 'systems', s + '.conf'))
        if c.get('HAS_LIGAND') == '1':
            holo.append(s)
    wanted = sys.argv[1:] or holo
    bad = [s for s in wanted if s not in holo]
    if bad:
        print("Not holo systems: %s" % ', '.join(bad))
        print("Holo systems are: %s" % ', '.join(holo))
        return 1

    # EG1 keeps ligand topologies under ligands/<name>_charmm/ instead of a raw
    # charmm_gui/ dump: the same CHARMM-GUI .itp, but one directory per ligand so
    # more than one ligand can coexist. LIGAND comes from systems/<name>.conf.
    lig = os.environ.get('LIGAND', 'X6')
    root = os.path.join(HERE, 'ligands', '%s_charmm' % lig)
    if not os.path.isdir(root):
        sys.stderr.write(
            "FATAL: no topology directory ligands/%s_charmm/.\n"
            "       Every oligosaccharide needs a CHARMM-GUI topology once;\n"
            "       see ligands/README.md.\n" % lig)
        sys.exit(1)
    print("=" * 74)
    print("CHARMM-GUI glycan topology")
    print("=" * 74)
    try:
        mol = find_glycan_itp(root)
    except Problem as exc:
        print("  *** %s" % exc)
        return 1

    counts = {}
    for a in mol['atoms']:
        e = element_of(a['name'], a['type'])
        counts[e] = counts.get(e, 0) + 1
    if (counts.get('C'), counts.get('H'), counts.get('O')) != (EXPECT_C, EXPECT_H, EXPECT_O):
        print("  *** the 105-atom moleculetype is %s, not C%d H%d O%d. That is "
              "not xylohexaose." % (counts, EXPECT_C, EXPECT_H, EXPECT_O))
        return 1
    qtot = sum(a['charge'] for a in mol['atoms'])
    if abs(qtot) > 0.01:
        print("  *** net charge %+.3f. Xylohexaose is neutral; a non-zero charge "
              "means the topology is not what it should be." % qtot)
        return 1

    resnames = sorted({a['resname'] for a in mol['atoms']})
    print("  source                   %s" % os.path.relpath(mol['path'], HERE))
    print("  moleculetype             %s" % mol['name'])
    print("  atoms                    %d  (C%d H%d O%d), net charge %+.3f"
          % (len(mol['atoms']), counts['C'], counts['H'], counts['O'], qtot))
    print("  residues                 %d  (%s)"
          % (len({a['resnr'] for a in mol['atoms']}), ', '.join(resnames)))
    print("  bonds                    %d" % len(mol['bonds']))

    found = find_template_coords(root, mol)
    tcoords = None
    if found:
        raw, _n, path = found
        base = min(a['resnr'] for a in mol['atoms']) - 1
        tcoords = {(r - base, n): xyz for (r, n), xyz in raw.items()}
        print("  coordinates              %s" % os.path.relpath(path, HERE))
        print("                           (used only for hydroxyl H torsions)")
    else:
        print("  coordinates              none found -- hydroxyl hydrogens will be")
        print("                           built anti. Harmless: they are free to")
        print("                           relax during the restrained NVT stage.")
    print("")

    failures = []
    for name in wanted:
        print("=" * 74)
        print("%s" % name)
        print("=" * 74)
        report = []
        try:
            build_one(name, mol, tcoords, report)
            for line in report:
                print("  " + line)
        except Problem as exc:
            print("  *** REFUSED: %s" % exc)
            failures.append((name, str(exc)))
        print("")

    print("=" * 74)
    if failures:
        print("%d of %d holo system(s) FAILED:" % (len(failures), len(wanted)))
        for n, w in failures:
            print("  %-12s %s" % (n, w.splitlines()[0]))
        return 1
    print("All %d holo system(s) built. Next:  ./prep_system.sh <name>" % len(wanted))
    return 0


if __name__ == '__main__':
    sys.exit(main())
