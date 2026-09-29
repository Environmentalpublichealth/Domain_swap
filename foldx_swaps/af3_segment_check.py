#!/usr/bin/env python3
"""Flag swap mutants whose swapped segment AF3 predicts to refold.

FoldX BuildModel keeps the host crystal backbone fixed, so its ddG is only
meaningful if the swapped segment keeps roughly the host conformation. AF3 is
used here only as a warning flag, never as a structure to score.

For every mutant in mutants.tsv, finds the AF3 job whose model_0 sequence
equals the FASTA sequence, then for each of its 5 models:
  moved   CA RMSD of the segment (+2 flanking residues) from the host crystal
          after superposing everything except segment +/-3 residues
  pLDDT   mean CA pLDDT over the segment
The AF3 wild-type model of the host gives the baseline for "moved".

Flag:  AF3_refold         every model moved > 2 A   -> do not trust ddG
       low_AF3_confidence segment pLDDT < 70 in any model -> backbone uncertain
       native_backbone    otherwise

Needs gemmi + numpy (vina env):
    /opt/anaconda3/envs/vina/bin/python af3_segment_check.py
Writes af3_segment_check.tsv (read by compare_mutants.py / compare_singles.py).
"""
import glob
import os
import sys
import zipfile

import gemmi
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EG1 = os.path.dirname(HERE)
FASTA = os.path.join(EG1, "1KS5_1XYN_region_swap_10_sequences.fasta")
AF3_WT = {"1XYN": os.path.join(EG1, "models", "1xyn_apo_model.pdb"),
          "1KS5": os.path.join(EG1, "models", "1ks5_WT_AF3_model_0.pdb")}
MOVED_MAX, PLDDT_MIN = 2.0, 70.0


def parse(st):
    ca, seq, pl = {}, {}, {}
    for r in st[0][0]:
        a = r.find_atom("CA", "*")
        if a:
            n = r.seqid.num
            ca[n] = np.array(a.pos.tolist())
            seq[n] = gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
            pl[n] = a.b_iso
    return ca, "".join(seq[k] for k in sorted(seq)), pl


def from_cif_text(text):
    st = gemmi.make_structure_from_block(gemmi.cif.read_string(text).sole_block())
    st.setup_entities()
    return parse(st)


# AF3 jobs: folders fold_*/ (EG1 and SwapMutants_AF3_JY) and fold_*.zip archives.
# Each job -> list of 5 loaders.
jobs = {}
for d in glob.glob(os.path.join(EG1, "fold_*")) + glob.glob(os.path.join(EG1, "SwapMutants_AF3_JY", "fold_*")):
    if os.path.isdir(d):
        files = sorted(glob.glob(os.path.join(d, "*_model_[0-4].cif")))
        if files:
            jobs[os.path.basename(d)] = [lambda f=f: parse(gemmi.read_structure(f)) for f in files]
    elif d.endswith(".zip"):
        z = zipfile.ZipFile(d)
        names = sorted(n for n in z.namelist() if n.split("/")[-1].split("_model_")[-1] in {f"{i}.cif" for i in range(5)})
        if names:
            jobs.setdefault(os.path.basename(d)[:-4],
                            [lambda z=z, n=n: from_cif_text(z.read(n).decode()) for n in names])
job_seq = {j: loaders[0]()[1] for j, loaders in jobs.items()}

fasta = {}
for blk in open(FASTA).read().split(">")[1:]:
    h, *s = blk.strip().split("\n")
    fasta[h.split("|")[0].strip()] = "".join(s)


def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - Pc).T @ (Q - Qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return lambda X: (R @ (X - Pc).T).T + Qc


def moved(model, xtal, s, e):
    ca = model[0]
    core = [k for k in xtal if (k < s - 3 or k > e + 3) and k in ca]
    fit = kabsch(np.array([ca[k] for k in core]), np.array([xtal[k] for k in core]))
    seg = [k for k in range(s - 2, e + 3) if k in ca and k in xtal]
    A, B = fit(np.array([ca[k] for k in seg])), np.array([xtal[k] for k in seg])
    return float(np.sqrt(((A - B) ** 2).sum(1).mean()))


xtal = {bb: parse(gemmi.read_structure(os.path.join(HERE, "repair", f"{bb}_Repair.pdb")))[0] for bb in AF3_WT}
wt = {bb: parse(gemmi.read_structure(p)) for bb, p in AF3_WT.items()}

out = ["mutant\taf3_job\tmoved_min\tmoved_max\taf3_wt_moved\tplddt_min\tplddt_max\tflag"]
for line in open(os.path.join(HERE, "mutants.tsv")).read().splitlines()[1:]:
    name, bb, _, span, _ = line.split("\t")
    s, e = map(int, span.split("-"))
    match = [j for j, sq in job_seq.items() if sq == fasta[name]]
    if not match:
        out.append(f"{name}\tNONE\tNA\tNA\tNA\tNA\tNA\tno_AF3_model")
        continue
    models = [ld() for ld in jobs[match[0]]]
    mv = [moved(m, xtal[bb], s, e) for m in models]
    pl = [float(np.mean([m[2][k] for k in range(s, e + 1)])) for m in models]
    flag = ("AF3_refold" if min(mv) > MOVED_MAX else
            "low_AF3_confidence" if min(pl) < PLDDT_MIN else "native_backbone")
    out.append(f"{name}\t{match[0]}\t{min(mv):.1f}\t{max(mv):.1f}\t{moved(wt[bb], xtal[bb], s, e):.1f}\t"
               f"{min(pl):.0f}\t{max(pl):.0f}\t{flag}")

open(os.path.join(HERE, "af3_segment_check.tsv"), "w").write("\n".join(out) + "\n")
for l in out:
    print(l.replace("\t", "  "))
