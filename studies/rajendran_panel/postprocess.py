#!/usr/bin/env python3
"""LaTeX tables of the panel study (Rajendran 2010, Secs. 5.5-5.7).

Reads results/patch.csv, results/convergence.csv (matlab/rp_run.m) and the
values transcribed from the paper (results/paper_values.csv,
results/paper_rates.csv); writes report/generated/*.tex.
"""
import csv
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
GEN = os.path.join(HERE, "report", "generated")
THETAS = [0, 90, 135, 175, 225, 270]
MESH = dict(zip(THETAS, "abcdef"))


def f4(v):
    """Displacement for the tables; NaN = element could not be built
    (coincident nodes or singular metric interpolation, error of the code
    base)."""
    v = float(v)
    return "err." if math.isnan(v) else f"{v:.4f}"


def read(name):
    with open(os.path.join(RES, name)) as f:
        return list(csv.DictReader(l for l in f if not l.startswith("#")))


def main():
    os.makedirs(GEN, exist_ok=True)
    P = {(r["variant"], int(r["theta"]), r["load"]): r for r in read("patch.csv")}
    C = {(r["variant"], int(r["theta"]), int(r["N"])): r for r in read("convergence.csv")}
    paper = {}
    for r in read("paper_values.csv"):
        paper[(int(r["table"]), int(r["theta"]), int(r["N"]))] = r

    # Tables 14-16: ours vs paper, p = 2 (9-node)
    for load, comp, tab, exact in (("tension", "uxA", 14, "1.0"), ("moment", "uyA", 15, "1.5"),
                                   ("shear", "uyA", 16, "3.3125")):
        L = ["\\begin{tabular}{l rr rr rr rr rr}", "\\toprule",
             " & \\multicolumn{2}{c}{PP9} & \\multicolumn{2}{c}{MM9} & \\multicolumn{2}{c}{PM9 (PM-QUAD9$^\\ast$)}"
             " & \\multicolumn{2}{c}{\\PM, $p=2$} & \\multicolumn{2}{c}{\\PH, $p=2$}\\\\",
             "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}",
             "mesh & present & paper & present & paper & present & paper & skew & code & centr. & code\\\\", "\\midrule"]
        for t in THETAS:
            pv = paper[(tab, t, 1)]
            cells = []
            for v, key in (("PP9", "PP9"), ("MM9", "MM9"), ("PM9", "PM9")):
                ours = float(P[(v, t, load)][comp])
                ref = float(pv[key])
                mark = "" if abs(ours - ref) <= 1.5e-4 else "$^\\ast$"
                cells += [f"{ours:.4f}", f"{pv[key]}{mark}"]
            cells += [f4(P[('PM9skew', t, load)][comp]), f4(P[('PMusem', t, load)][comp]),
                      f4(P[('PH9', t, load)][comp]), f4(P[('PHusem', t, load)][comp])]
            L.append(f"{MESH[t]} ({t}$^\\circ$) & " + " & ".join(cells) + "\\\\")
        L += ["\\bottomrule", "\\end{tabular}"]
        with open(os.path.join(GEN, f"patch_{load}.tex"), "w") as f:
            f.write("\n".join(L) + "\n")

    # higher p, shear: u_A and nodal error
    L = ["\\begin{tabular}{l" + "r" * len(THETAS) + "}", "\\toprule",
         " & " + " & ".join(f"{MESH[t]} ({t}$^\\circ$)" for t in THETAS) + "\\\\", "\\midrule"]
    for p in range(3, 9):
        for v, name in ((f"PP{p}", f"\\PP, $p={p}$"), (f"PM{p}", f"\\PM, $p={p}$"),
                        (f"PM{p}usem", f"\\PM$_\\mathrm{{code}}$, $p={p}$"),
                        (f"PH{p}", f"\\PH, $p={p}$"), (f"PH{p}usem", f"\\PH$_\\mathrm{{code}}$, $p={p}$")):
            cells = []
            for t in THETAS:
                e = float(P[(v, t, "shear")]["nodal"])
                cells.append("err." if math.isnan(e) else f"{e:.0e}".replace("e-0", "e-").replace("e+0", "e"))
            L.append(name + " & " + " & ".join(cells) + "\\\\")
        if p < 8:
            L.append("\\addlinespace[2pt]")
    L += ["\\bottomrule", "\\end{tabular}"]
    with open(os.path.join(GEN, "patch_highp.tex"), "w") as f:
        f.write("\n".join(L) + "\n")

    # Tables 17-19 and rates (Fig. 15)
    rates = {r["sequence"]: r for r in read("paper_rates.csv")}
    L = ["\\begin{tabular}{ll rr rr rr rr rr}", "\\toprule",
         " & & \\multicolumn{2}{c}{PP9} & \\multicolumn{2}{c}{MM9} & \\multicolumn{2}{c}{PM9}"
         " & \\multicolumn{2}{c}{\\PM, $p=2$} & \\multicolumn{2}{c}{\\PH, $p=2$}\\\\",
         "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}\\cmidrule(lr){11-12}",
         "seq. & mesh & present & paper & present & paper & present & paper & skew & code & centr. & code\\\\", "\\midrule"]
    R = []
    for t, seq, tab in ((0, "A", 17), (90, "B", 18), (135, "C", 19)):
        for N in (1, 2, 4):
            pv = paper[(tab, t, N)]
            cells = []
            for v in ("PP9", "MM9", "PM9"):
                ours = float(C[(v, t, N)]["uyA"])
                ref = float(pv[v])
                mark = "" if abs(ours - ref) <= 1.5e-4 else "$^\\ast$"
                cells += [f"{ours:.4f}", f"{pv[v]}{mark}"]
            cells += [f4(C[('PM9skew', t, N)]['uyA']), f4(C[('PMusem', t, N)]['uyA']),
                      f4(C[('PH9', t, N)]['uyA']), f4(C[('PHusem', t, N)]['uyA'])]
            L.append((seq if N == 1 else "") + f" & ${N}\\times{N}$ & " + " & ".join(cells) + "\\\\")
        L.append("\\midrule" if seq != "C" else "\\bottomrule")
        # least-squares rates of the energy-norm error (Eq. (24))
        row = [seq]
        for v in ("PP9", "PM9", "PM9skew", "PMusem", "PH9", "PHusem", "PP3", "PM3", "PH3"):
            e = [math.sqrt(abs(float(C[(v, t, N)]["energy"]))) for N in (1, 2, 4)]
            if any(math.isnan(x) for x in e):
                row.append("err.")
            elif min(e) <= 0 or min(e) < 1e-12:
                row.append("--")
            else:
                row.append(f"{np.polyfit(np.log([1, .5, .25]), np.log(e), 1)[0]:.2f}")
        row.insert(2, rates[seq]["PP9"])
        row.insert(4, rates[seq]["PM9"])
        R.append(row)
    L.append("\\end{tabular}")
    with open(os.path.join(GEN, "convergence.tex"), "w") as f:
        f.write("\n".join(L) + "\n")
    L = ["\\begin{tabular}{l rr rr rr rr rrr}", "\\toprule",
         " & \\multicolumn{2}{c}{PP9} & \\multicolumn{2}{c}{PM9} & \\multicolumn{2}{c}{\\PM, $p=2$}"
         " & \\multicolumn{2}{c}{\\PH, $p=2$} & \\multicolumn{3}{c}{$p=3$}\\\\",
         "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-12}",
         "seq. & present & paper & present & paper & skew & code & centr. & code & \\PP & \\PM & \\PH\\\\", "\\midrule"]
    for row in R:
        L.append(" & ".join(row) + "\\\\")
    L += ["\\bottomrule", "\\end{tabular}"]
    with open(os.path.join(GEN, "rates.tex"), "w") as f:
        f.write("\n".join(L) + "\n")

    # mesh figure (gmsh meshes), full text width: Meshes a-f in 3 x 2 panels,
    # sequences B and C in a row of 4 panels
    def meshgroup(cols, rows, w, sep, panels):
        G = ["\\begin{tikzpicture}",
             f"\\begin{{groupplot}}[group style={{group size={cols} by {rows}, horizontal sep={sep}\\textwidth,"
             " vertical sep=7mm}, scale only axis,"
             f" width={w}\\textwidth, height={w}\\textwidth, hide axis, xmin=-0.03, xmax=2.03, ymin=-1.03, ymax=1.03,"
             " title style={font=\\small, yshift=-0.5ex}]"]
        for N, t, ttl in panels:
            G.append(f"\\nextgroupplot[title={{{ttl}}}]")
            G.append(f"\\addplot[black, line width=0.4pt] table {{../mesh/plotdata/mesh_N{N}_t{t:03d}.dat}};")
        return G + ["\\end{groupplot}", "\\end{tikzpicture}"]
    L = meshgroup(3, 2, 0.3, 0.05, [(1, t, f"Mesh {MESH[t]}: $\\theta={t}^\\circ$") for t in THETAS])
    L += ["", "\\medskip", ""]
    L += meshgroup(4, 1, 0.22, 0.04, [(2, 90, "B, $2\\times2$"), (4, 90, "B, $4\\times4$"),
                                      (2, 135, "C, $2\\times2$"), (4, 135, "C, $4\\times4$")])
    with open(os.path.join(GEN, "fig_meshes.tex"), "w") as f:
        f.write("\n".join(L) + "\n")
    print("tables written to", GEN)


if __name__ == "__main__":
    main()
