#!/usr/bin/env python3
"""Data and tables of the manuscript on metric trial functions.

Collects the results of
  docs/figures/data/metric-*.dat            element- and mesh-level studies
                                            (studies/metric_basis/*.m)
  studies/rajendran_panel/results/*.csv     Jacobian-unlimited patch tests
  studies/thick_cylinder/results/            distortion and convergence study
  results/unisolvence.txt                    matlab/unisolvence_check.m
and writes data/*.dat (plot data) and generated/*.tex (tables).
"""
import csv
import math
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
DOCS = os.path.join(ROOT, "docs", "figures", "data")
CYL = os.path.join(ROOT, "studies", "thick_cylinder", "results")
PANEL = os.path.join(ROOT, "studies", "rajendran_panel", "results")
DATA = os.path.join(HERE, "data")
GEN = os.path.join(HERE, "generated")
FORMS = ("PM", "PH", "PP")
NAME = {"PM": "$Q_p(\\hat x,\\hat y)$", "PH": "$T_p$ (hybrid)", "PP": "parametric"}


def sci(v, nd=1):
    """1.2e-3 -> 1.2\\cdot10^{-3} (math mode)."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "rej."
    if v == 0:
        return "$0$"
    e = math.floor(math.log10(abs(v)))
    m = v / 10 ** e
    if round(m, nd) >= 10:
        m, e = m / 10, e + 1
    return f"${m:.{nd}f}\\cdot10^{{{e}}}$" if e != 0 else f"${m:.{nd}f}$"


def write(name, lines):
    with open(os.path.join(GEN, name), "w") as f:
        f.write("\n".join(lines) + "\n")


def read_cyl(alpha):
    with open(os.path.join(CYL, "raw", f"study_a{alpha:03d}.csv")) as f:
        return {(r["problem"], r["form"], int(r["N"]), int(r["p"])): r for r in csv.DictReader(f)}


def read_interp(alpha):
    with open(os.path.join(CYL, "raw", f"interp_a{alpha:03d}.csv")) as f:
        return {(r["problem"], r["form"], int(r["N"]), int(r["p"])): float(r["energy"])
                for r in csv.DictReader(f)}


def main():
    os.makedirs(DATA, exist_ok=True)
    os.makedirs(GEN, exist_ok=True)
    # element- and mesh-level data (copied unchanged)
    for f in sorted(os.listdir(DOCS)):
        if f.startswith("metric-") and f.endswith(".dat"):
            shutil.copy(os.path.join(DOCS, f), os.path.join(DATA, f))
    # thick cylinder: h- and p-convergence tables (copied unchanged)
    for pr in ("smooth", "lame"):
        for fo in FORMS:
            for a in (0, 15, 30, 45):
                for kind in ("h", "p"):
                    src = os.path.join(CYL, "tables", f"{kind}_{pr}_{fo}_a{a:03d}.dat")
                    shutil.copy(src, os.path.join(DATA, f"cyl_{kind}_{pr}_{fo}_a{a:03d}.dat"))

    # mesh plots (gmsh meshes of the two studies)
    for a in (0, 15, 30, 45, 60, 90):
        shutil.copy(os.path.join(ROOT, "studies", "thick_cylinder", "mesh", "plotdata", f"mesh_N2_a{a:03d}.dat"),
                    os.path.join(DATA, f"cyl_mesh_N2_a{a:03d}.dat"))
    for t in (0, 90, 135, 175, 225, 270):
        shutil.copy(os.path.join(ROOT, "studies", "rajendran_panel", "mesh", "plotdata", f"mesh_N1_t{t:03d}.dat"),
                    os.path.join(DATA, f"panel_mesh_N1_t{t:03d}.dat"))

    # Table: asymmetry of K_e on the parallelogram (affine element)
    rows = {}
    with open(os.path.join(DOCS, "metric-basis-parallelogram.dat")) as f:
        for r in csv.DictReader(f, delimiter=" ", skipinitialspace=True):
            r = {k: v for k, v in r.items() if k}
            rows[int(float(r["p"]))] = r
    ps = (4, 6, 8, 10, 12, 16)
    L = ["\\begin{tabular}{@{}l" + "c" * len(ps) + "@{}}", "\\toprule",
         "$p$ & " + " & ".join(str(p) for p in ps) + "\\\\", "\\midrule"]
    for key, name in (("mon", "monomial"), ("lag", "Lagrange (scaled)")):
        L.append(name + " & " + " & ".join(sci(float(rows[p][f"asym_{key}"])) for p in ps) + "\\\\")
    L.append("\\midrule")
    for key, name in (("mon", "monomial"), ("lag", "Lagrange (scaled)")):
        L.append("$\\kappa_2(\\bm V)$, " + name + " & " + " & ".join(sci(float(rows[p][f"kappa_{key}"])) for p in ps) + "\\\\")
    L += ["\\bottomrule", "\\end{tabular}"]
    write("tab_affine.tex", L)

    # Table: Jacobian-unlimited patch test (shear field, cubic), max over Meshes b-f
    with open(os.path.join(PANEL, "patch.csv")) as f:
        P = [r for r in csv.DictReader(l for l in f if not l.startswith("#"))]
    L = ["\\begin{tabular}{@{}l" + "c" * 6 + "@{}}", "\\toprule",
         " & " + " & ".join(f"$p={p}$" for p in range(3, 9)) + "\\\\", "\\midrule"]
    for fo, name in (("PP", "parametric (\\PP)"), ("PM", "$Q_p(\\hat x,\\hat y)$ (\\PM)"), ("PH", "$T_p$ (\\PH)")):
        cells = []
        for p in range(3, 9):
            v = [float(r["nodal"]) for r in P if r["variant"] == f"{fo}{p}" and r["load"] == "shear"
                 and int(r["theta"]) != 0]
            cells.append(sci(max(v)) if not any(math.isnan(x) for x in v) else "rej.")
        L.append(name + " & " + " & ".join(cells) + "\\\\")
    L += ["\\bottomrule", "\\end{tabular}"]
    write("tab_panel.tex", L)

    # Thick cylinder tables
    D = {a: read_cyl(a) for a in (0, 15, 30, 45, 60, 90)}
    I = {a: read_interp(a) for a in (0, 15, 30, 45, 60, 90)}

    def g(a, pr, fo, N, p, k="energy"):
        return float(D[a][(pr, fo, N, p)][k])

    def floor(a, fo, N, p, k="energy"):
        # empirical round-off level: error of the polynomial field (p >= 5)
        f = 1e-12 if k == "energy" else 1e-12
        if fo != "PP" and p >= 5:
            f = max(f, 10 * g(a, "poly", fo, N, p, k))
        return f

    # energy rates N = 4 -> 8
    L = ["\\begin{tabular}{@{}ll" + "r" * 7 + "@{}}", "\\toprule",
         " & & " + " & ".join(f"$p={p}$" for p in range(2, 9)) + "\\\\", "\\midrule"]
    for pr, prn in (("smooth", "analytic"), ("lame", "Lam\\'e")):
        for a in (15, 30):
            for fo in FORMS:
                cells = []
                for p in range(2, 9):
                    e1, e2 = g(a, pr, fo, 4, p), g(a, pr, fo, 8, p)
                    if min(e1, e2) < floor(a, fo, 8, p) or math.isnan(e1) or math.isnan(e2):
                        cells.append("--")
                    else:
                        cells.append(f"{math.log2(e1 / e2):.2f}")
                lab = f"{prn}, $\\alpha={a}^\\circ$" if fo == "PM" else ""
                L.append(lab + " & " + "\\" + fo + " & " + " & ".join(cells) + "\\\\")
            L.append("\\addlinespace[2pt]")
    L.pop()
    L += ["\\bottomrule", "\\end{tabular}"]
    write("tab_rates.tex", L)

    # quasi-optimality: e_E(discrete)/e_E(interpolant), N = 8; min/max over all
    # non-folded cases above round-off
    L = ["\\begin{tabular}{@{}ll" + "r" * 7 + "@{}}", "\\toprule",
         " & & " + " & ".join(f"$p={p}$" for p in range(2, 9)) + "\\\\", "\\midrule"]
    for pr, prn in (("smooth", "analytic"), ("lame", "Lam\\'e")):
        for a in (15, 30):
            for fo in FORMS:
                cells = []
                for p in range(2, 9):
                    ei, ed = I[a][(pr, fo, 8, p)], g(a, pr, fo, 8, p)
                    if ei < floor(a, fo, 8, p):
                        cells.append("--")
                    else:
                        r = ed / ei
                        cells.append(f"{r:.2f}" if r < 100 else sci(r, 0))
                lab = f"{prn}, $\\alpha={a}^\\circ$" if fo == "PM" else ""
                L.append(lab + " & \\" + fo + " & " + " & ".join(cells) + "\\\\")
            L.append("\\addlinespace[2pt]")
    L.pop()
    L += ["\\bottomrule", "\\end{tabular}"]
    write("tab_quasiopt.tex", L)
    stats = {}
    for fo in FORMS:
        rs = []
        for a in (0, 15, 30):
            for N in (1, 2, 4, 8):
                for p in range(2, 9):
                    for pr in ("smooth", "lame"):
                        ei, ed = I[a][(pr, fo, N, p)], g(a, pr, fo, N, p)
                        if ei > floor(a, fo, N, p) and ei > 1e-11:
                            rs.append(ed / ei)
        stats[fo] = (len(rs), min(rs), max(rs))
    with open(os.path.join(GEN, "quasiopt_stats.txt"), "w") as f:
        for fo, (n, lo, hi) in stats.items():
            f.write(f"{fo}: {n} cases, ratio {lo:.3f} .. {hi:.3g}\n")

    # Lebesgue constants and kappa(V), N = 8 (max over the elements)
    L = ["\\begin{tabular}{@{}ll" + "r" * 7 + "@{}}", "\\toprule",
         " & & " + " & ".join(f"$p={p}$" for p in range(2, 9)) + "\\\\", "\\midrule"]
    for key, kname in (("lebesgue", "$\\Lambda$"), ("condV", "$\\kappa_2(\\bm V)$")):
        for a in (15, 30, 45):
            for fo in ("PM", "PH"):
                cells = [sci(g(a, "lame", fo, 8, p, key)) for p in range(2, 9)]
                lab = f"{kname}, $\\alpha={a}^\\circ$" if fo == "PM" else ""
                L.append(lab + " & \\" + fo + " & " + " & ".join(cells) + "\\\\")
        if key == "lebesgue":
            L.append(" & GLL & " + " & ".join(sci(g(0, "lame", "PP", 8, p, "lebesgueGLL"), 2) for p in range(2, 9)) + "\\\\")
            L.append("\\midrule")
    L += ["\\bottomrule", "\\end{tabular}"]
    write("tab_lebesgue.tex", L)

    # distortion sensitivity and accuracy at N = 8
    L = ["\\begin{tabular}{@{}ll" + "r" * 7 + "@{}}", "\\toprule",
         " & & " + " & ".join(f"$p={p}$" for p in range(2, 9)) + "\\\\", "\\midrule"]
    for pr, prn in (("smooth", "analytic"), ("lame", "Lam\\'e")):
        for a in (30, 45):
            for fo in FORMS:
                cells = [sci(g(a, pr, fo, 8, p, "nodal")) for p in range(2, 9)]
                lab = f"{prn}, $\\alpha={a}^\\circ$" if fo == "PM" else ""
                L.append(lab + " & \\" + fo + " & " + " & ".join(cells) + "\\\\")
            L.append("\\addlinespace[2pt]")
    L.pop()
    L += ["\\bottomrule", "\\end{tabular}"]
    write("tab_nodal.tex", L)

    # polynomial patch test on the cylinder meshes: max nodal error over N for p >= 5
    L = ["\\begin{tabular}{@{}l" + "r" * 6 + "@{}}", "\\toprule",
         " & " + " & ".join(f"$\\alpha={a}^\\circ$" for a in (0, 15, 30, 45, 60, 90)) + "\\\\", "\\midrule"]
    for fo in ("PM", "PH"):
        cells = []
        for a in (0, 15, 30, 45, 60, 90):
            v = [g(a, "poly", fo, N, p, "nodal") for N in (1, 2, 4, 8) for p in range(5, 9)]
            vv = [x for x in v if not math.isnan(x)]
            nrej = len(v) - len(vv)
            cells.append(sci(max(vv)) + (f" ({nrej})" if nrej else ""))
        L.append("\\" + fo + " & " + " & ".join(cells) + "\\\\")
    L += ["\\bottomrule", "\\end{tabular}"]
    write("tab_cylpoly.tex", L)

    # unisolvence check
    with open(os.path.join(HERE, "results", "unisolvence.txt")) as f:
        U = list(csv.DictReader(f, delimiter=" "))
    L = ["\\begin{tabular}{@{}rrrcccc@{}}", "\\toprule",
         "$\\alpha$ & $N$ & $p$ & $\\min\\mathrm{rcond}(\\bm V_T)$ & $\\min\\mathrm{rcond}(\\bm A)$ & "
         "$\\min\\mathrm{rcond}(\\bm V_Q)$ & singular\\\\", "\\midrule"]
    for r in U:
        L.append(f"${r['alpha']}^\\circ$ & {r['N']} & {r['p']} & {sci(float(r['rcondVT']))} & "
                 f"{sci(float(r['rcondA']))} & {sci(float(r['rcondVQ']))} & {r['nsingular']}\\\\")
    L += ["\\bottomrule", "\\end{tabular}"]
    write("tab_unisolvence.tex", L)
    print("data written to", DATA, "and", GEN)


if __name__ == "__main__":
    main()
