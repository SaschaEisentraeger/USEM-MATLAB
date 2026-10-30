#!/usr/bin/env python3
"""Thick-cylinder tables of the manuscript.

Collects the results of
  studies/thick_cylinder/results/            distortion and convergence study
  data/unisolvence.txt                       studies/unisolvence_check.m
and writes data/*.dat (plot data) and generated/*.tex (tables).
"""
import csv
import math
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))  # repository root
CYL = os.path.join(ROOT, "studies", "thick_cylinder", "results")
# rerun of the thick-cylinder study with the calibrated Gauss rules (ms_cylinder.m)
CYL_RAW = os.path.join(CYL, "raw_calibrated")
DATA = os.path.join(HERE, "..", "data")
GEN = os.path.join(HERE, "..", "generated")
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


def num(v, nd=2):
    """Fixed-point number; negative values in math mode (minus sign, not hyphen)."""
    s = f"{v:.{nd}f}"
    return f"${s}$" if s.startswith("-") else s


def write(name, lines):
    with open(os.path.join(GEN, name), "w") as f:
        f.write("\n".join(lines) + "\n")


def read_cyl(alpha):
    with open(os.path.join(CYL_RAW, f"study_a{alpha:03d}.csv")) as f:
        return {(r["problem"], r["form"], int(r["N"]), int(r["p"])): r for r in csv.DictReader(f)}


def read_interp(alpha):
    with open(os.path.join(CYL_RAW, f"interp_a{alpha:03d}.csv")) as f:
        return {(r["problem"], r["form"], int(r["N"]), int(r["p"])): float(r["energy"])
                for r in csv.DictReader(f)}


def write_convergence_tables(pr, fo, a):
    """data/cyl_h_* (rows: N) and data/cyl_p_* (rows: p), NaN: case rejected."""
    Ns, ps = (1, 2, 4, 8), range(2, 9)
    D = read_cyl(a)
    get = lambda N, p, k: float(D[(pr, fo, N, p)][k])
    line = lambda vals: " ".join("nan" if math.isnan(v) else f"{v:.6e}" for v in vals) + "\n"
    with open(os.path.join(DATA, f"cyl_h_{pr}_{fo}_a{a:03d}.dat"), "w") as f:
        f.write(" ".join(["N", "h"] + [f"{k}{p}" for k in ("eE", "eL", "eN") for p in ps]) + "\n")
        for N in Ns:
            f.write(line([N, 1.0 / N] + [get(N, p, k) for k in ("energy", "L2", "nodal") for p in ps]))
    with open(os.path.join(DATA, f"cyl_p_{pr}_{fo}_a{a:03d}.dat"), "w") as f:
        f.write(" ".join(["p"] + [f"{k}N{N}" for k in ("eE", "eL", "eN", "dof") for N in Ns]) + "\n")
        for p in ps:
            f.write(line([p] + [get(N, p, k) for k in ("energy", "L2", "nodal", "ndof") for N in Ns]))


def main():
    os.makedirs(DATA, exist_ok=True)
    os.makedirs(GEN, exist_ok=True)
    # thick cylinder: h- and p-convergence tables
    for pr in ("smooth", "lame"):
        for fo in FORMS:
            for a in (0, 15, 30, 45):
                write_convergence_tables(pr, fo, a)

    # mesh plots (gmsh meshes of the two studies)
    for a in (0, 15, 30, 45, 60, 90):
        shutil.copy(os.path.join(ROOT, "studies", "thick_cylinder", "mesh", "plotdata", f"mesh_N2_a{a:03d}.dat"),
                    os.path.join(DATA, f"cyl_mesh_N2_a{a:03d}.dat"))
    for t in (0, 90, 135, 175, 225, 270):
        shutil.copy(os.path.join(ROOT, "studies", "rajendran_panel", "mesh", "plotdata", f"mesh_N1_t{t:03d}.dat"),
                    os.path.join(DATA, f"panel_mesh_N1_t{t:03d}.dat"))

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
                        cells.append(num(math.log2(e1 / e2)))
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
                        cells.append(num(r) if r < 100 else sci(r, 0))
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

    # unisolvence check
    with open(os.path.join(DATA, "unisolvence.txt")) as f:
        U = list(csv.DictReader(f, delimiter=" "))
    L = ["\\begin{tabular}{@{}rrrcccc@{}}", "\\toprule",
         "$\\alpha$ & $N$ & $p$ & $\\min\\mathrm{rcond}(\\bm V_{\\mathrm{T}})$ & $\\min\\mathrm{rcond}(\\bm S)$ & "
         "$\\min\\mathrm{rcond}(\\bm V_{\\mathrm{Q}})$ & singular\\\\", "\\midrule"]
    for r in U:
        L.append(f"${r['alpha']}^\\circ$ & {r['N']} & {r['p']} & {sci(float(r['rcondVT']))} & "
                 f"{sci(float(r['rcondA']))} & {sci(float(r['rcondVQ']))} & {r['nsingular']}\\\\")
    L += ["\\bottomrule", "\\end{tabular}"]
    write("tab_unisolvence.tex", L)
    print("data written to", DATA, "and", GEN)


if __name__ == "__main__":
    main()
