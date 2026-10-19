#!/usr/bin/env python3
"""Post-processing of the thick-cylinder study.

Reads results/raw/study_aXXX.csv (written by matlab/run_study.m) and writes
  results/tables/h_<problem>_<form>_a<alpha>.dat   h-convergence (rows: N)
  results/tables/p_<problem>_<form>_a<alpha>.dat   p-convergence (rows: p)
  results/tables/alpha_<problem>_<form>_N<N>.dat   error vs distortion angle
  results/tables/diag_N<N>.dat                     metric-function diagnostics
  results/tables/rates.csv                         observed convergence rates
  report/generated/*.tex                           LaTeX tables for the report

Rates are computed from consecutive meshes, r = log(e_N/e_2N)/log(2) (the
element size scales with 1/N). A rate is only reported if both errors are
above FLOOR and, for the PM element with p >= 5, above ROUNDOFF_FACTOR times
the error of the polynomial field on the same mesh (empirical round-off
level), otherwise '--'.
"""
import csv
import glob
import math
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "results", "raw")
TAB = os.path.join(HERE, "results", "tables")
GEN = os.path.join(HERE, "report", "generated")
FLOOR = {"energy": 1e-12, "L2": 1e-13, "nodal": 1e-12}
# For p >= 5 the polynomial field lies in the PM trial space: its computed
# error is the round-off level of that mesh and degree. Rates of the PM
# element are only reported if both errors exceed ROUNDOFF_FACTOR times it.
ROUNDOFF_FACTOR = 10.0
PROBLEMS = ["lame", "poly", "smooth"]
FORMS = ["PM", "PP"]


def load():
    rows = []
    for f in sorted(glob.glob(os.path.join(RAW, "study_a*.csv"))):
        with open(f) as fh:
            for r in csv.DictReader(fh):
                for k in r:
                    if k not in ("problem", "form"):
                        r[k] = float(r[k])
                for k in ("alpha", "N", "p", "ndof", "negDetJ"):
                    r[k] = int(r[k])
                rows.append(r)
    return rows


def index(rows):
    D = {}
    for r in rows:
        D[(r["problem"], r["form"], r["alpha"], r["N"], r["p"])] = r
    return D


def rate(e1, e2, key, f1=0.0, f2=0.0):
    """Observed rate; None if an error is below the round-off level.
    f1, f2: empirical round-off floors of the two meshes."""
    if e1 is None or e2 is None:
        return None
    if min(e1, e2) < FLOOR[key] or e1 < f1 or e2 < f2:
        return None
    return math.log(e1 / e2) / math.log(2.0)


def sci(v, small=1e3):
    """Compact number for tables: 3 significant digits below small, else
    mantissa with one digit and exponent (e.g. 1.4e4)."""
    if v < small and v >= 1e-2:
        return f"{v:.3g}"
    m, e = f"{v:.1e}".split("e")
    return f"{m}e{int(e)}"


def fmt(v, nd=2):
    return "nan" if v is None else f"{v:.{nd}f}"


def main():
    rows = load()
    D = index(rows)
    alphas = sorted({r["alpha"] for r in rows})
    Ns = sorted({r["N"] for r in rows})
    ps = sorted({r["p"] for r in rows})
    os.makedirs(TAB, exist_ok=True)
    os.makedirs(GEN, exist_ok=True)
    get = lambda pr, fo, a, N, p, k: D[(pr, fo, a, N, p)][k] if (pr, fo, a, N, p) in D else None

    # h- and p-convergence tables
    for pr in PROBLEMS:
        for fo in FORMS:
            for a in alphas:
                with open(os.path.join(TAB, f"h_{pr}_{fo}_a{a:03d}.dat"), "w") as f:
                    cols = ["N", "h"] + [f"{k}{p}" for k in ("eE", "eL", "eN") for p in ps]
                    f.write(" ".join(cols) + "\n")
                    for N in Ns:
                        vals = [N, 1.0 / N]
                        for k in ("energy", "L2", "nodal"):
                            vals += [get(pr, fo, a, N, p, k) for p in ps]
                        f.write(" ".join("nan" if v is None else f"{v:.6e}" for v in vals) + "\n")
                with open(os.path.join(TAB, f"p_{pr}_{fo}_a{a:03d}.dat"), "w") as f:
                    cols = ["p"] + [f"{k}N{N}" for k in ("eE", "eL", "eN", "dof") for N in Ns]
                    f.write(" ".join(cols) + "\n")
                    for p in ps:
                        vals = [p]
                        for k in ("energy", "L2", "nodal", "ndof"):
                            vals += [get(pr, fo, a, N, p, k) for N in Ns]
                        f.write(" ".join("nan" if v is None else f"{v:.6e}" for v in vals) + "\n")
            for N in Ns:
                with open(os.path.join(TAB, f"alpha_{pr}_{fo}_N{N}.dat"), "w") as f:
                    f.write("alpha " + " ".join(f"eE{p}" for p in ps) + " " + " ".join(f"eN{p}" for p in ps) + "\n")
                    for a in alphas:
                        f.write(f"{a} " + " ".join(
                            "nan" if get(pr, fo, a, N, p, k) is None
                            else f"{get(pr, fo, a, N, p, k):.6e}" for k in ("energy", "nodal") for p in ps) + "\n")

    # diagnostics (independent of the problem)
    for N in Ns:
        with open(os.path.join(TAB, f"diag_N{N}.dat"), "w") as f:
            cols = ["alpha"] + [f"{k}{p}" for k in ("leb", "lebGLL", "cKePM", "cKePP", "cKPM", "cKPP", "neg", "minJ") for p in ps]
            f.write(" ".join(cols) + "\n")
            for a in alphas:
                vals = [a]
                vals += [get("lame", "PM", a, N, p, "lebesgueM") for p in ps]
                vals += [get("lame", "PM", a, N, p, "lebesgueGLL") for p in ps]
                vals += [get("lame", "PM", a, N, p, "condKe") for p in ps]
                vals += [get("lame", "PP", a, N, p, "condKe") for p in ps]
                vals += [get("lame", "PM", a, N, p, "condK") for p in ps]
                vals += [get("lame", "PP", a, N, p, "condK") for p in ps]
                vals += [get("lame", "PM", a, N, p, "negDetJ") for p in ps]
                vals += [get("lame", "PM", a, N, p, "minDetJ") for p in ps]
                f.write(" ".join("nan" if (v is None or (isinstance(v, float) and math.isnan(v)))
                                 else f"{v:.6e}" for v in vals) + "\n")

    # rates
    R = {}
    with open(os.path.join(TAB, "rates.csv"), "w") as f:
        f.write("problem,form,alpha,p,key," + ",".join(f"r{Ns[i]}_{Ns[i+1]}" for i in range(len(Ns) - 1)) + "\n")
        for pr in PROBLEMS:
            for fo in FORMS:
                for a in alphas:
                    for p in ps:
                        for k in ("energy", "L2", "nodal"):
                            def floor(N):
                                if fo != "PM" or p < 5 or pr == "poly":
                                    return 0.0
                                e = get("poly", "PM", a, N, p, k)
                                return 0.0 if e is None else ROUNDOFF_FACTOR * e
                            rr = [rate(get(pr, fo, a, Ns[i], p, k), get(pr, fo, a, Ns[i + 1], p, k), k,
                                       floor(Ns[i]), floor(Ns[i + 1]))
                                  for i in range(len(Ns) - 1)]
                            R[(pr, fo, a, p, k)] = rr
                            f.write(f"{pr},{fo},{a},{p},{k}," + ",".join(fmt(x) for x in rr) + "\n")

    # LaTeX tables ---------------------------------------------------------
    def neg(a, N):
        return get("lame", "PM", a, N, ps[0], "negDetJ") or 0

    def head(ncol):
        return ["\\begin{tabular}{l" + "r" * ncol + "}", "\\toprule",
                "& \\multicolumn{%d}{c}{$p$}\\\\" % ncol, "\\cmidrule(lr){2-%d}" % (ncol + 1),
                " & " + " & ".join(str(p) for p in ps) + "\\\\", "\\midrule"]

    def alabel(a, Nlist):
        # dagger: folded elements (det J <= 0) on one of the meshes
        return f"$\\alpha={a}^\\circ$" + ("$^\\dagger$" if any(neg(a, N) for N in Nlist) else "")

    FORM_NAME = {"PM": "\\PM\\ (unsymmetric)", "PP": "\\PP\\ (conventional)"}
    for pr in PROBLEMS:
        for k, kname, alist in (("energy", "E", [a for a in alphas if a <= 30]),
                                ("nodal", "N", alphas)):
            lines = head(len(ps))
            for fo in FORMS:
                lines.append("\\multicolumn{%d}{l}{%s}\\\\" % (len(ps) + 1, FORM_NAME[fo]))
                for a in alist:
                    cells = []
                    for p in ps:
                        x = R[(pr, fo, a, p, k)][-1]       # finest pair N = 4 -> 8
                        if pr == "poly" and fo == "PM" and p >= 5:
                            cells.append("ex.")
                        elif x is None:
                            cells.append("--")
                        else:
                            cells.append(f"{x:.2f}")
                    lines.append(alabel(a, Ns[-2:]) + " & " + " & ".join(cells) + "\\\\")
                lines.append("\\midrule" if fo == "PM" else "\\bottomrule")
            lines.append("\\end{tabular}")
            with open(os.path.join(GEN, f"rates_{pr}_{kname}.tex"), "w") as f:
                f.write("\n".join(lines) + "\n")

    # distortion sensitivity e_n(alpha)/e_n(0) on the finest mesh
    Nf = Ns[-1]
    for pr in PROBLEMS:
        lines = head(len(ps))
        for fo in FORMS:
            lines.append("\\multicolumn{%d}{l}{%s}\\\\" % (len(ps) + 1, FORM_NAME[fo]))
            for a in alphas[1:]:
                cells = []
                for p in ps:
                    e0 = get(pr, fo, 0, Nf, p, "nodal")
                    e1 = get(pr, fo, a, Nf, p, "nodal")
                    r = e1 / e0
                    cells.append(sci(r))
                lines.append(alabel(a, [Nf]) + " & " + " & ".join(cells) + "\\\\")
            lines.append("\\midrule" if fo == "PM" else "\\bottomrule")
        lines.append("\\end{tabular}")
        with open(os.path.join(GEN, f"ratio_{pr}.tex"), "w") as f:
            f.write("\n".join(lines) + "\n")

    # polynomial field, PM, p >= 5: nodal and energy error (round-off level)
    pe = [p for p in ps if p >= 5]
    lines = ["\\begin{tabular}{l" + "r" * (2 * len(pe)) + "}", "\\toprule",
             " & \\multicolumn{%d}{c}{$N=1$} & \\multicolumn{%d}{c}{$N=8$}\\\\" % (len(pe), len(pe)),
             "\\cmidrule(lr){2-%d}\\cmidrule(lr){%d-%d}" % (len(pe) + 1, len(pe) + 2, 2 * len(pe) + 1),
             " & " + " & ".join(f"$p={p}$" for p in pe) + " & " + " & ".join(f"$p={p}$" for p in pe) + "\\\\",
             "\\midrule"]
    for a in alphas:
        cells = [sci(get('poly', 'PM', a, N, p, 'nodal'), small=0)
                 for N in (Ns[0], Ns[-1]) for p in pe]
        lines.append(f"$\\alpha={a}^\\circ$ & " + " & ".join(cells) + "\\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    with open(os.path.join(GEN, "poly_exact.tex"), "w") as f:
        f.write("\n".join(lines) + "\n")

    # diagnostics on the finest mesh: Lebesgue constant of the metric
    # interpolation and element condition number ratio PM/PP
    for key, fname, fmtf in (("lebesgueM", "diag_leb.tex", lambda v: sci(v, 1e3)),
                             ("condratio", "diag_cond.tex", lambda v: sci(v, 1e3))):
        lines = head(len(ps))
        for a in alphas:
            cells = []
            for p in ps:
                if key == "condratio":
                    v = get("lame", "PM", a, Nf, p, "condKe") / get("lame", "PP", a, Nf, p, "condKe")
                else:
                    v = get("lame", "PM", a, Nf, p, key)
                cells.append(fmtf(v))
            lines.append(alabel(a, [Nf]) + " & " + " & ".join(cells) + "\\\\")
        if key == "lebesgueM":
            lines.append("\\midrule")
            lines.append("GLL & " + " & ".join(f"{get('lame', 'PM', 0, Nf, p, 'lebesgueGLL'):.3g}" for p in ps) + "\\\\")
        lines += ["\\bottomrule", "\\end{tabular}"]
        with open(os.path.join(GEN, fname), "w") as f:
            f.write("\n".join(lines) + "\n")

    # folded elements per mesh
    lines = ["\\begin{tabular}{l" + "r" * len(Ns) + "}", "\\toprule",
             " & " + " & ".join(f"$N={N}$" for N in Ns) + "\\\\", "\\midrule"]
    for a in alphas:
        lines.append(f"$\\alpha={a}^\\circ$ & " + " & ".join(
            f"{neg(a, N)}/{4 * N * N}" for N in Ns) + "\\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    with open(os.path.join(GEN, "folded.tex"), "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"{len(rows)} rows, alphas {alphas}, N {Ns}, p {ps}")


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------
# TikZ/PGFPlots figure bodies (explicit code; \foreach does not work reliably
# inside groupplot environments)
ALPHAS_ALL = [0, 15, 30, 45, 60, 90]
PS_ALL = list(range(2, 9))


def legend_matrix(anchor_cmd, ps=PS_ALL, extra=""):
    cells = " ".join(f"\\draw[pline={p}] (0,0) -- (4mm,0) node[black, right] {{$p={p}$}}; \\pgfmatrixnextcell"
                     for p in ps)
    return (f"\\matrix[draw=gray!50, fill=white, font=\\scriptsize, column sep=1mm, inner sep=1.5pt, "
            f"anchor=south] at {anchor_cmd} {{ {cells} {extra} \\\\ }};\n")


def fig_meshes():
    L = ["\\begin{tikzpicture}",
         "\\begin{groupplot}[group style={group size=6 by 2, horizontal sep=1mm, vertical sep=3mm}, meshplot]"]
    for N in (1, 2):
        for a in ALPHAS_ALL:
            L.append(f"\\nextgroupplot[title={{$\\alpha={a}^\\circ$}}]" if N == 1 else "\\nextgroupplot")
            L.append(f"\\addplot[black, line width=0.35pt] table {{\\pdat/mesh_N{N}_a{a:03d}.dat}};")
    L += ["\\end{groupplot}",
          "\\node[font=\\scriptsize, rotate=90, anchor=south] at (group c1r1.west) {$N=1$};",
          "\\node[font=\\scriptsize, rotate=90, anchor=south] at (group c1r2.west) {$N=2$};",
          "\\end{tikzpicture}"]
    return "\n".join(L) + "\n"


def fig_hconv(problem, key, ylabel, alphas=ALPHAS_ALL, ymin="1e-15", ymax="1e2"):
    """Rows PM/PP; columns: alphas (one problem) or a list of
    (problem, alpha, title prefix) tuples."""
    cols = [c if isinstance(c, tuple) else (problem, c, "") for c in alphas]
    n = len(cols)
    L = ["\\begin{tikzpicture}",
         f"\\begin{{groupplot}}[group style={{group size={n} by 2, horizontal sep=1.5mm, vertical sep=8mm, "
         "yticklabels at=edge left, xticklabels at=edge bottom}, "
         f"conv, width={21.6 / n + 0.3:.2f}cm, height=4.6cm, ymin={ymin}, ymax={ymax}]"]
    for fo, fl in (("PM", "\\PM"), ("PP", "\\PP")):
        for i, (pr, a, pre) in enumerate(cols):
            if pre:      # several problems: form in the y-label, problem in the title
                opts = f"title={{{pre}$\\alpha={a}^\\circ$}}"
                if i == 0:
                    opts += f", ylabel={{{fl}: {ylabel}}}"
            else:
                opts = f"title={{{fl}, $\\alpha={a}^\\circ$}}"
                if i == 0:
                    opts += f", ylabel={{{ylabel}}}"
            L.append(f"\\nextgroupplot[{opts}]")
            for p in PS_ALL:
                L.append(f"\\addplot[pline={p}] table[x=h, y={key}{p}] {{\\tab/h_{pr}_{fo}_a{a:03d}.dat}};")
    L.append("\\end{groupplot}")
    L.append(f"\\node[font=\\small] at ($(group c1r2.south)!0.5!(group c{n}r2.south)+(0,-8mm)$) "
             "{element size $h\\propto 1/N$};")
    L.append(legend_matrix(f"($(group c1r1.north)!0.5!(group c{n}r1.north)+(0,7mm)$)"))
    L.append("\\end{tikzpicture}")
    return "\n".join(L) + "\n"


def fig_alpha(problems, N, titles):
    n = len(problems)
    L = ["\\begin{tikzpicture}",
         f"\\begin{{groupplot}}[group style={{group size={n} by 1, horizontal sep=14mm}}, "
         "width=6.2cm, height=5.2cm, ymode=log, xtick={0,15,30,45,60,90}, xlabel={$\\alpha$ [$^\\circ$]}, "
         "tick label style={font=\\scriptsize}, label style={font=\\small}, title style={font=\\small}, "
         "unbounded coords=jump, every axis plot/.append style={line width=0.7pt, mark size=1.3pt}, grid=major, grid style={gray!20}]"]
    for i, (pr, t) in enumerate(zip(problems, titles)):
        L.append(f"\\nextgroupplot[title={{{t}}}" + (", ylabel={$e_n$}" if i == 0 else "") + "]")
        for p in (2, 4, 6, 8):
            L.append(f"\\addplot[pline={p}] table[x=alpha, y=eN{p}] {{\\tab/alpha_{pr}_PM_N{N}.dat}};")
            L.append(f"\\addplot[pdash={p}] table[x=alpha, y=eN{p}] {{\\tab/alpha_{pr}_PP_N{N}.dat}};")
    L.append("\\end{groupplot}")
    cells = " ".join(f"\\draw[pline={p}] (0,0) -- (4mm,0) node[black, right] {{$p={p}$}}; \\pgfmatrixnextcell" for p in (2, 4, 6, 8))
    cells += " \\draw[black, solid] (0,0) -- (4mm,0) node[black, right] {\\PM}; \\pgfmatrixnextcell"
    cells += " \\draw[black, dashed] (0,0) -- (4mm,0) node[black, right] {\\PP}; \\pgfmatrixnextcell"
    L.append(f"\\matrix[draw=gray!50, fill=white, font=\\scriptsize, column sep=1mm, inner sep=1.5pt, anchor=south] "
             f"at ($(group c1r1.north)!0.5!(group c{n}r1.north)+(0,9mm)$) {{ {cells} \\\\ }};")
    L.append("\\end{tikzpicture}")
    return "\n".join(L) + "\n"


def fig_pconv(problems, titles, N):
    n = len(problems)
    sty = {0: "black", 15: "p3", 30: "p5", 45: "p7"}
    L = ["\\begin{tikzpicture}",
         f"\\begin{{groupplot}}[group style={{group size={n} by 1, horizontal sep=14mm}}, "
         "width=6.2cm, height=5.2cm, ymode=log, xtick={2,...,8}, xlabel={$p$}, "
         "tick label style={font=\\scriptsize}, label style={font=\\small}, title style={font=\\small}, "
         "unbounded coords=jump, every axis plot/.append style={line width=0.7pt, mark size=1.3pt}, grid=major, grid style={gray!20}]"]
    for i, (pr, t) in enumerate(zip(problems, titles)):
        L.append(f"\\nextgroupplot[title={{{t}}}" + (", ylabel={$e_n$}" if i == 0 else "") + "]")
        for a, c in sty.items():
            L.append(f"\\addplot[color={c}, mark=*, mark options={{solid, fill={c}}}] table[x=p, y=eNN{N}] {{\\tab/p_{pr}_PM_a{a:03d}.dat}};")
            L.append(f"\\addplot[color={c}, mark=o, dashed, mark options={{solid}}] table[x=p, y=eNN{N}] {{\\tab/p_{pr}_PP_a{a:03d}.dat}};")
    L.append("\\end{groupplot}")
    cells = " ".join(f"\\draw[color={c}, line width=0.7pt] (0,0) -- (4mm,0) node[black, right] {{$\\alpha={a}^\\circ$}}; \\pgfmatrixnextcell"
                     for a, c in sty.items())
    cells += " \\draw[black, solid] (0,0) -- (4mm,0) node[black, right] {\\PM}; \\pgfmatrixnextcell"
    cells += " \\draw[black, dashed] (0,0) -- (4mm,0) node[black, right] {\\PP}; \\pgfmatrixnextcell"
    L.append(f"\\matrix[draw=gray!50, fill=white, font=\\scriptsize, column sep=1mm, inner sep=1.5pt, anchor=south] "
             f"at ($(group c1r1.north)!0.5!(group c{n}r1.north)+(0,9mm)$) {{ {cells} \\\\ }};")
    L.append("\\end{tikzpicture}")
    return "\n".join(L) + "\n"


def write_figures():
    os.makedirs(GEN, exist_ok=True)
    out = {
        "fig_meshes.tex": fig_meshes(),
        "fig_h_lame_eN.tex": fig_hconv("lame", "eN", "$e_n$"),
        "fig_h_smooth_eN.tex": fig_hconv("smooth", "eN", "$e_n$"),
        "fig_h_poly_eN.tex": fig_hconv("poly", "eN", "$e_n$"),
        "fig_h_eE.tex": fig_hconv(None, "eE", "$\\eE$",
                                  alphas=[("smooth", a, "MMS ") for a in (0, 15, 30)] +
                                         [("lame", a, "Lam\\'e ") for a in (0, 15, 30)]),
        "fig_alpha_N4.tex": fig_alpha(["lame", "smooth"], 4, ["Lam\\'e, $N=4$", "smooth MMS, $N=4$"]),
        "fig_alpha_N8.tex": fig_alpha(["lame", "smooth"], 8, ["Lam\\'e, $N=8$", "smooth MMS, $N=8$"]),
        "fig_p_N2.tex": fig_pconv(["lame", "smooth"], ["Lam\\'e, $N=2$", "smooth MMS, $N=2$"], 2),
        "fig_p_N8.tex": fig_pconv(["lame", "smooth"], ["Lam\\'e, $N=8$", "smooth MMS, $N=8$"], 8),
    }
    for k, v in out.items():
        with open(os.path.join(GEN, k), "w") as f:
            f.write(v)


if __name__ == "__main__":
    write_figures()
