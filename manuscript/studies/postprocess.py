#!/usr/bin/env python3
"""Post-processing of the manuscript studies (no third-party packages).

- splits manuscript/data/conv-<tag>.dat into conv-<tag>-p<p>.dat (one file
  per degree, column h = 1/N added) for the convergence plots;
- writes LaTeX tables to manuscript/generated/:
    rates-<tag>.tex      observed energy-norm rates, finest pair of meshes
    patch.tex            patch-test errors on the distorted 3 x 2 meshes
    panel.tex            patch tests on the folded panel meshes (PM vs PH)
    spectrum.tex         spurious eigenvalues (count, max |Im|/|lambda|)
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
GEN = os.path.join(HERE, '..', 'generated')
os.makedirs(GEN, exist_ok=True)


def read_table(path):
    with open(path) as f:
        head = f.readline().split()
        rows = [[float(v) for v in line.split()] for line in f if line.strip()]
    return {k: [r[i] for r in rows] for i, k in enumerate(head)}


def write_table(path, cols, names):
    with open(path, 'w') as f:
        f.write(' '.join(names) + '\n')
        for i in range(len(cols[names[0]])):
            f.write(' '.join('%.6e' % cols[k][i] for k in names) + '\n')


def sci(v, digits=1):
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return '--'
    if v == 0:
        return '$0$'
    e = int(math.floor(math.log10(abs(v))))
    m = v / 10 ** e
    if round(m, digits) >= 10:
        m /= 10
        e += 1
    return '$%.*f{\\cdot}10^{%d}$' % (digits, m, e)


def num(v, fmt):
    """Fixed-point number; negative values in math mode (minus sign, not hyphen)."""
    s = fmt % v
    return '$%s$' % s if s.startswith('-') else s


def convergence(tag):
    path = os.path.join(DATA, 'conv-%s.dat' % tag)
    if not os.path.exists(path):
        return
    T = read_table(path)
    ps = sorted(set(int(p) for p in T['p']))
    names = list(T.keys())
    rates = {}
    for p in ps:
        idx = [i for i, q in enumerate(T['p']) if int(q) == p]
        cols = {k: [T[k][i] for i in idx] for k in names}
        cols['h'] = [1.0 / n for n in cols['N']]
        write_table(os.path.join(DATA, 'conv-%s-p%d.dat' % (tag, p)), cols, ['h'] + names)
        for f in ('PP', 'PM', 'PH'):
            e = cols['energy_' + f]
            r = [math.log2(e[i] / e[i + 1]) if e[i] > 0 and e[i + 1] > 0 else float('nan')
                 for i in range(len(e) - 1)]
            rates[(p, f)] = (r, e[-1])
    if tag == 'd030':  # folded meshes: convergence plot only
        return
    with open(os.path.join(GEN, 'rates-%s.tex' % tag), 'w') as f:
        f.write('\\begin{tabular}{@{}c' + 'rrr' * 3 + '@{}}\n\\toprule\n')
        f.write('& \\multicolumn{3}{c}{\\PP} & \\multicolumn{3}{c}{\\PM} '
                '& \\multicolumn{3}{c}{\\PH}\\\\\n')
        f.write('\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\n')
        f.write('$p$' + ' & $r_{4,8}$ & $r_{8,16}$ & $e_{\\mathrm{E}}^{16}$' * 3 + '\\\\\n\\midrule\n')
        for p in ps:
            line = '%d' % p
            for fm in ('PP', 'PM', 'PH'):
                r, e = rates[(p, fm)]
                line += ' & %s & %s & %s' % (num(r[-2], '%.2f'), num(r[-1], '%.2f'), sci(e))
            f.write(line + '\\\\\n')
        f.write('\\bottomrule\n\\end{tabular}\n')


def patch():
    meshes = [('curved15', 'curved $0.15$'), ('curved20', 'curved $0.20$'),
              ('combined30', 'angular + curved $0.30$')]
    have = [m for m in meshes if os.path.exists(os.path.join(DATA, 'spectrum-%s.dat' % m[0]))]
    if not have:
        return
    tabs = {m: read_table(os.path.join(DATA, 'spectrum-%s.dat' % m)) for m, _ in have}
    ps = [int(p) for p in tabs[have[0][0]]['p']]
    show = [p for p in ps if p in (2, 4, 6, 8, 10)]
    with open(os.path.join(GEN, 'patch.tex'), 'w') as f:
        f.write('\\begin{tabular}{@{}ll' + 'c' * len(show) + '@{}}\n\\toprule\n')
        f.write('mesh & & ' + ' & '.join('$p=%d$' % p for p in show) + '\\\\\n\\midrule\n')
        for m, label in have:
            T = tabs[m]
            for k, (col, name) in enumerate((('patch_cla', '\\PP'), ('patch_skew', '\\PM'),
                                             ('patch_hyb', '\\PH'))):
                first = ('\\multirow{3}{*}{%s}' % label) if k == 0 else ''
                vals = [sci(T[col][ps.index(p)]) for p in show]
                f.write('%s & %s & %s\\\\\n' % (first, name, ' & '.join(vals)))
            f.write('\\midrule\n' if m != have[-1][0] else '')
        f.write('\\bottomrule\n\\end{tabular}\n')
    with open(os.path.join(GEN, 'spectrum.tex'), 'w') as f:
        f.write('\\begin{tabular}{@{}ll' + 'c' * len(show) + '@{}}\n\\toprule\n')
        f.write('mesh & & ' + ' & '.join('$p=%d$' % p for p in show) + '\\\\\n\\midrule\n')
        for m, label in have:
            T = tabs[m]
            for k, (fm, name) in enumerate((('skew', '\\PM'), ('hyb', '\\PH'))):
                first = ('\\multirow{2}{*}{%s}' % label) if k == 0 else ''
                vals = ['%d / %s' % (int(T['neg_' + fm][ps.index(p)]),
                                     sci(T['imag_' + fm][ps.index(p)])) for p in show]
                f.write('%s & %s & %s\\\\\n' % (first, name, ' & '.join(vals)))
            f.write('\\midrule\n' if m != have[-1][0] else '')
        f.write('\\bottomrule\n\\end{tabular}\n')


def panel():
    path = os.path.join(DATA, 'panel', 'patch.csv')
    if not os.path.exists(path):
        return
    rows = list(csv.DictReader(open(path)))
    thetas = [0, 90, 135, 175, 225, 270]
    ps = [3, 5, 8]
    with open(os.path.join(GEN, 'panel.tex'), 'w') as f:
        f.write('\\begin{tabular}{@{}cl' + 'c' * len(thetas) + '@{}}\n\\toprule\n')
        f.write('$p$ & & ' + ' & '.join('$%d^\\circ$' % t for t in thetas) + '\\\\\n\\midrule\n')
        for p in ps:
            for k, (v, name) in enumerate((('PM%d' % p, '\\PM'), ('PH%d' % p, '\\PH'))):
                vals = []
                for t in thetas:
                    e = [float(r['nodal']) for r in rows
                         if r['variant'] == v and int(r['theta']) == t]
                    e = [x for x in e if not math.isnan(x)]
                    vals.append(sci(max(e)) if e else '--')
                first = ('\\multirow{2}{*}{%d}' % p) if k == 0 else ''
                f.write('%s & %s & %s\\\\\n' % (first, name, ' & '.join(vals)))
            if p != ps[-1]:
                f.write('\\midrule\n')
        f.write('\\bottomrule\n\\end{tabular}\n')


def growth():
    meshes = [('curved15', 'curved $0.15$'), ('curved20', 'curved $0.20$'),
              ('combined30', 'angular + curved $0.30$')]
    have = [m for m in meshes if os.path.exists(os.path.join(DATA, 'growth-%s.dat' % m[0]))]
    if not have:
        return
    tabs = {m: read_table(os.path.join(DATA, 'growth-%s.dat' % m)) for m, _ in have}
    ps = [int(p) for p in tabs[have[0][0]]['p']]
    show = [p for p in ps if p in (2, 4, 6, 8)]
    with open(os.path.join(GEN, 'growth.tex'), 'w') as f:
        f.write('\\begin{tabular}{@{}ll' + 'cc' * len(show) + '@{}}\n\\toprule\n')
        f.write('& & ' + ' & '.join('\\multicolumn{2}{c}{$p=%d$}' % p for p in show) + '\\\\\n')
        f.write(''.join('\\cmidrule(lr){%d-%d}' % (3 + 2 * i, 4 + 2 * i) for i in range(len(show))) + '\n')
        f.write('mesh & & ' + ' & '.join(['$\\mu$ & $\\mu_\\Delta$'] * len(show)) + '\\\\\n\\midrule\n')
        for m, label in have:
            T = tabs[m]
            for k, (fm, name) in enumerate((('PM', '\\PM'), ('PH', '\\PH'))):
                first = ('\\multirow{2}{*}{%s}' % label) if k == 0 else ''
                vals = []
                for p in show:
                    i = ps.index(p) if p in ps else None
                    vals += [sci(T['gSemi_' + fm][i]) if i is not None else '--',
                             sci(T['gTrap_' + fm][i]) if i is not None else '--']
                f.write('%s & %s & %s\\\\\n' % (first, name, ' & '.join(vals)))
            f.write('\\midrule\n' if m != have[-1][0] else '')
        f.write('\\bottomrule\n\\end{tabular}\n')


CYL = os.path.join(HERE, '..', '..', 'studies', 'thick_cylinder', 'results', 'raw_calibrated')


def cylinder(raw=CYL, problem='lame'):
    """Thick cylinder (study of the code base): errors on the finest mesh versus the
    distortion angle, observed rates, Lebesgue constants."""
    alphas = [0, 15, 30, 45, 60, 90]
    rows = []
    for a in alphas:
        f = os.path.join(raw, 'study_a%03d.csv' % a)
        if os.path.exists(f):
            rows += list(csv.DictReader(open(f)))
    if not rows:
        return
    val = lambda r, k: float(r[k]) if r[k] not in ('', 'NaN', 'nan') else float('nan')
    def get(a, N, p, form, key, prob=problem):
        for r in rows:
            if (int(r['alpha']) == a and int(r['N']) == N and int(r['p']) == p
                    and r['form'] == form and r['problem'] == prob):
                return val(r, key)
        return float('nan')
    for p in (4, 6, 8):
        cols = {'alpha': [], 'PP': [], 'PM': [], 'PH': [], 'nPP': [], 'nPM': [], 'nPH': [],
                'lebPM': [], 'lebPH': []}
        for a in alphas:
            cols['alpha'].append(a)
            for fm in ('PP', 'PM', 'PH'):
                cols[fm].append(get(a, 8, p, fm, 'energy'))
                cols['n' + fm].append(get(a, 8, p, fm, 'nodal'))
            cols['lebPM'].append(get(a, 8, p, 'PM', 'lebesgue'))
            cols['lebPH'].append(get(a, 8, p, 'PH', 'lebesgue'))
        write_table(os.path.join(DATA, 'cyl-%s-N8-p%d.dat' % (problem, p)), cols, list(cols))
    with open(os.path.join(GEN, 'cyl-%s.tex' % problem), 'w') as f:
        ps = [4, 6, 8]
        f.write('\\begin{tabular}{@{}rl' + 'ccc' * len(ps) + '@{}}\n\\toprule\n')
        f.write('& & ' + ' & '.join('\\multicolumn{3}{c}{$p=%d$}' % p for p in ps) + '\\\\\n')
        f.write(''.join('\\cmidrule(lr){%d-%d}' % (3 + 3 * i, 5 + 3 * i) for i in range(len(ps))) + '\n')
        f.write('$\\alpha$ & & ' + ' & '.join(['$e$ & $r$ & $\\Lambda$'] * len(ps)) + '\\\\\n\\midrule\n')
        for a in alphas:
            for k, (fm, name) in enumerate((('PP', '\\PP'), ('PM', '\\PM'), ('PH', '\\PH'))):
                first = ('\\multirow{3}{*}{$%d^{\\circ%s}$}' % (a, '\\ast' if a >= 45 else '')) if k == 0 else ''
                vals = []
                key = 'energy' if a < 45 else 'nodal'      # alpha >= 45: folded at N = 8
                for p in ps:
                    e8, e4 = get(a, 8, p, fm, key), get(a, 4, p, fm, key)
                    r = math.log2(e4 / e8) if e4 > 0 and e8 > 0 else float('nan')
                    leb = get(a, 8, p, fm, 'lebesgue')
                    vals += [sci(e8), '--' if math.isnan(r) else num(r, '%.1f'),
                             '--' if fm == 'PP' or math.isnan(leb) else sci(leb, 1)]
                f.write('%s & %s & %s\\\\\n' % (first, name, ' & '.join(vals)))
            f.write('\\midrule\n' if a != alphas[-1] else '')
        f.write('\\bottomrule\n\\end{tabular}\n')


if __name__ == '__main__':
    tags = sys.argv[1:] or ['d010', 'd020', 'd030']
    for t in tags:
        convergence(t)
    patch()
    panel()
    growth()
    cylinder()
    print('done')
