#!/usr/bin/env python3
"""Table of the calibrated Gauss rules (manuscript, Sect. 5).

Reads data/quadrature/<family>.dat (ms_quadrature_calibration.m) and writes
  generated/tab_quadrule.tex   rows: mesh family, columns: p = 2..8, entries
                               nG of SEM / USEM / USEM-H ('n.c.': no rule
                               below the reference rule converges), and the
                               largest relative deviation of the former 2p
                               rule from the reference over p = 2..8 and
                               the formulations with converged rules
  data/quadrature/summary.txt  the same numbers in plain text
For the thick cylinder, the rule of an angle is the maximum over the meshes
N = 1, 2, 4, 8 (the study uses the rule of each mesh).
"""
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
QDIR = os.path.join(HERE, '..', 'data', 'quadrature')
GEN = os.path.join(HERE, '..', 'generated')
FORMS = ('PP', 'PM', 'PH')
PS = range(2, 9)


def read(family):
    path = os.path.join(QDIR, family + '.dat')
    with open(path) as f:
        head = f.readline().split()
        rows = [[float(v) for v in line.split()] for line in f if line.strip()]
    return {int(r[0]): dict(zip(head, r)) for r in rows}


def merge(families):
    """Maximum rule over several families (meshes) per p and form."""
    out = {}
    for fam in families:
        for p, r in read(fam).items():
            o = out.setdefault(p, {})
            for f in FORMS:
                if 'nG_' + f not in r:
                    continue
                nG, conv, d2p, rej = r['nG_' + f], r['conv_' + f], r['dev2p_' + f], r['rej_' + f]
                n0, c0, d0, r0 = o.get(f, (0, 1, 0, 0))
                if math.isnan(nG):
                    o[f] = (n0, c0, d0, 1)
                    continue
                o[f] = (max(n0, nG), min(c0, conv), max(d0, d2p), max(r0, rej))
    return out


def sci(v):
    e = int(math.floor(math.log10(v)))
    m = v / 10 ** e
    if round(m, 0) >= 10:
        m, e = m / 10, e + 1
    return '$%.0f{\\cdot}10^{%d}$' % (m, e)


ROWS = [
    ('square, $\\delta=0.1$', ['conv010']),
    ('square, $\\delta=0.2$', ['conv020']),
    ('square, $\\delta=0.3^\\ast$', ['conv030']),
    ('$3\\times2$, curved $0.15$', ['curved15']),
    ('$3\\times2$, curved $0.20$', ['curved20']),
    ('$3\\times2$, angular + curved $0.30$', ['combined30']),
] + [('cylinder, $\\alpha=%d^{\\circ%s}$' % (a, '\\ast' if a >= 45 else ''),
      ['cyl_a%03d_N%d' % (a, N) for N in (1, 2, 4, 8)]) for a in (0, 15, 30, 45, 60, 90)]


def cell(o, p):
    if p not in o:
        return '--'
    parts = []
    for f in FORMS:
        n, c, _, r = o[p].get(f, (math.nan, 0, 0, 0))
        if r and n == 0:
            parts.append('rej.')
        elif not c:
            parts.append('n.c.')
        else:
            parts.append('%d' % n)
    return '/'.join(parts)


def main():
    lines, txt = [], []
    for label, fams in ROWS:
        if not all(os.path.exists(os.path.join(QDIR, f + '.dat')) for f in fams):
            continue
        o = merge(fams)
        d2p = [o[p][f][2] for p in PS if p in o for f in FORMS
               if f in o[p] and o[p][f][1] and o[p][f][0] > 0]
        dmax = max(d2p) if d2p else math.nan
        lines.append(label + ' & ' + ' & '.join(cell(o, p) for p in PS) +
                     ' & ' + (sci(dmax) if d2p else '--') + '\\\\')
        txt.append('%-40s %s  max dev(2p) = %.2e' % (fams[0], ' '.join(cell(o, p) for p in PS), dmax))
        # deviation of the 2p rule per form (for the text)
        for f in FORMS:
            v = [o[p][f][2] for p in PS if p in o and f in o[p] and o[p][f][1] and o[p][f][0] > 0]
            if v:
                txt.append('    %s: dev(2p) %.2e .. %.2e' % (f, min(v), max(v)))
    head = ['\\begin{tabular}{@{}l' + 'c' * (len(PS) + 1) + '@{}}', '\\toprule',
            'mesh & ' + ' & '.join('$p=%d$' % p for p in PS) + ' & dev.\\ $2p$\\\\', '\\midrule']
    with open(os.path.join(GEN, 'tab_quadrule.tex'), 'w', newline='\n') as f:
        f.write('\n'.join(head + lines + ['\\bottomrule', '\\end{tabular}']) + '\n')
    with open(os.path.join(QDIR, 'summary.txt'), 'w', newline='\n') as f:
        f.write('\n'.join(txt) + '\n')
    print('\n'.join(txt))


if __name__ == '__main__':
    main()
