#!/usr/bin/env python3
"""Digitise Legendre, Lauga & Magnaudet (2009) fig. 4 (omega_max/(U/a) at the
cylinder surface vs Re, log-log) from the vector data of the arXiv PDF (page 8).
Calibration from the tick marks: Re = 10*10**((x - 251.3)/76.7) (the axis spans
Re = 2..1000; decades at x = 251.3 and 328.0), w = 10**((293.4 - y)/107.7)
(ticks 1 at y = 293.4 and 10 at y = 185.7).  Curves: drawing 76 thick solid
Kn = 0, 78 dash-dotted Kn = inf, 77 + 79 thin lines Kn = 0.05, 0.1, 0.2, 0.4, 1,
2.2, 5 (ordered by value), 74 long dashes: separation, 75 short dashes: shedding.
Output: legendre_fig4_curves.dat (curve Re w)"""
import sys, numpy as np
import pymupdf as fitz

pdf = sys.argv[1] if len(sys.argv) > 1 else 'Legendre_JFM_2009.pdf'
dr = fitz.open(pdf)[7].get_drawings()
RE = lambda x: 10.0*10**((x - 251.3)/76.7)
W = lambda y: 10**((293.4 - y)/107.7)

def polylines(k):
    out, cur = [], []
    for it in dr[k]['items']:
        if it[0] != 'l': continue
        a, b = it[1], it[2]
        if cur and (abs(cur[-1][0] - a.x) > 0.05 or abs(cur[-1][1] - a.y) > 0.05):
            out.append(cur); cur = []
        if not cur: cur = [(a.x, a.y)]
        cur.append((b.x, b.y))
    if cur: out.append(cur)
    return [np.array([(RE(x), W(y)) for x, y in l]) for l in out]

curves = {'Kn0': polylines(76)[0], 'Kninf': polylines(78)[0], 'sep': polylines(74)[0], 'shed': polylines(75)[0]}
thin = polylines(77) + polylines(79)
assert len(thin) == 7, len(thin)
thin.sort(key=lambda c: -np.interp(np.log(100.0), np.log(c[:, 0]), c[:, 1]))   # by value at Re = 100
for kn, c in zip(('0.05', '0.1', '0.2', '0.4', '1', '2.2', '5'), thin): curves['Kn' + kn] = c
with open('legendre_fig4_curves.dat', 'w') as f:
    f.write('# curve Re omega_max*a/U   Legendre et al. 2009 fig 4 (vector data of the arXiv PDF, page 8)\n')
    for name, c in curves.items():
        for re, w in c: f.write('%-6s %9.3f %8.4f\n' % (name, re, w))
        print('%-6s %2d points  Re %6.1f..%6.1f  w %.3f..%.3f  (w at Re 100: %s)' % (name, len(c), c[0, 0], c[-1, 0], c[0, 1], c[-1, 1],
              '%.3f' % np.interp(np.log(100.0), np.log(c[:, 0]), c[:, 1]) if c[0, 0] <= 100 <= c[-1, 0] else '-'))
