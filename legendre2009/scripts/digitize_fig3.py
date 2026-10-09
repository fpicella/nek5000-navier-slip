#!/usr/bin/env python3
"""Digitise Legendre, Lauga & Magnaudet (2009) fig. 3 (St* = St(Kn)/St(0) vs Kn/Knc,
Re = 50, 100, 200, 500, 800) from the vector data of the arXiv PDF: the polylines
joining the markers are the data.  Page 7: x = 218.20 + 161.30*Kn/Knc (ticks every
0.05 = 8.065 pt), St* = 1 + (239.57 - y)/38.39 (ticks every 0.5 = 19.195 pt).
Series identified by their end value (Re = 800 highest).
Output: legendre_fig3_all.dat (Re Kn/Knc St*)"""
import sys, math, numpy as np
import pymupdf as fitz

pdf = sys.argv[1] if len(sys.argv) > 1 else 'Legendre_JFM_2009.pdf'
dr = fitz.open(pdf)[6].get_drawings()
X = lambda x: (x - 218.20)/161.30
Y = lambda y: 1.0 + (239.57 - y)/38.39
lines = []
for k in (5, 6):                            # data drawings: polylines + markers
    cur = []
    for it in dr[k]['items']:
        if it[0] != 'l': continue
        a, b = it[1], it[2]
        if math.hypot(b.x - a.x, b.y - a.y) <= 3.5: continue          # marker strokes
        if cur and (abs(cur[-1][0] - a.x) > 0.05 or abs(cur[-1][1] - a.y) > 0.05):
            lines.append(cur); cur = []
        if not cur: cur = [(a.x, a.y)]
        cur.append((b.x, b.y))
    if cur: lines.append(cur)
lines = [l for l in lines if len(l) >= 6]
lines.sort(key=lambda l: -l[-1][1])         # increasing end value of St*
assert len(lines) == 5, len(lines)
with open('legendre_fig3_all.dat', 'w') as f:
    f.write('# Re Kn/Knc St*   Legendre et al. 2009 fig 3 (vector data of the arXiv PDF); their Knc = 0.09, 0.38, 0.47, 0.49, 0.53\n')
    for re, l in zip((50, 100, 200, 500, 800), lines):
        for x, y in l: f.write('%5d %8.4f %8.4f\n' % (re, X(x), Y(y)))
        print('Re %3d: %2d points, Kn/Knc %.3f..%.3f, St* %.3f..%.3f' % (re, len(l), X(l[0][0]), X(l[-1][0]), Y(l[0][1]), Y(l[-1][1])))
