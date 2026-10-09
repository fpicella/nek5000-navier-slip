#!/usr/bin/env python3
"""Polar O-grid spectral-element mesh around a circular cylinder -> Nek5000 .re2

Cylinder of radius a centred at the origin, outer circular boundary at r_inf.
ntheta uniform sectors, nr rings with geometric radial stretching starting
from a first-ring thickness dr1.  All ring sides are exact circular arcs
(Nek 'C' curved sides).  Boundary conditions written in the file:
  cylinder (r=a)         : 'W  '  (changed to 'shl' at run time by the .usr)
  outer circle, x<0 half : 'v  '  (uniform inflow)
  outer circle, x>0 half : 'O  '  (outflow)
Vertex ordering per element (counter-clockwise, as in Nek):
  v1=(r_k,th_i) v2=(r_k+1,th_i) v3=(r_k+1,th_i+1) v4=(r_k,th_i+1)
  side 2 = outer arc (radius +r_k+1), side 4 = inner arc (radius -r_k).
3D (--nz NZ >= 3, --lz LZ): the 2D mesh extruded in NZ uniform layers,
0 <= z <= LZ, periodic in z.  Hex vertices 1-4 = 2D element at z_k, 5-8 at
z_k+1; faces 1-4 keep the 2D sides' BCs, the arcs sit on edges 2/6 (outer)
and 4/8 (inner); faces 5 (bottom layer) and 6 (top layer) are 'P', connected
to the same 2D element in the top / bottom layer.
"""
import argparse, struct, numpy as np

def radii(a, rinf, nr, dr1):
    """ring radii with geometric growth q solved so that sum(dr) = rinf-a"""
    L = rinf - a
    f = lambda q: dr1*(q**nr - 1.0)/(q - 1.0) - L
    lo, hi = 1.0 + 1e-9, 4.0
    for _ in range(200):
        mid = 0.5*(lo + hi)
        if f(mid) > 0: hi = mid
        else: lo = mid
    q = 0.5*(lo + hi)
    dr = dr1*q**np.arange(nr)
    r = np.concatenate(([a], a + np.cumsum(dr)))
    r[-1] = rinf
    return r, q

def pack_char(s, n=8):
    b = s.encode('ascii'); b = b + b'\0'*(n - len(b))
    return np.frombuffer(b, dtype='<f8')[0]

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-o', '--out', default='cyl')
    p.add_argument('--a', type=float, default=0.5)
    p.add_argument('--rinf', type=float, default=40.0)
    p.add_argument('--ntheta', type=int, default=40)
    p.add_argument('--nr', type=int, default=26)
    p.add_argument('--dr1', type=float, default=0.08)
    p.add_argument('--nz', type=int, default=0, help='0: 2D; >= 3: extruded, periodic in z')
    p.add_argument('--lz', type=float, default=1.0, help='span (3D)')
    args = p.parse_args()
    a, rinf, nt, nr, dr1 = args.a, args.rinf, args.ntheta, args.nr, args.dr1
    assert nt % 4 == 0, 'ntheta must be a multiple of 4 (v/O split at x=0)'
    r, q = radii(a, rinf, nr, dr1)
    th = 2.0*np.pi*np.arange(nt)/nt
    X = np.outer(r, np.cos(th)); Y = np.outer(r, np.sin(th))   # exact closure in theta
    nel = nt*nr
    xc = np.zeros((nel, 4)); yc = np.zeros((nel, 4))
    curves = []; bcs = []
    e = 0
    for k in range(nr):            # rings, inner to outer
        for i in range(nt):        # sectors
            i1 = (i + 1) % nt
            xc[e] = [X[k, i], X[k+1, i], X[k+1, i1], X[k, i1]]
            yc[e] = [Y[k, i], Y[k+1, i], Y[k+1, i1], Y[k, i1]]
            curves.append((e+1, 2,  r[k+1]))
            curves.append((e+1, 4, -r[k]))
            if k == 0:
                bcs.append((e+1, 4, 'W  '))
            if k == nr - 1:
                thm = 0.5*(th[i] + th[i] + 2.0*np.pi/nt)
                bcs.append((e+1, 2, 'v  ' if np.cos(thm) < 0 else 'O  '))
            e += 1
    # sanity: counter-clockwise (positive shoelace area) for every element
    area = 0.5*np.sum(xc*np.roll(yc, -1, axis=1) - np.roll(xc, -1, axis=1)*yc, axis=1)
    assert (area > 0).all(), 'element orientation error'
    if args.nz > 0:
        write3d(args, nel, xc, yc, curves, bcs)
        return
    with open(args.out + '.re2', 'wb') as f:
        hdr = '#v002%9d%3d%9d this is the hdr' % (nel, 2, nel)
        f.write(hdr.ljust(80).encode('ascii'))
        f.write(struct.pack('<f', 6.54321))
        rec = np.zeros((nel, 9))
        rec[:, 1:5] = xc; rec[:, 5:9] = yc
        rec.astype('<f8').tofile(f)
        np.array([len(curves)], dtype='<f8').tofile(f)
        cr = np.zeros((len(curves), 8))
        for j, (ie, s, rad) in enumerate(curves):
            cr[j, 0] = ie; cr[j, 1] = s; cr[j, 2] = rad; cr[j, 7] = pack_char('C')
        cr.astype('<f8').tofile(f)
        np.array([len(bcs)], dtype='<f8').tofile(f)
        br = np.zeros((len(bcs), 8))
        for j, (ie, s, cb) in enumerate(bcs):
            br[j, 0] = ie; br[j, 1] = s; br[j, 7] = pack_char(cb)
        br.astype('<f8').tofile(f)
    from collections import Counter
    print('wrote %s.re2: nel=%d (ntheta=%d x nr=%d), a=%g rinf=%g, dr1=%g growth q=%.4f, last dr=%.3f, outer arc=%.3f'
          % (args.out, nel, nt, nr, a, rinf, dr1, q, r[-1]-r[-2], 2*np.pi*rinf/nt))
    print('curved sides', len(curves), ' BC faces', Counter(cb for _, _, cb in bcs))
    print('ring radii:', np.array2string(r, precision=3, max_line_width=120))

def write3d(args, nel, xc, yc, curves, bcs):
    nz, lz = args.nz, args.lz
    assert nz >= 3, 'nz >= 3 (periodic connectivity)'
    zl = lz*np.arange(nz + 1)/nz
    nel3 = nel*nz
    rec = np.zeros((nel3, 25))                     # igroup, x(8), y(8), z(8)
    c3, b3 = [], []
    for kz in range(nz):
        o = kz*nel
        rec[o:o+nel, 1:5] = xc; rec[o:o+nel, 5:9] = xc
        rec[o:o+nel, 9:13] = yc; rec[o:o+nel, 13:17] = yc
        rec[o:o+nel, 17:21] = zl[kz]; rec[o:o+nel, 21:25] = zl[kz+1]
        for ie, sd, rad in curves:                 # 2D side s -> edges s (bottom), s+4 (top)
            c3.append((o + ie, sd, rad)); c3.append((o + ie, sd + 4, rad))
        for ie, sd, cb in bcs:                     # side faces keep their numbers
            b3.append((o + ie, sd, cb, 0, 0))
    for e2 in range(1, nel + 1):                   # periodic z faces
        bot, top = e2, (nz - 1)*nel + e2
        b3.append((bot, 5, 'P  ', top, 6)); b3.append((top, 6, 'P  ', bot, 5))
    with open(args.out + '.re2', 'wb') as f:
        f.write(('#v002%9d%3d%9d this is the hdr' % (nel3, 3, nel3)).ljust(80).encode('ascii'))
        f.write(struct.pack('<f', 6.54321))
        rec.astype('<f8').tofile(f)
        np.array([len(c3)], dtype='<f8').tofile(f)
        cr = np.zeros((len(c3), 8))
        for j, (ie, ed, rad) in enumerate(c3):
            cr[j, 0] = ie; cr[j, 1] = ed; cr[j, 2] = rad; cr[j, 7] = pack_char('C')
        cr.astype('<f8').tofile(f)
        np.array([len(b3)], dtype='<f8').tofile(f)
        br = np.zeros((len(b3), 8))
        for j, (ie, fc, cb, p1, p2) in enumerate(b3):
            br[j, 0] = ie; br[j, 1] = fc; br[j, 2] = p1; br[j, 3] = p2; br[j, 7] = pack_char(cb)
        br.astype('<f8').tofile(f)
    from collections import Counter
    print('wrote %s.re2 (3D): nel=%d = %d x nz=%d, span lz=%g, dz=%g' % (args.out, nel3, nel, nz, lz, lz/nz))
    print('curved edges', len(c3), ' BC faces', Counter(b[2] for b in b3))

if __name__ == '__main__':
    main()
