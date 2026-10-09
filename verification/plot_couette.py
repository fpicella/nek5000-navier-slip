#!/usr/bin/env python3
"""plot_couette.py [RUNS_DIR] -- the circular-Couette validation in one figure
(RUNS_DIR/c3_*/profile.dat and p3_*/verif_summary.dat, default data/).
Inner cylinder r = a = 0.5 at rest with the Navier slip condition, outer
cylinder r = R = 1.5 rotating at omega = 1 (no slip).  Exact solution
u_th = A r + B/r with u_th = lambda (du_th/dr - u_th/r) at r = a.
(a) azimuthal velocity u_th(r): Nek5000 at all GLL points against the exact solution
(b) slip velocity on the inner cylinder versus Kn, Kn = 1e-3 to 1e3
(c) relative error of that slip velocity, Nek5000 at N = 7
(d) Kn = 1: error against the exact solution versus polynomial order"""
import glob, os, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
runs = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'data')
a, R, om = 0.5, 1.5, 1.0

def coef(kn, flat=False):
    if kn == 0:
        B = -om*a/(1/a - a/R**2)
    elif kn < 0:
        B = om/(1/a**2 + 1/R**2) if flat else 0.0
    else:
        lam = kn*a
        B = (om*(lam - a)/(1/a - a/R**2 + lam/R**2 + lam/a**2) if flat
             else -om*a/(1/a - a/R**2 + 2*lam/a**2))
    return om - B/R**2, B

def uth(kn, r, flat=False):
    A, B = coef(kn, flat); return A*r + B/r

def nek(kn):
    f = os.path.join(runs, 'c3_kn%g_N7' % kn, 'profile.dat')
    if not os.path.exists(f): return None
    x, y, ux, uy, uxe, uye = np.loadtxt(f, unpack=True)
    r = np.hypot(x, y); ut = (-ux*y + uy*x)/r
    return r, ut

lab = lambda kn: 'no-slip' if kn == 0 else ('shear-free' if kn < 0 else 'Kn = %g' % kn)
fig, ax = plt.subplots(2, 2, figsize=(12, 9)); ax = ax.ravel()
rr = np.linspace(a, R, 200)
for i, kn in enumerate((0.0, 0.1, 1.0, -1.0)):
    c = 'C%d' % i
    ax[0].plot(rr, uth(kn, rr), '-', color=c, lw=1.5, label=lab(kn))
    d = nek(kn)
    if d is not None:
        ax[0].plot(d[0][::3], d[1][::3], 'o', color=c, ms=3.5, mfc='none')
ax[0].plot([], [], 'k-', label='exact solution')
ax[0].plot([], [], 'ko', mfc='none', ms=4, label='Nek5000, N = 7')
ax[0].set(xlabel='r   (inner cylinder r = 0.5, outer cylinder r = 1.5)', ylabel=r'azimuthal velocity $u_\theta$',
          title='(a) velocity across the gap', xlim=(a, R))
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)

K = np.logspace(-3, 3, 300)
ax[1].loglog(K, [uth(k, a) for k in K], 'k-', label='exact solution')
ax[1].axhline(uth(-1, a), color='C3', ls=':', label=r'shear-free limit: rigid rotation, $\omega a$')
print('   Kn     u_wall Nek       exact            rel. error')
kk, uw = [], []
for kn in (0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0):
    d = nek(kn)
    if d is None: continue
    w = np.abs(d[0] - a) < 1e-9; u = d[1][w].mean()
    kk.append(kn); uw.append(u)
    print('%7g  %.10f  %.10f  %.2e' % (kn, u, uth(kn, a), abs(u - uth(kn, a))/uth(kn, a)))
ax[1].loglog(kk, uw, 'o', color='C0', ms=7, mfc='none', mew=1.5, label='Nek5000, N = 7')
ax[1].set(xlabel=r'Kn = $\lambda/a$', ylabel=r'slip velocity $u_\theta(a)$', title='(b) slip velocity on the fixed inner cylinder',
          xlim=(5e-4, 2e3), ylim=(5e-4, 1))
ax[1].legend(fontsize=8, loc='lower right'); ax[1].grid(alpha=.3, which='both')
rel = [abs(u - uth(k, a))/uth(k, a) for k, u in zip(kk, uw)]
ax[2].loglog(kk, rel, 'o-', color='C0', ms=7, mfc='none', mew=1.5, label='Nek5000, N = 7, 8 x 2 elements')
ax[2].set(xlabel=r'Kn = $\lambda/a$', ylabel=r'$|u_\theta(a) - u_{exact}| \, / \, u_{exact}$', title='(c) relative error of the slip velocity',
          xlim=(5e-4, 2e3), ylim=(1e-12, 1e-6))
ax[2].legend(fontsize=8); ax[2].grid(alpha=.3, which='both')

N, e, q = [], [], []
for f in sorted(glob.glob(os.path.join(runs, 'p3_kn1_N*', 'verif_summary.dat'))):
    s = np.atleast_2d(np.loadtxt(f))[0]; N.append(s[1]); e.append(s[6]); q.append(s[9])
o = np.argsort(N); N, e, q = np.array(N)[o], np.array(e)[o], np.array(q)[o]
ax[3].semilogy(N, e, 'ks-', mfc='none', label=r'max $|u - u_{exact}|$ over the domain')
ax[3].semilogy(N, q, 'k^--', mfc='none', label='relative error of the torque')
ax[3].set(xlabel='polynomial order N', ylabel='error', title='(d) Kn = 1: error versus resolution')
ax[3].legend(fontsize=8); ax[3].grid(alpha=.3, which='both')
fig.suptitle('Circular Couette flow: inner cylinder (r = 0.5) at rest with Navier slip,\nouter cylinder (r = 1.5) rotating at $\\omega$ = 1', fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'couette_validation.png'), dpi=140)
