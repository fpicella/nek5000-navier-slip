#!/usr/bin/env python3
"""Collect the steady drag of the Re=20 Kn sweep (runs/re20_kn*/cyl_drag.dat)
and compare the normalised drag Cd*(Kn) = (Cd(Kn)-Cd(inf))/(Cd(0)-Cd(inf))
with Legendre, Lauga & Magnaudet (2009) figure 2a (Re=20 circles, digitised
from the arXiv PDF: legendre_fig2a_re20.dat).
Usage: plot_fig2a.py RUNS_DIR  -> results/re20_sweep.dat, results/fig2a_re20.png
"""
import sys, glob, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

runs = sys.argv[1] if len(sys.argv) > 1 else '../runs'

def load(prefix):
    out = []
    for d in sorted(glob.glob(os.path.join(runs, prefix + '_kn*'))):
        f = os.path.join(d, 'cyl_drag.dat')
        if not os.path.exists(f): continue
        kn = float(d.split(prefix + '_kn')[-1])
        a = np.loadtxt(f)
        if a.ndim == 1 or len(a) < 2: continue
        t, cd, cdp, cdv, cl, usmax, usmin, res, utmax, dudt = a[-1, 1:11]
        m = a[:, 1] >= t - 10.0          # Cd drift over the last 10 time units
        drift = a[m, 2].max() - a[m, 2].min()
        out.append((kn, t, cd, cdp, cdv, usmax, usmin, res, dudt, drift))
    out = np.array(out)
    return out[np.argsort(out[:, 0])] if len(out) else out   # numeric Kn order

rows = load('re20')
fine = load('fine')
hdr = 'Kn(<0=shear-free) t_end Cd Cd_p Cd_v us_max us_min res_navier dudt_max Cd_drift_last10'
np.savetxt('re20_sweep.dat', rows, header=hdr, fmt='%12.5e')
print(hdr); print(rows)

cd0 = rows[rows[:, 0] == 0.0, 2]; cdi = rows[rows[:, 0] < 0.0, 2]
if len(cd0) and len(cdi):
    cd0, cdi = cd0[0], cdi[0]
    print('Cd(Kn=0)   = %.4f   (Legendre 2.035, Dennis & Chang 2.045)' % cd0)
    print('Cd(Kn=inf) = %.4f   (Legendre 1.33)' % cdi)
    sl = rows[rows[:, 0] > 0.0]
    cds = (sl[:, 2] - cdi)/(cd0 - cdi)
    ref = np.loadtxt('legendre_fig2a_re20.dat')
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].semilogx(ref[:, 0], ref[:, 1], 'o', mfc='none', ms=8, color='k', label='Legendre et al. 2009, fig. 2a (Re=20)')
    ax[0].semilogx(sl[:, 0], cds, 's', ms=5, color='C3', label='Nek5000, implicit Robin (this work)')
    ax[0].set_xlabel('Kn = $\\lambda/a$'); ax[0].set_ylabel('$C_D^* = (C_D(Kn)-C_D(\\infty))/(C_D(0)-C_D(\\infty))$')
    ax[0].set_xlim(0.008, 30); ax[0].set_ylim(0, 1); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)
    ax[0].set_title('Re = 20: normalised drag vs slip')
    ax[1].semilogx(sl[:, 0], sl[:, 2], 's-', color='C3', label='$C_D$ total')
    ax[1].semilogx(sl[:, 0], sl[:, 3], '^--', color='C0', label='pressure')
    ax[1].semilogx(sl[:, 0], sl[:, 4], 'v--', color='C2', label='viscous')
    ax[1].axhline(cd0, color='k', ls=':', label='no-slip %.3f' % cd0)
    ax[1].axhline(cdi, color='gray', ls=':', label='shear-free %.3f' % cdi)
    ax[1].set_xlabel('Kn'); ax[1].set_ylabel('$C_D$'); ax[1].grid(alpha=.3); ax[1].legend(fontsize=8)
    ax[1].set_title('Re = 20: drag decomposition')
    if len(fine):
        f0 = fine[fine[:, 0] == 0.0, 2]; fi = fine[fine[:, 0] < 0.0, 2]
        if len(f0) and len(fi):
            fs = fine[fine[:, 0] > 0.0]
            ax[0].semilogx(fs[:, 0], (fs[:, 2] - fi[0])/(f0[0] - fi[0]), 'x', ms=9, color='C0',
                           label='Nek5000, refined mesh (1904 el.)')
            ax[0].legend(fontsize=8)
        print('\nmesh refinement (1040 el., dr1=0.08  vs  1904 el., dr1=0.05), N=7:')
        print('   Kn      Cd_1040     Cd_1904     rel.diff')
        for r in fine:
            c = rows[rows[:, 0] == r[0], 2]
            if len(c):
                print('%7.3f  %10.6f  %10.6f  %+9.2e' % (r[0], c[0], r[2], (c[0] - r[2])/r[2]))
        np.savetxt('re20_sweep_fine.dat', fine, header=hdr, fmt='%12.5e')
    fig.tight_layout(); fig.savefig('fig2a_re20.png', dpi=150)
    # interpolate reference at our Kn for a numeric comparison
    refi = np.interp(np.log(sl[:, 0]), np.log(ref[:, 0]), ref[:, 1])
    print('\n   Kn      Cd      Cd*_nek   Cd*_Legendre   diff')
    for k, c, cs, r in zip(sl[:, 0], sl[:, 2], cds, refi):
        print('%7.3f  %7.4f  %8.4f  %8.4f  %+8.4f' % (k, c, cs, r, cs - r))
