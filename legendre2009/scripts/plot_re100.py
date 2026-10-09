#!/usr/bin/env python3
"""plot_re100.py [RUNS_DATA] -- Re = 100 with slip: mean drag, lift amplitude
and Strouhal number from runs_data/re100_kn*/cyl_drag.dat (last 100 time
units), compared with Legendre et al. (2009) fig 2a, 2b and 3 (Re = 100
diamonds, Knc = 0.38).  Writes re100_sweep.dat and re100_vs_legendre.png."""
import glob, os, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # legendre2009/
DATA, FIGS = os.path.join(HERE, 'data'), os.path.join(HERE, 'figures')
runs = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'runs')
KNC = 0.38

def stats(f, twin=100.0):
    a = np.loadtxt(f); t, cd, cl = a[:, 1], a[:, 2], a[:, 5]
    m = t >= t[-1] - twin
    t, cd, cl = t[m], cd[m], cl[m]
    amp = 0.5*(cl.max() - cl.min())
    c = cl - cl.mean(); up = np.where((c[:-1] < 0) & (c[1:] >= 0))[0]
    tz = t[up] - c[up]*(t[up + 1] - t[up])/(c[up + 1] - c[up])
    st = 1.0/np.mean(np.diff(tz)) if len(tz) > 2 and amp > 1e-3 else 0.0
    h = len(cl)//2        # growth check: lift amplitude in the two halves of the window
    g = (cl[h:].max() - cl[h:].min())/max(cl[:h].max() - cl[:h].min(), 1e-30)
    return cd.mean(), cl.max() if amp > 1e-3 else 0.0, st, amp, g, t[-1]

def growth(f, tmin=40.0):
    # linear growth (or decay) rate of the lift envelope, fitted on the lift
    # maxima that are still in the linear regime (below 30% of the final amplitude)
    a = np.loadtxt(f); t, cl = a[:, 1], a[:, 5]
    m = t > tmin; t, cl = t[m], cl[m]
    i = np.where((cl[1:-1] > cl[:-2]) & (cl[1:-1] >= cl[2:]))[0] + 1
    tp, ap = t[i], cl[i]
    amax = np.abs(cl[t >= t[-1] - 50]).max()
    if amax > 0.01:        # saturated: fit the early linear phase only
        sel = (ap > 1e-9) & (ap < 0.3*amax)
    else:                  # still linear at the end (near threshold): fit the
        sel = (ap > 1e-9) & (tp > 0.5*t[-1])   # second half, transients gone
    if sel.sum() < 4:
        return np.nan
    return np.polyfit(tp[sel], np.log(ap[sel]), 1)[0]

rows = {}
for pre in ('re100', 're100fine'):
    for d in glob.glob(os.path.join(runs, pre + '_kn*')):
        f = os.path.join(d, 'cyl_drag.dat')
        if os.path.exists(f):
            rows[(pre, float(d.split('_kn')[-1]))] = stats(f)
R = {k[1]: v for k, v in rows.items() if k[0] == 're100'}
kns = sorted(R)
print('   Kn     t_end   Cd_mean   CL_max    CL_amp   St       amp ratio (2nd/1st half)')
for k in kns:
    cd, clm, st, amp, g, te = R[k]
    print('%7.3f  %6.1f  %.4f   %.4f   %.4f  %.4f   %.4f' % (k, te, cd, clm, amp, st, g))
for k, v in sorted(rows.items()):
    if k[0] == 're100fine':
        print('fine mesh Kn=%g: Cd %.4f  CL_max %.4f  St %.4f   (production mesh: %.4f  %.4f  %.4f)'
              % (k[1], v[0], v[1], v[2], R[k[1]][0], R[k[1]][1], R[k[1]][2]))
G = {}
for k in kns:
    if k >= 0.2:
        G[k] = growth(os.path.join(runs, 're100_kn%g' % k, 'cyl_drag.dat'))
print('\nlift-envelope growth rate (linear regime):')
for k in sorted(G):
    print('   Kn = %5.3f   growth rate = %+.5f' % (k, G[k]))
gk = sorted(k for k in G if np.isfinite(G[k]))
knc = np.nan
for k1, k2 in zip(gk[:-1], gk[1:]):
    if G[k1] > 0 >= G[k2]:
        knc = k1 + (k2 - k1)*G[k1]/(G[k1] - G[k2])
print('   threshold (zero growth, linear interpolation): Knc = %.3f   (Legendre et al.: 0.38)' % knc)
np.savetxt(os.path.join(DATA, 're100_sweep.dat'), np.array([[k] + list(R[k]) for k in kns]),
           header='Kn Cd_mean CL_max St CL_amp amp_growth_ratio t_end', fmt='%12.6f')
cd0, cl0, st0 = R[0.0][0], R[0.0][1], R[0.0][2]
cdi = R[-1.0][0] if -1.0 in R else np.nan
print('\nno-slip: Cd %.3f  CL %.3f  St %.4f | Legendre 1.350, 0.334, 0.176 | lit. Cd 1.33-1.35, CL ~0.33, St 0.164-0.167' % (cd0, cl0, st0))
print('shear-free: Cd %.4f | Legendre 0.415' % cdi)
ks = [k for k in kns if k > 0]
fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
L = np.loadtxt(os.path.join(DATA, 'legendre_re100_fig2a.dat'))
ax[0].semilogx(L[:, 0], L[:, 1], 'D', mfc='none', ms=7, color='k', label='Legendre et al. 2009 fig 2a')
ax[0].semilogx(ks, [(R[k][0] - cdi)/(cd0 - cdi) for k in ks], 's', color='C3', label='Nek5000')
ax[0].set(xlabel='Kn', ylabel='$C_D^*$', title='(a) Re = 100: normalised mean drag', xlim=(0.008, 30), ylim=(0, 1))
L = np.loadtxt(os.path.join(DATA, 'legendre_re100_fig2b.dat'))
ax[1].plot(L[:, 0], L[:, 1], 'D', mfc='none', ms=7, color='k', label='Legendre fig 2b')
ax[1].plot([k/KNC for k in ks], [R[k][1]/cl0 for k in ks], 's', mfc='none', color='C3', label='Nek5000, Kn$_c$ = 0.38 (Legendre)')
if np.isfinite(knc):
    ax[1].plot([k/knc for k in ks], [R[k][1]/cl0 for k in ks], 's', color='C3', label='Nek5000, own Kn$_c$ = %.3f' % knc)
ax[1].set(xlabel='Kn / Kn$_c$', ylabel='$C_L^* = C_L/C_L(0)$', title='(b) lift amplitude', xlim=(-0.02, 1.4), ylim=(-0.02, 1.05))
L = np.loadtxt(os.path.join(DATA, 'legendre_re100_fig3.dat'))
ax[2].plot(L[:, 0], L[:, 1], 'D', mfc='none', ms=7, color='k', label='Legendre fig 3')
kk = [k for k in ks if R[k][2] > 0]
ax[2].plot([k/KNC for k in kk], [R[k][2]/st0 for k in kk], 's', mfc='none', color='C3', label='Nek5000, Kn$_c$ = 0.38')
if np.isfinite(knc):
    ax[2].plot([k/knc for k in kk], [R[k][2]/st0 for k in kk], 's', color='C3', label='Nek5000, own Kn$_c$')
ax[2].set(xlabel='Kn / Kn$_c$', ylabel='$St^* = St/St(0)$', title='(c) Strouhal number', xlim=(-0.02, 1.05))
for a in ax: a.grid(alpha=.3); a.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(FIGS, 're100_vs_legendre.png'), dpi=140)
