#!/usr/bin/env python3
"""plot_legendre.py [RUNS_DATA] -- Legendre, Lauga & Magnaudet (2009)
figures 1 and 2 from the Nek5000 runs (runs_data/{re20,re100,L<Re>,Lfine<Re>}_kn*).
  legendre_fig1.png         stability diagram: separation and shedding boundaries
  legendre_fig1_insets.png  vorticity snapshots (Re = 200 row, Kn = 0.2 column)
  legendre_fig2.png         Cd*(Kn) for Re = 20-800 and CL*(Kn/Knc) for Re = 50-800
  legendre_campaign.txt     per-run statistics and the derived boundaries
Separation: sign of the slope cw of the wall vorticity at the rear stagnation
point (steady runs).  Shedding: sign of the growth rate of the lift envelope."""
import glob, os, re, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
runs = sys.argv[1] if len(sys.argv) > 1 else '../runs_data'
OUT = []
def P(s=''):
    print(s); OUT.append(s)

# ---------------------------------------------------------------- data
def growth(t, cl, tmin=40.0):
    m = t > tmin; t, cl = t[m], cl[m]
    if len(t) < 20: return np.nan
    i = np.where((cl[1:-1] > cl[:-2]) & (cl[1:-1] >= cl[2:]))[0] + 1
    tp, ap = t[i], cl[i]
    amax = np.abs(cl[t >= t[-1] - 50]).max()
    sel = ((ap > 1e-9) & (ap < 0.3*amax)) if amax > 0.01 else ((ap > 1e-9) & (tp > 0.5*t[-1]))
    return np.polyfit(tp[sel], np.log(ap[sel]), 1)[0] if sel.sum() >= 4 else np.nan

def stats(a):
    a = a[a[:, 1] <= a[-1, 1] - 0.5]          # drop the last 0.5 time units (blow-up tail of a failed run)
    t, cd, cl = a[:, 1], a[:, 2], a[:, 5]
    win = t >= t[-1] - min(100.0, t[-1]/3)
    clw = cl[win]; amp = 0.5*(clw.max() - clw.min())
    r = dict(t=t[-1], growth=growth(t, cl), unsteady=amp > 1e-3, amp=amp,
             cw=a[-1, 12] if a.shape[1] > 13 else np.nan, ths=a[-1, 13] if a.shape[1] > 13 else np.nan)
    if r['unsteady']:
        c = clw - clw.mean(); tw = t[win]
        up = np.where((c[:-1] < 0) & (c[1:] >= 0))[0]
        tz = tw[up] - c[up]*(tw[up + 1] - tw[up])/(c[up + 1] - c[up])
        h = len(clw)//2
        r.update(cd=cd[win].mean(), cl=clw.max(), st=1/np.mean(np.diff(tz)) if len(tz) > 2 else np.nan,
                 sat=(clw[h:].max() - clw[h:].min())/max(clw[:h].max() - clw[:h].min(), 1e-30))
    else:
        r.update(cd=cd[-1], cl=0.0, st=np.nan, sat=np.nan)
    return r

def wallsep(d):
    """local slope of the wall vorticity at the rear stagnation point (fit
    over |th| < 5 deg) and upper-surface separation angle (zero crossing of
    the wall vorticity, linear interpolation), from wall_*.dat"""
    fs = glob.glob(os.path.join(d, 'diag', 'wall_*.dat')) or glob.glob(os.path.join(d, 'wall_*.dat'))
    fs = [f for f in fs if os.path.getsize(f) > 0]          # ranks without wall nodes write empty files
    if not fs: return np.nan, np.nan
    W = np.vstack([np.loadtxt(f, ndmin=2) for f in fs]); W = W[np.argsort(W[:, 0])]
    th, om = W[:, 0], W[:, 2]
    k = np.abs(th) < np.radians(5.0)
    A = np.c_[th[k], th[k]**3]; c = np.linalg.lstsq(A, om[k], rcond=None)[0][0]
    u = (th > 0) & (th < np.pi/2); tu, ou = th[u], om[u]
    ths = 0.0
    if ou[0] > 0:
        j = np.where((ou[:-1] > 0) & (ou[1:] <= 0))[0]
        if len(j): j = j[0]; ths = np.degrees(tu[j] + (tu[j + 1] - tu[j])*ou[j]/(ou[j] - ou[j + 1]))
    return c, ths

D = {}                                   # (fine, Re, Kn) -> stats ; L runs override re20/re100
H = {}                                   # hysteresis runs (restart from shedding at Kn = 0.3)
for d in sorted(glob.glob(os.path.join(runs, '*_kn*')), key=lambda s: os.path.basename(s).startswith('L')):
    m = re.match(r'^(re20|re100|Lfine|L|H)([0-9.]*)_kn(-?[0-9.]+)$', os.path.basename(d))
    f = os.path.join(d, 'cyl_drag.dat')
    if not m or not os.path.exists(f): continue
    fam, reS, knS = m.groups()
    Re = 20.0 if fam == 're20' else 100.0 if fam == 're100' else float(reS)
    a = np.loadtxt(f, ndmin=2)
    if len(a) <= 5: continue
    st = stats(a); st['wc'], st['wths'] = wallsep(d)
    if fam == 'H': H[(Re, float(knS))] = st
    else: D[(fam == 'Lfine', Re, float(knS))] = st
S = {(Re, kn): v for (fine, Re, kn), v in D.items() if not fine}
P('   Re      Kn      t_end  unsteady  Cd        CL_max   St       growth     cw(20deg) th_sep  | wall: slope(5deg) th_sep')
for (Re, kn) in sorted(S):
    v = S[(Re, kn)]
    P('%6g  %7g  %6.0f  %-5s  %8.5f  %7.4f  %7.4f  %+9.5f  %+8.4f  %5.1f  |  %+9.4f  %6.2f' % (Re, kn, v['t'], v['unsteady'],
      v['cd'], v['cl'], v['st'], v['growth'], v['cw'], v['ths'], v['wc'], v['wths']))

# ---------------------------------------------------------------- boundaries
def zero(x1, y1, x2, y2): return x1 + (x2 - x1)*y1/(y1 - y2)
sep = []; P('\nseparation boundary (steady runs with a wall distribution; zero of the slope of the wall')
P('vorticity at the rear stagnation point, fit over |th| < 5 deg; check: extrapolation of th_sep^2):')
for Re in sorted({r for r, k in S}):
    pts = sorted((k, S[(Re, k)]['wc'], S[(Re, k)]['wths']) for k in [k for r, k in S if r == Re]
                 if k > 0 and np.isfinite(S[(Re, k)]['wc']))
    for i in range(len(pts) - 1):
        (k1, c1, t1), (k2, c2, t2) = pts[i], pts[i + 1]
        if c1 > 0 >= c2:
            ks = zero(k1, c1, k2, c2)
            if (k2 - k1)/k1 > 0.5:           # bracket too wide for a linear interpolation
                P('   Re = %5g   Kn_sep in [%g, %g] (bracket too wide, not plotted)' % (Re, k1, k2)); continue
            sep.append((ks, Re))
            th2 = [(k, t**2) for k, c, t in pts[:i + 1] if t > 0][-2:]
            ext = np.nan                     # th_sep^2 is linear in Kn close to the onset
            if len(th2) == 2 and th2[0][1] > th2[1][1]:
                ext = th2[1][0] + th2[1][1]*(th2[1][0] - th2[0][0])/(th2[0][1] - th2[1][1])
            P('   Re = %5g   Kn_sep = %.3f   (slope %+.4f at Kn %g, %+.4f at Kn %g; th_sep^2 extrapolation %.3f)'
              % (Re, ks, c1, k1, c2, k2, ext))
k0 = sorted((Re, v['cw']) for (Re, k), v in S.items() if k == 0 and Re < 10 and np.isfinite(v['cw']))
for (r1, c1), (r2, c2) in zip(k0[:-1], k0[1:]):
    if c1 <= 0 < c2:
        re1 = zero(r1, c1, r2, c2); sep.append((0.0, re1)); P('   Kn = 0      Re1 = %.2f   (Legendre 6.5 in the text, 6.24 on the figure; lit. 6.2-6.3)' % re1)
shed = []; shed_err = []
P('\nshedding boundary: bracket [last shedding Kn, first decaying Kn]; estimate = zero of the lift growth rate')
P('when a growing rate is measured, otherwise zero of the squared saturated lift amplitude (supercritical onset,')
P('confirmed by the hysteresis restarts) extrapolated from the two nearest saturated runs, clamped to the bracket:')
for Re in (50.0, 100.0, 200.0, 500.0, 800.0):
    kns = sorted(k for r, k in S if r == Re and 0 < k <= 1.0)          # thresholds lie below Kn = 1
    sat = [k for k in kns if S[(Re, k)]['unsteady'] and S[(Re, k)]['amp'] > 1e-3
           and np.isfinite(S[(Re, k)]['sat']) and abs(S[(Re, k)]['sat'] - 1) < 0.05]
    growing = [k for k in kns if k in sat or (np.isfinite(S[(Re, k)]['growth']) and S[(Re, k)]['growth'] > 0
                                              and S[(Re, k)]['amp'] > 1e-5)]    # ignore noise-level lift
    decay = [k for k in kns if k not in growing]
    if not growing or not decay: continue
    kd = min(k for k in decay if k > max(growing)) if any(k > max(growing) for k in decay) else None
    if kd is None: continue
    kp = max(k for k in growing if k < kd)
    gp, gd = S[(Re, kp)]['growth'], S[(Re, kd)]['growth']
    if np.isfinite(gp) and gp > 0 and np.isfinite(gd) and gd < 0:
        kc = zero(kp, gp, kd, gd); how = 'growth %+.4f at Kn %g, %+.4f at Kn %g' % (gp, kp, gd, kd)
    else:
        s2 = [k for k in sat if k <= kp][-2:]
        if len(s2) < 2: continue
        a1, a2 = S[(Re, s2[0])]['amp']**2, S[(Re, s2[1])]['amp']**2
        kc = s2[1] + a2*(s2[1] - s2[0])/(a1 - a2) if a1 > a2 else np.nan
        kc = min(max(kc, kp), kd); how = 'amplitude^2 from Kn %g and %g' % tuple(s2)
    shed.append((kc, Re)); shed_err.append((kp, kd, Re))
    P('   Re = %5g   Knc = %.3f   bracket [%g, %g]   (%s)' % (Re, kc, kp, kd, how))
g0 = sorted((Re, v['growth']) for (Re, k), v in S.items() if k == 0 and 40 < Re < 60 and np.isfinite(v['growth']))
for (r1, g1), (r2, g2) in zip(g0[:-1], g0[1:]):
    if g1 <= 0 < g2:
        re2 = zero(r1, g1, r2, g2); shed.append((0.0, re2)); P('   Kn = 0      Re2 = %.2f   (Legendre 47.5; lit. 46-47.5)' % re2)
KNC = {Re: kc for kc, Re in shed if kc > 0}
if H:
    P('\nhysteresis runs (restart from the saturated shedding state at Kn = 0.3):')
    for (Re, kn), v in sorted(H.items()):
        P('   Re = %5g  Kn = %5g   t_end %5.0f   CL_max %.4f   lift amplitude ratio (2nd/1st half of last window) %.3f   growth %+.4f'
          % (Re, kn, v['t'], v['cl'], v['sat'] if np.isfinite(v['sat']) else 0.0, v['growth']))

# ---------------------------------------------------------------- figure 1
Ls = np.loadtxt('legendre_fig1_separation.dat'); Lh = np.loadtxt('legendre_fig1_shedding.dat')
fig, ax = plt.subplots(figsize=(7.5, 5.5))
ax.semilogy(Ls[:, 0], Ls[:, 1], 'k-o', ms=3.5, lw=1, label='Legendre et al. 2009, separation')
ax.semilogy(Lh[:, 0], Lh[:, 1], 'k--o', ms=3.5, lw=1, label='Legendre et al. 2009, vortex shedding')
if sep:  s = np.array(sorted(sep, key=lambda p: p[1]));  ax.semilogy(s[:, 0], s[:, 1], 's', color='C0', ms=8, mfc='none', mew=1.8, label='Nek5000, separation')
if shed: s = np.array(sorted(shed, key=lambda p: p[1])); ax.semilogy(s[:, 0], s[:, 1], '^', color='C3', ms=8, mfc='none', mew=1.8, label='Nek5000, shedding')
for k1, k2, Re in shed_err: ax.plot([k1, k2], [Re, Re], '-', color='C3', lw=1.2)
ax.text(0.05, 700, 'vortex\nshedding', fontsize=9); ax.text(1.0, 250, 'steady separated wake', fontsize=9); ax.text(0.75, 7.0, 'unseparated wake', fontsize=9)
ax.set(xlabel='Kn = $\\lambda/a$', ylabel='Re', xlim=(0, 2.6), ylim=(5, 2000), title='Stability diagram (Legendre et al. 2009, fig. 1)')
ax.grid(alpha=.3, which='both'); ax.legend(fontsize=8, loc='lower right'); fig.tight_layout(); fig.savefig('legendre_fig1.png', dpi=150)

# ---------------------------------------------------------------- figure 2
L2a = np.loadtxt('legendre_fig2a_all.dat'); L2b = np.loadtxt('legendre_fig2b_all.dat')
LKNC = {50: 0.07, 100: 0.38, 200: 0.47, 500: 0.49, 800: 0.53}           # values used in their fig 2b
fig, ax = plt.subplots(2, 6, figsize=(24, 7.5))
P('\nnormalisation values (Kn = 0 and Kn = inf):')
for j, Re in enumerate((20.0, 50.0, 100.0, 200.0, 500.0, 800.0)):
    a = ax[0, j]; ks = sorted(k for r, k in S if r == Re and k > 0)
    if (Re, 0.0) in S and (Re, -1.0) in S and ks:
        c0, ci = S[(Re, 0.0)]['cd'], S[(Re, -1.0)]['cd']
        P('   Re = %5g   Cd(0) = %.4f   Cd(inf) = %.4f   CL(0) = %.4f   St(0) = %.4f' % (Re, c0, ci, S[(Re, 0.0)]['cl'], S[(Re, 0.0)]['st']))
        L = L2a[L2a[:, 0] == Re]
        a.semilogx(L[:, 1], L[:, 2], 'o', mfc='none', ms=7, color='k', label='Legendre et al.')
        a.semilogx(ks, [(S[(Re, k)]['cd'] - ci)/(c0 - ci) for k in ks], 's-', color='C3', ms=4, label='Nek5000')
    a.set(title='Re = %g' % Re, xlabel='Kn', xlim=(0.007, 30), ylim=(0, 1.02)); a.grid(alpha=.3)
    if j == 0: a.set_ylabel('$C_D^* = (C_D - C_D(\\infty))/(C_D(0) - C_D(\\infty))$'); a.legend(fontsize=8)
    b = ax[1, j]
    if Re == 20.0:
        b.axis('off'); continue
    L = L2b[L2b[:, 0] == Re]
    b.plot(L[:, 1], L[:, 2], 'o', mfc='none', ms=7, color='k', label='Legendre et al. (their Kn$_c$ = %g)' % LKNC[int(Re)])
    if (Re, 0.0) in S and S[(Re, 0.0)]['cl'] > 0:
        ku = [k for k in ks if S[(Re, k)]['unsteady']] + [k for k in ks if not S[(Re, k)]['unsteady']][:2]
        ku = sorted(ku); cl0 = S[(Re, 0.0)]['cl']
        kc = KNC.get(Re, np.nan)
        if np.isfinite(kc):
            b.plot([k/kc for k in ku], [S[(Re, k)]['cl']/cl0 for k in ku], 's-', color='C3', ms=4, label='Nek5000 (own Kn$_c$ = %.3f)' % kc)
        b.plot([k/LKNC[int(Re)] for k in ku], [S[(Re, k)]['cl']/cl0 for k in ku], 's', mfc='none', color='C3', ms=5, alpha=.6, label='Nek5000 with their Kn$_c$')
        hk = sorted(k for r, k in H if r == Re)
        if hk and np.isfinite(KNC.get(Re, np.nan)):
            b.plot([k/KNC[Re] for k in hk], [H[(Re, k)]['cl']/cl0 for k in hk], 'D', color='C2', ms=6,
                   label='Nek5000, restarted from shedding at Kn = 0.3')
    b.set(xlabel='Kn / Kn$_c$', xlim=(-0.02, 1.6), ylim=(-0.02, 1.05)); b.grid(alpha=.3); b.legend(fontsize=7)
    if j == 1: b.set_ylabel('$C_L^* = C_L/C_L(0)$')
fig.suptitle('Legendre et al. (2009) fig. 2: normalised drag (top) and lift amplitude (bottom)', fontsize=13)
fig.tight_layout(); fig.savefig('legendre_fig2.png', dpi=120)

# ---------------------------------------------------------------- figure 1 insets
sel = [('L200_kn0', 'Re = 200, Kn = 0'), ('L200_kn0.1', 'Re = 200, Kn = 0.1'), ('L200_kn0.2', 'Re = 200, Kn = 0.2'),
       ('L200_kn0.4', 'Re = 200, Kn = 0.4'), ('L200_kn0.5', 'Re = 200, Kn = 0.5'),
       ('L800_kn0.2', 'Re = 800, Kn = 0.2'), ('L500_kn0.2', 'Re = 500, Kn = 0.2'), ('L100_kn0.2', 'Re = 100, Kn = 0.2'), ('L50_kn0.2', 'Re = 50, Kn = 0.2')]
fig, ax = plt.subplots(3, 3, figsize=(15, 6.5)); ax = ax.ravel()
for i, (d, title) in enumerate(sel):
    fs = glob.glob(os.path.join(runs, d, 'vort_*.dat'))
    if not fs: ax[i].axis('off'); continue
    v = np.vstack([np.loadtxt(f, ndmin=2) for f in fs])
    m = (v[:, 0] < 10.5) & (np.abs(v[:, 1]) < 2.5)
    ax[i].tricontourf(v[m, 0], v[m, 1], np.clip(v[m, 2], -4.99, 4.99), levels=np.linspace(-5, 5, 41), cmap='RdBu_r')
    ax[i].add_patch(plt.Circle((0, 0), 0.5, color='0.4'))
    ax[i].set(aspect='equal', xlim=(-1, 10.5), ylim=(-2.5, 2.5), title=title); ax[i].set_xticks([]); ax[i].set_yticks([])
fig.suptitle('Vorticity (red positive, colour range $\\pm 5 U/D$) at the end of the runs, as in the insets of Legendre et al. fig. 1', fontsize=11)
fig.tight_layout(); fig.savefig('legendre_fig1_insets.png', dpi=120)
open('legendre_campaign.txt', 'w').write('\n'.join(OUT) + '\n')
