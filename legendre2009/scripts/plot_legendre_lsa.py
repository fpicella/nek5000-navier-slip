#!/usr/bin/env python3
"""Legendre, Lauga & Magnaudet (2009) figures 1, 3 and 4 against this work:
DNS (legendre_campaign.txt, wall dumps wall_dns/) and nekStab linear stability
(../nekstab/results: lsa_summary.txt, wall_lsa/ = wall distributions of the
Newton base flows, incl. the unstable steady branch of the bf_* continuation).
Digitised: legendre_fig1_{separation,shedding}.dat, legendre_fig3_all.dat,
legendre_fig4_curves.dat (digitize_fig3.py, digitize_fig4.py).
Outputs: legendre_fig1_lsa.png, legendre_fig3.png, legendre_fig4.png,
legendre_fig134_summary.txt"""
import os, re, glob, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

D = os.path.dirname(os.path.abspath(__file__)); NS = os.path.join(os.path.dirname(D), 'nekstab', 'results')
txt = open(os.path.join(D, 'legendre_campaign.txt')).read()
dns = {}
for line in txt.splitlines():
    t = line.split()
    if len(t) >= 8 and t[3] in ('True', 'False'):
        dns[(float(t[0]), float(t[1]))] = dict(uns=t[3] == 'True', st=float(t[6]))
sep = {int(m.group(1)): (float(m.group(2)), float(m.group(3)), float(m.group(4))) for m in
       re.finditer(r'Re =\s+(\d+)\s+Kn_sep = ([\d.]+)\s+\(slope \S+ at Kn ([\d.]+), \S+ at Kn ([\d.]+)', txt)}
m = re.search(r'Re1 = ([\d.]+)', txt); re1 = float(m.group(1)) if m else 6.28
kd = {int(m.group(1)): (float(m.group(2)), float(m.group(3)), float(m.group(4))) for m in
      re.finditer(r'Re =\s+(\d+)\s+Knc = ([\d.]+)\s+bracket \[([\d.]+), ([\d.]+)\]', txt)}
m = re.search(r'Re2 = ([\d.]+)', txt); re2 = float(m.group(1)) if m else np.nan
ls = open(os.path.join(NS, 'lsa_summary.txt')).read()
knc = {float(m.group(1)): (float(m.group(2)), float(m.group(3))) for m in
       re.finditer(r'Re =\s+([\d.]+):\s+Knc = ([\d.]+) \(quadratic [^)]*\)\s+St_c = ([\d.]+)', ls)}
m = re.search(r'Re_c = ([\d.]+) \(quadratic [^)]*\)\s+St_c = ([\d.]+)', ls); rec, stc0 = float(m.group(1)), float(m.group(2))

# surface vorticity omega_max a/U (a = 0.5, U = 1) from the wall dumps
def wmax(files):
    W = np.vstack([np.loadtxt(f, ndmin=2) for f in files if os.path.getsize(f) > 0])
    return 0.5*np.abs(W[:, 2]).max()
wd, wb = {}, {}
for d in glob.glob(os.path.join(D, 'wall_dns', '*')):
    m = re.match(r'(?:L|re)([\d.]+)_kn(-?[\d.]+)$', os.path.basename(d))
    if m and not dns.get((float(m.group(1)), float(m.group(2))), {}).get('uns', False):
        wd[(float(m.group(1)), float(m.group(2)))] = wmax(glob.glob(d + '/diag/wall_*.dat'))
for d in glob.glob(os.path.join(NS, 'wall_lsa', '*Re*_kn*')):
    m = re.match(r'(?:bf_)?Re([\d.]+)_kn(-?[\d.]+)$', os.path.basename(d))
    if m: wb[(float(m.group(1)), float(m.group(2)))] = wmax(glob.glob(d + '/wall/wall_*.dat'))
wall = dict(wd); wall.update(wb)                     # steady states: DNS where steady, Newton base flows
snap = {}                                          # shedding DNS: instantaneous values (wall_snap.txt)
if os.path.exists(os.path.join(D, 'wall_snap.txt')):
    for line in open(os.path.join(D, 'wall_snap.txt')):
        if line.startswith('#') or not line.strip(): continue
        n, _, w = line.split()
        m = re.match(r'(?:L|re)([\d.]+)_kn(-?[\d.]+)$', n)
        if m: snap.setdefault((float(m.group(1)), float(m.group(2))), []).append(float(w))

def interp_kn(Re, kn0, src):
    p = sorted((k, w) for (r, k), w in src.items() if r == Re and k >= 0)
    for (k1, w1), (k2, w2) in zip(p[:-1], p[1:]):
        if k1 <= kn0 <= k2: return w1 + (w2 - w1)*(kn0 - k1)/(k2 - k1)
    return np.nan
def interp_re(re0, kn, src):
    p = sorted((r, w) for (r, k), w in src.items() if k == kn)
    for (r1, w1), (r2, w2) in zip(p[:-1], p[1:]):
        if r1 <= re0 <= r2: return np.exp(np.interp(np.log(re0), [np.log(r1), np.log(r2)], [np.log(w1), np.log(w2)]))
    return np.nan

S = ['Legendre, Lauga & Magnaudet (2009) figures 1, 3, 4 vs this work (DNS + nekStab linear stability)\n']
# ---------------- figure 1: stability diagram
lsep = np.loadtxt(os.path.join(D, 'legendre_fig1_separation.dat')); lshd = np.loadtxt(os.path.join(D, 'legendre_fig1_shedding.dat'))
fig, ax = plt.subplots(figsize=(6.6, 5.0))
ax.plot(lsep[:, 0], lsep[:, 1], '--', color='k', lw=1, label='Legendre et al.: separation')
ax.plot(lshd[:, 0], lshd[:, 1], '-', color='k', lw=1, label='Legendre et al.: vortex shedding')
rs = sorted(sep); ax.plot([0] + [sep[r][0] for r in rs], [re1] + rs, 'o', color='C1', mfc='none', ms=6, label='DNS: separation')
if kd:
    r = np.array(sorted(kd)); v = np.array([kd[x] for x in r])
    ax.errorbar(v[:, 0], r, xerr=[v[:, 0] - v[:, 1], v[:, 2] - v[:, 0]], fmt='s', color='C3', ms=5, capsize=3, label='DNS: shedding (bracket)')
kr = sorted(knc); ax.plot([0] + [knc[r][0] for r in kr], [rec] + kr, 'D-', color='C0', ms=5, lw=1.5, label='linear stability (nekStab): shedding')
ax.set_yscale('log'); ax.set_xlim(-0.05, 2.7); ax.set_ylim(4, 1000); ax.set_xlabel('Kn = $\\lambda/a$'); ax.set_ylabel('Re')
ax.text(1.75, 28, 'unseparated', fontsize=8); ax.text(0.9, 30, 'steady\nseparated', fontsize=8); ax.text(0.03, 400, 'vortex\nshedding', fontsize=8)
ax.grid(alpha=.3, which='both'); ax.legend(fontsize=7.5, loc='lower right'); ax.set_title('Figure 1: stability diagram', fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(D, 'legendre_fig1_lsa.png'), dpi=150)

# ---------------- figure 3: St* = St(Kn)/St(0)
l3 = np.loadtxt(os.path.join(D, 'legendre_fig3_all.dat'))
lk = {50: 0.09, 100: 0.38, 200: 0.47, 500: 0.49, 800: 0.53}           # their Knc (fig. 2 caption)
mk = {50: 's', 100: 'D', 200: '^', 500: 'v', 800: '*'}
fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
S.append('figure 3: St* at the last shedding point (DNS) and at onset (linear stability), vs Legendre\n'
         '   Re   St(0)_DNS  Kn_last  St*_DNS  Kn/Knc   St*_c(LSA)  | Legendre St* at Kn/Knc~1')
for i, Re in enumerate((50, 100, 200, 500, 800)):
    c = 'C%d' % i; st0 = dns[(Re, 0.0)]['st']; kc = knc.get(Re, (np.nan, np.nan))[0]
    L = l3[l3[:, 0] == Re]
    for ax, xs in ((axs[0], L[:, 1]), (axs[1], L[:, 1]*lk[Re])):
        ax.plot(xs, L[:, 2], mk[Re] + '-', color=c, mfc='none', ms=6, lw=0.8, label='Legendre, Re = %d' % Re if ax is axs[0] else None)
    p = sorted((k, v['st']/st0) for (r, k), v in dns.items() if r == Re and v['uns'])
    k, s = np.array(p).T
    for ax, xs in ((axs[0], k/kc), (axs[1], k)):
        ax.plot(xs, s, mk[Re], color=c, ms=7, label='this work (DNS), Re = %d' % Re if ax is axs[0] else None)
    if Re in knc:
        for ax, x in ((axs[0], 1.0), (axs[1], kc)):
            ax.plot([x], [knc[Re][1]/st0], 'o', color=c, mec='k', ms=9, zorder=5)
    S.append('  %4d  %8.4f  %7.3f  %7.3f  %6.3f   %9.3f   | %6.3f'
             % (Re, st0, k[-1], s[-1], k[-1]/kc, knc[Re][1]/st0 if Re in knc else np.nan, L[-1, 2]))
axs[0].plot([], [], 'o', color='w', mec='k', ms=9, label='linear stability: onset frequency')
axs[0].set_xlabel('Kn / Kn$_c$ (each work with its own Kn$_c$)'); axs[1].set_xlabel('Kn')
for ax in axs: ax.set_ylabel('St* = St(Kn)/St(0)'); ax.grid(alpha=.3)
axs[0].legend(fontsize=6.5, ncol=2); axs[0].set_title('Figure 3, as in the paper', fontsize=10)
axs[1].set_title('same data against Kn', fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(D, 'legendre_fig3.png'), dpi=150)

# ---------------- figure 4: surface vorticity
l4 = {}
for line in open(os.path.join(D, 'legendre_fig4_curves.dat')):
    if line.startswith('#'): continue
    c, r, w = line.split(); l4.setdefault(c, []).append((float(r), float(w)))
kns = [(0.0, 'Kn0'), (0.05, 'Kn0.05'), (0.1, 'Kn0.1'), (0.2, 'Kn0.2'), (0.4, 'Kn0.4'), (1.0, 'Kn1'),
       (2.2, 'Kn2.2'), (5.0, 'Kn5'), (-1.0, 'Kninf')]
cmap = plt.get_cmap('viridis')
fig, ax = plt.subplots(figsize=(7.0, 5.6))
S.append('\nfigure 4: omega_max a/U, this work / Legendre on their iso-Kn curves (steady states only)')
for j, (kn, key) in enumerate(kns):
    col = cmap(j/(len(kns) - 1)); c = np.array(l4[key])
    lw, stl = (1.8, '-') if kn == 0 else ((1.4, '-.') if kn < 0 else (0.9, '-'))
    ax.plot(c[:, 0], c[:, 1], stl, color=col, lw=lw, label='Kn = %s' % ('$\\infty$' if kn < 0 else '%g' % kn))
    for src, mkr in ((wd, 's'), (wb, 'D')):
        pts = sorted((r, w) for (r, k), w in src.items() if k == kn)
        if not pts: continue
        r, w = np.array(pts).T
        ax.plot(r, w, mkr, color=col, mec='k', mew=.4, ms=5.5)
        lw_ = np.exp(np.interp(np.log(r), np.log(c[:, 0]), np.log(c[:, 1]), left=np.nan, right=np.nan))
        ok = np.isfinite(lw_)
        if ok.any(): S.append('  Kn = %-5s %-4s Re %s: ratio %s' % ('inf' if kn < 0 else '%g' % kn, 'DNS' if src is wd else 'BF',
                              ' '.join('%g' % x for x in r[ok]), ' '.join('%.3f' % x for x in (w/lw_)[ok])))
    sp = sorted((r, np.mean(v), min(v), max(v)) for (r, k), v in snap.items() if k == kn)
    if sp:
        r, mu, lo, hi = np.array(sp).T
        ax.errorbar(r, mu, yerr=[mu - lo, hi - mu], fmt='o', color=col, mfc='none', mew=1.6, ms=8, capsize=2, zorder=6)
        lw_ = np.exp(np.interp(np.log(r), np.log(c[:, 0]), np.log(c[:, 1]), left=np.nan, right=np.nan))
        S.append('  Kn = %-5s DNS shedding, instantaneous (snapshots) Re %s: ratio %s' % ('%g' % kn,
                 ' '.join('%g' % x for x in r), ' '.join('%.3f' % x for x in mu/lw_)))
cs = np.array(l4['sep']); ch = np.array(l4['shed'])
ax.plot(cs[:, 0], cs[:, 1], '--', color='k', lw=1.3, label='Legendre: separation')
ax.plot(ch[:, 0], ch[:, 1], ':', color='k', lw=1.8, label='Legendre: shedding')
ws = [(re1, interp_re(re1, 0.0, wd))] + [(r, interp_kn(r, sep[r][0], wd)) for r in rs]
wc = [(rec, interp_re(rec, 0.0, wb))] + [(r, interp_kn(r, knc[r][0], wb)) for r in kr]
ws = np.array([p for p in ws if np.isfinite(p[1])]); wc = np.array([p for p in wc if np.isfinite(p[1])])
ax.plot(ws[:, 0], ws[:, 1], 'o--', color='C1', mfc='none', ms=6, lw=1, label='this work: separation (DNS)')
ax.plot(wc[:, 0], wc[:, 1], 'o:', color='C3', mec='k', ms=6, lw=1.5, label='this work: shedding (linear stability)')
ax.plot([], [], 's', color='0.6', mec='k', label='this work: steady DNS'); ax.plot([], [], 'D', color='0.6', mec='k', label='this work: Newton base flow')
ax.plot([], [], 'o', color='0.3', mfc='none', mew=1.6, ms=8, label='this work: shedding DNS, instantaneous')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(2, 1000); ax.set_ylim(1, 40)
ax.set_xlabel('Re'); ax.set_ylabel('$\\omega_{max}/(U/a)$'); ax.grid(alpha=.3, which='both')
ax.legend(fontsize=6.5, ncol=2, loc='upper left'); ax.set_title('Figure 4: maximum surface vorticity', fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(D, 'legendre_fig4.png'), dpi=150)
S.append('\ncritical surface vorticity omega_max a/U (this work):')
S.append('  separation (DNS, at Kn_sep): ' + '  '.join('Re %g: %.2f' % tuple(p) for p in ws))
S.append('  shedding (linear stability, at Knc): ' + '  '.join('Re %g: %.2f' % tuple(p) for p in wc))
S.append('  Legendre: separation %.2f..%.2f, shedding %.2f..%.2f (their "5 < omega_max a/U < 7")'
         % (cs[:, 1].min(), cs[:, 1].max(), ch[:, 1].min(), ch[:, 1].max()))
open(os.path.join(D, 'legendre_fig134_summary.txt'), 'w').write('\n'.join(S) + '\n')
print('\n'.join(S))
