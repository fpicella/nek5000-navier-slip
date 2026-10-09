#!/usr/bin/env python3
"""rerun_analysis.py -- analysis of the re-run of 2026-10-08 (new file; reads, never writes outside this folder).

Inputs  (this folder):  runs/couette_kn*/{run.log|run.log.gz, cyl_drag.dat},  runs/re20_kn*/{cyl_drag.dat, wall_*.dat, run.log*}
        (read only)  ../data/legendre_fig2a_all.dat   digitised Legendre, Lauga & Magnaudet (2009) fig. 2a (Re = 20 rows)
                     ../data/re20_sweep.dat           the earlier (2026-09-11) Re = 20 sweep
                     ../data/legendre_campaign.txt    the earlier runs at Re = 20 Kn = 0.322, 0.391, 0.529 and wall analysis
Outputs (this folder):  rerun_results.json, rerun_results.txt, legendre_campaign_rerun.txt (campaign file with the Re = 20 rows replaced)

Definitions (identical to the original Re = 20 comparison):
  Cd*(Kn) = (Cd(Kn) - Cd(inf)) / (Cd(0) - Cd(inf)),  Cd(0) = no-slip run, Cd(inf) = shear-free run (Kn = -1), both re-run;
  each Legendre marker is compared at the NOMINAL Kn nearest (in log) to its digitised Kn;
  max_abs, mean_abs, rms of (own - Legendre) over the markers, max_abs_in_cd = max_abs * (Cd(0) - Cd(inf)).
Extra (not in the earlier analysis): the same comparison with the own curve read at the DIGITISED Kn of each marker
  (monotone cubic interpolation in log Kn through all own Re = 20 runs), because the digitised Kn differ from the
  nominal ones by up to 1.3 %.
"""
import os, re, gzip, glob, json, sys
import numpy as np
from scipy.interpolate import PchipInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OLDRUNS = os.path.join(HERE, '..', 'runs')          # the original campaign runs (not in the repository)
RUNS = os.path.join(HERE, 'runs')
LEG = os.path.join(DATA, 'legendre_fig2a_all.dat')
OLDSW = os.path.join(DATA, 're20_sweep.dat')
OLDCAMP = os.path.join(DATA, 'legendre_campaign.txt')
NOMINAL = np.array([0.01, 0.02, 0.05, 0.07, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 1, 2, 2.4, 5, 10, 20])
OUT = []


def P(s=''):
    print(s)
    OUT.append(s)


def opentxt(path):
    if os.path.exists(path):
        return open(path, errors='replace')
    if os.path.exists(path + '.gz'):
        return gzip.open(path + '.gz', 'rt', errors='replace')
    return None


def runname(prefix, kn):
    return '%s_kn%g' % (prefix, kn)


# ------------------------------------------------------------------ (i) circular Couette
def couette_exact(kn, a=0.5, R=5.0, om=1.0, mu=1.0):
    lam = kn * a
    if kn > 0:
        B = -om * a / (1.0 / a - a / R**2 + 2.0 * lam / a**2)
        u = -2.0 * lam * B / a**2
    elif kn < 0:
        B, u = 0.0, om * a
    else:
        B, u = -om * a / (1.0 / a - a / R**2), 0.0
    return u, -4.0 * np.pi * mu * B


res = {}
P('(i) Circular-Couette annulus a = 0.5, R = 5, outer wall omega = 1, mu = 1 (steady state); slip on the inner cylinder')
P('    exact: u_th = A r + B/r, B = -omega a / (1/a - a/R^2 + 2 lambda/a^2), u_th(a) = -2 lambda B / a^2, torque = -4 pi mu B')
P('')
P('  Kn        wall velocity: Nek (cyl_drag.dat, 9 digits)   exact       rel. diff      torque: Nek (log, 7 digits)  exact       diff        Navier residual (final)  t_end  steps')
res['couette'] = {}
for kn, lab in ((1, '1'), (0, '0 (no slip)'), (-1, '-1 (shear-free)')):
    d = os.path.join(RUNS, runname('couette', kn))
    if not os.path.exists(os.path.join(d, 'cyl_drag.dat')):
        P('  %-6s  missing' % lab)
        continue
    a = np.loadtxt(os.path.join(d, 'cyl_drag.dat'))
    last = a[-1]
    f = opentxt(os.path.join(d, 'run.log'))
    txt = f.read().split('\n')
    f.close()
    if not any('run successful' in l for l in txt):
        P('  %-6s  not finished' % lab)
        continue
    cl = [l.split() for l in txt if l.rstrip().endswith('couette')]
    c = [float(x) for x in cl[-1][:8]]          # step time usmax utha torque_num torque_exact rmax rout
    ue, te = couette_exact(kn)
    un = last[6]                                # us_max in cyl_drag.dat
    tn = c[4]
    rel_u = (un - ue) / ue if ue != 0 else un - ue
    P('  %-15s %.9f                  %.9f   %+.2e        %.6f                  %.6f   %+.1e      %.2e         %5.1f  %d'
      % (lab, un, ue, rel_u, tn, te, tn - te, last[8], last[1], last[0]))
    res['couette'][str(kn)] = dict(wall_velocity_nek=un, wall_velocity_exact=ue, rel_diff=rel_u, torque_nek=tn,
                                   torque_exact=te, torque_log_exact=c[5], res_navier_last=last[8], t_end=last[1],
                                   steps=int(last[0]), rout=c[7])
P('')

# ------------------------------------------------------------------ (ii) Re = 20 cylinder
def load_run(kn):
    d = os.path.join(RUNS, runname('re20', kn))
    f = os.path.join(d, 'cyl_drag.dat')
    if not os.path.exists(f):
        return None
    a = np.loadtxt(f)
    if a.ndim == 1 or len(a) < 2:
        return None
    t = a[-1, 1]
    m = a[:, 1] >= t - 10.0
    ok = False
    lg = opentxt(os.path.join(d, 'run.log'))
    if lg is not None:
        ok = any('run successful' in l for l in lg)
        lg.close()
    return dict(kn=kn, step=int(a[-1, 0]), t=t, cd=a[-1, 2], cdp=a[-1, 3], cdv=a[-1, 4], usmin=a[-1, 7], res=a[-1, 8],
                dudt=a[-1, 10], drift=a[m, 2].max() - a[m, 2].min(), ok=ok, dir=d)


def wallsep(d):
    """slope of the wall vorticity at the rear stagnation point (fit over |th| < 5 deg) and separation angle (zero crossing
    on the upper surface) -- same method as scripts/plot_legendre.py"""
    fs = [f for f in glob.glob(os.path.join(d, 'wall_*.dat')) if os.path.getsize(f) > 0]
    if not fs:
        return np.nan, np.nan
    W = np.vstack([np.loadtxt(f, ndmin=2) for f in fs])
    W = W[np.argsort(W[:, 0])]
    th, om = W[:, 0], W[:, 2]
    k = np.abs(th) < np.radians(5.0)
    A = np.c_[th[k], th[k]**3]
    c = np.linalg.lstsq(A, om[k], rcond=None)[0][0]
    u = (th > 0) & (th < np.pi / 2)
    tu, ou = th[u], om[u]
    ths = 0.0
    if ou[0] > 0:
        j = np.where((ou[:-1] > 0) & (ou[1:] <= 0))[0]
        if len(j):
            j = j[0]
            ths = np.degrees(tu[j] + (tu[j + 1] - tu[j]) * ou[j] / (ou[j] - ou[j + 1]))
    return c, ths


allkn = sorted({float(re.sub(r'^re20_kn', '', os.path.basename(p))) for p in glob.glob(os.path.join(RUNS, 're20_kn*'))})
R = {kn: load_run(kn) for kn in allkn}
R = {k: v for k, v in R.items() if v is not None}
for k in sorted(R):
    if not R[k]['ok']:
        print('  (run Kn = %g is not finished: left out of every table)' % k)
R = {k: v for k, v in R.items() if v['ok']}      # only runs that ended with 'run successful'
old = np.loadtxt(OLDSW)                      # Kn t Cd Cd_p Cd_v us_max us_min res dudt drift (6 significant digits only)
oldcd = {round(r[0], 6): r for r in old}
# earlier runs at full precision: runs_data/re20_kn*/cyl_drag.dat (sweep) and runs_data/L20_kn*/cyl_drag.dat (campaign)
OLDHIST = {}
for pth in glob.glob(os.path.join(OLDRUNS, 're20_kn*')) + glob.glob(os.path.join(OLDRUNS, 'L20_kn*')):
    kk = float(re.sub(r'^.*_kn', '', os.path.basename(pth)))
    f_ = os.path.join(pth, 'cyl_drag.dat')
    if os.path.exists(f_) and 'L20' not in pth or (os.path.exists(f_) and round(kk, 6) not in OLDHIST):
        OLDHIST[round(kk, 6)] = np.loadtxt(f_)
# the three extra Kn of the separation bracket: earlier values from the campaign file
oldx = {}
for line in open(OLDCAMP):
    s = line.split('|')[0].split()
    if len(s) >= 10 and s[0] == '20' and s[3] in ('True', 'False'):
        oldx[round(float(s[1]), 6)] = float(s[4])

P('(ii) 2D cylinder, Re = 20, steady, 1040 elements, N = 7, dt = 5e-3, BDF2 + OIFS; stop at max|du/dt| < 1e-6 or t = 150')
P('')
P('  Kn         steps   t_end     Cd         Cd_p       Cd_v      us_min     res_Navier(final)  dCd/dt(last)  Cd drift(10 t.u.)  run ok   Cd earlier     Cd new - old')
res['runs'] = {}
for kn in sorted(R):
    r = R[kn]
    h = OLDHIST.get(round(kn, 6))
    cdold = h[-1, 2] if h is not None else oldx.get(round(kn, 6), np.nan)
    hist = ''
    if h is not None:
        n_ = min(len(h), len(np.loadtxt(os.path.join(r['dir'], 'cyl_drag.dat'))))
        a_ = np.loadtxt(os.path.join(r['dir'], 'cyl_drag.dat'))
        common = np.intersect1d(h[:, 0], a_[:, 0])
        dh = np.abs(h[np.isin(h[:, 0], common), 2] - a_[np.isin(a_[:, 0], common), 2]).max() if len(common) else np.nan
        hist = 'steps %d vs %d, max|dCd| over %d common diagnostic steps = %.1e' % (len(a_), len(h), len(common), dh)
        r['hist'] = (len(a_), len(h), len(common), float(dh), int(a_[-1, 0]), int(h[-1, 0]))
    P('  %-9g %6d  %6.2f  %.6f  %.6f  %.6f  %+.2e  %.2e       %.2e      %.2e          %s     %.8f     %+.2e'
      % (kn, r['step'], r['t'], r['cd'], r['cdp'], r['cdv'], r['usmin'], r['res'], r['dudt'], r['drift'],
         'yes' if r['ok'] else 'NO', cdold, r['cd'] - cdold if np.isfinite(cdold) else np.nan))
    if hist:
        P('             history vs earlier run: ' + hist)
    res['runs'][str(kn)] = dict(steps=r['step'], t_end=r['t'], cd=r['cd'], cd_p=r['cdp'], cd_v=r['cdv'], us_min=r['usmin'],
                                res_navier_final=r['res'], dudt_last=r['dudt'], cd_drift_last10=r['drift'], run_ok=r['ok'],
                                cd_earlier=None if not np.isfinite(cdold) else float(cdold), history=r.get('hist'))
P('')

have = 0.0 in R and -1.0 in R
if have:
    cd0, cdi = R[0.0]['cd'], R[-1.0]['cd']
    P('  Cd(Kn = 0)   = %.5f   (earlier %.5f; Legendre table 2: 2.035; Dennis & Chang 1970: 2.045)' % (cd0, OLDHIST[0.0][-1, 2]))
    P('  Cd(Kn = inf) = %.5f   (earlier %.5f; Legendre table 2: 1.33)' % (cdi, OLDHIST[-1.0][-1, 2]))
    P('')
    kn_pos = np.array(sorted(k for k in R if k > 0))
    cds = {k: (R[k]['cd'] - cdi) / (cd0 - cdi) for k in kn_pos}
    old0, oldi = OLDHIST[0.0][-1, 2], OLDHIST[-1.0][-1, 2]
    L = np.loadtxt(LEG)
    L = L[L[:, 0] == 20]
    L = L[np.argsort(L[:, 1])]
    # own curve, monotone cubic in log Kn, through all own runs (for reading at the digitised Kn)
    f = PchipInterpolator(np.log(kn_pos), np.array([cds[k] for k in kn_pos]))
    P('  Cd* = (Cd - Cd(inf)) / (Cd(0) - Cd(inf)) against the digitised Legendre et al. (2009) fig. 2a, Re = 20 (11 markers)')
    P('  nominal Kn   digitised Kn   Legendre Cd*   Cd* re-run   diff (nominal Kn)   Cd* earlier   re-run - earlier   own Cd* at digitised Kn   diff (digitised Kn)')
    rows = []
    for kd, cl in zip(L[:, 1], L[:, 2]):
        k = NOMINAL[np.argmin(np.abs(np.log(NOMINAL / kd)))]
        if k not in cds:
            continue
        cs = cds[k]
        csold = (OLDHIST[round(k, 6)][-1, 2] - oldi) / (old0 - oldi)
        csd = float(f(np.log(kd)))
        rows.append((k, kd, cl, cs, cs - cl, csold, cs - csold, csd, csd - cl))
        P('  %-10g   %-12.5f   %.4f         %.4f       %+.4f             %.4f        %+.5f            %.4f                    %+.4f'
          % rows[-1])
    A = np.array(rows)
    d = A[:, 4]
    i = int(np.argmax(np.abs(d)))
    dd = A[:, 8]
    j = int(np.argmax(np.abs(dd)))
    cdrng = cd0 - cdi
    stat = dict(n=len(d), max_abs=float(np.abs(d).max()), at_kn=float(A[i, 0]), mean_abs=float(np.abs(d).mean()),
                rms=float(np.sqrt((d**2).mean())), cd0=float(cd0), cdinf=float(cdi), max_abs_in_cd=float(np.abs(d).max() * cdrng),
                max_abs_digitised_kn=float(np.abs(dd).max()), at_kn_digitised=float(A[j, 0]),
                mean_abs_digitised_kn=float(np.abs(dd).mean()), rms_digitised_kn=float(np.sqrt((dd**2).mean())),
                max_abs_vs_earlier=float(np.abs(A[:, 6]).max()),
                max_abs_cd_vs_earlier=float(max(abs(R[k]['cd'] - OLDHIST[round(k, 6)][-1, 2]) for k in R if round(k, 6) in OLDHIST)))
    P('')
    P('  RESULT (nominal Kn, convention of the earlier analysis): %d markers, max |dCd*| = %.4f at Kn = %g, mean |dCd*| = %.4f, rms = %.4f;'
      % (stat['n'], stat['max_abs'], stat['at_kn'], stat['mean_abs'], stat['rms']))
    P('           in Cd units the maximum is %.4f (%.2f %% of Cd(0)).' % (stat['max_abs_in_cd'], 100 * stat['max_abs_in_cd'] / cd0))
    P('  RESULT (own curve read at the digitised Kn): max |dCd*| = %.4f at Kn = %g, mean = %.4f, rms = %.4f.'
      % (stat['max_abs_digitised_kn'], stat['at_kn_digitised'], stat['mean_abs_digitised_kn'], stat['rms_digitised_kn']))
    P('  Earlier analysis (fig_validation_numbers.json): max 0.0035898 at Kn = 1, mean 0.000856, rms 0.001336.')
    P('  Re-run minus earlier: max |dCd*| = %.2e, max |dCd| over all runs of both sets = %.2e.' % (stat['max_abs_vs_earlier'], stat['max_abs_cd_vs_earlier']))
    res['cdstar_vs_legendre'] = stat
    res['cdstar_table'] = [dict(kn_nominal=a[0], kn_digitised=a[1], legendre=a[2], rerun=a[3], diff=a[4], earlier=a[5],
                                rerun_minus_earlier=a[6], own_at_digitised_kn=a[7], diff_digitised_kn=a[8]) for a in A]
    P('')

# ------------------------------------------------------------------ (iii) steady separation
P('(iii) Steady separation at Re = 20 (wall vorticity at the rear stagnation point, fit over |th| < 5 deg; wall_*.dat at the last step)')
P('  Kn         slope of wall vorticity   separation angle (deg)   us_min    (earlier campaign: slope / angle)')
oldwall = {}
for line in open(OLDCAMP):
    if '|' in line:
        s = line.split('|')
        t = s[0].split()
        w = s[1].split()
        if len(t) >= 10 and t[0] == '20' and t[3] in ('True', 'False') and len(w) >= 2:
            try:
                oldwall[round(float(t[1]), 6)] = (float(w[0]), float(w[1]))
            except ValueError:
                pass
sl = {}
res['separation'] = {}
for kn in sorted(R):
    c, th = wallsep(R[kn]['dir'])
    sl[kn] = (c, th)
    o = oldwall.get(round(kn, 6))
    P('  %-9g  %+9.4f                %7.2f                  %+.2e   %s'
      % (kn, c, th, R[kn]['usmin'], ('%+.4f / %.2f' % o) if o else '-'))
    res['separation'][str(kn)] = dict(slope=float(c), angle_deg=float(th), earlier=o)
pos = sorted(k for k in sl if k > 0 and np.isfinite(sl[k][0]))
sep = None
for k1, k2 in zip(pos[:-1], pos[1:]):
    if sl[k1][0] > 0 >= sl[k2][0]:
        sep = k1 + (k2 - k1) * sl[k1][0] / (sl[k1][0] - sl[k2][0])
        res['kn_sep_bracket'] = [k1, k2]
        break
if sep is not None:
    P('  Kn_sep(Re = 20) = %.3f  (linear zero of the slope between Kn = %g and %g); earlier 0.476; Legendre et al. fig. 1 (digitised, interpolated) 0.470'
      % (sep, res['kn_sep_bracket'][0], res['kn_sep_bracket'][1]))
    res['kn_sep'] = float(sep)
    # same bracket as the earlier analysis (0.391 / 0.5), for a like-for-like value
    if 0.391 in sl and 0.5 in sl:
        s2 = 0.391 + (0.5 - 0.391) * sl[0.391][0] / (sl[0.391][0] - sl[0.5][0])
        P('  Kn_sep from the earlier bracket (0.391 / 0.5): %.3f' % s2)
        res['kn_sep_391_05'] = float(s2)
if sep is not None and all(k in sl for k in (0.322, 0.391, 0.5, 0.529)):
    kk_ = np.array([0.322, 0.391, 0.5, 0.529]); ss_ = np.array([sl[k][0] for k in kk_])
    alt = {}
    for deg, sel in ((2, slice(0, 4)), (3, slice(0, 4)), (2, slice(1, 4))):
        rt = np.roots(np.polyfit(kk_[sel], ss_[sel], deg)); rt = [x.real for x in rt if abs(x.imag) < 1e-12 and 0.35 < x.real < 0.55]
        alt['deg%d_%s' % (deg, 'Kn0.322-0.529' if sel == slice(0, 4) else 'Kn0.391-0.529')] = rt[0] if rt else None
    P('  Sensitivity of Kn_sep to the interpolation (polynomial through the slopes): ' + ', '.join('%s: %.4f' % (a, b) for a, b in alt.items()))
    res['kn_sep_polynomial_alternatives'] = alt
if 0.0 in sl:
    P('  No-slip separation angle: %.2f deg (earlier campaign 43.68; Dennis & Chang 1970: 43.7, as quoted in notes/README.md)' % sl[0.0][1])
P('  Separated rear (us_min < 0) for Kn <= %s, attached for Kn >= %s'
  % (max([k for k in R if k > 0 and R[k]['usmin'] < -1e-6] or [np.nan]), min([k for k in R if k > 0 and R[k]['usmin'] >= -1e-6] or [np.nan])))

# ------------------------------------------------------------------ campaign-format file for the figure script (Re = 20 rows replaced)
if have:
    lines = open(OLDCAMP).read().split('\n')
    new = []
    done = set()
    for line in lines:
        t = line.split('|')[0].split()
        if len(t) >= 10 and t[0] == '20' and t[3] in ('True', 'False') and '|' in line:
            kn = float(t[1])
            if round(kn, 6) in [round(k, 6) for k in R]:
                kk = [k for k in R if round(k, 6) == round(kn, 6)][0]
                r = R[kk]
                c, th = sl[kk]
                new.append('%5g %8g %6.0f  False   %.5f   %6.4f   %6s   %+.5f   %+6s  %5s  |  %+9.4f %7.2f'
                           % (20, kk, r['t'], r['cd'], 0.0, 'nan', float(t[7]), 'nan', 'nan', c, th))
                done.add(kk)
                continue
        new.append(line)
    open(os.path.join(HERE, 'legendre_campaign_rerun.txt'), 'w').write('\n'.join(new))
    P('')
    P('  legendre_campaign_rerun.txt: copy of data/legendre_campaign.txt with the Re = 20 rows of the runs %s replaced by the re-run Cd and wall slope'
      % sorted(done))
    P('  (unused columns CL_max/St/growth/cw copied or set to the same placeholders; Re >= 50 rows are the earlier ones, unchanged).')
    res['replaced_rows'] = sorted(done)

json.dump(res, open(os.path.join(HERE, 'rerun_results.json'), 'w'), indent=1, default=float)
open(os.path.join(HERE, 'rerun_results.txt'), 'w').write('\n'.join(OUT) + '\n')
