#!/usr/bin/env python3
"""Summarise the timing test (scripts/timing_test.sh): wall time, time per
step and Nek's runstat breakdown for runs/timing/<round>_<case>/run.log.
Usage: timing_parse.py [runs/timing]"""
import sys, os, re, glob

T = sys.argv[1] if len(sys.argv) > 1 else 'runs/timing'
CASE = {'A': 'no-slip, laplacian form (standard)',
        'B': 'no-slip, stress form',
        'C': 'Kn = 1e-3, stress form + Robin',
        'D': 'Kn = 1e-3, diagnostics every step',
        'E': 'Kn = 1e-4, stress form + Robin',
        'F': 'Kn = 1e-2, stress form + Robin',
        'G': 'Kn = 1, stress form + Robin',
        'H': 'shear-free (shl, no Robin term)'}

def parse(log):
    s = open(log, errors='replace').read()
    r = {}
    m = re.findall(r'^Step\s+(\d+), t=\s*([^,\s]+)', s, re.M)
    r['steps'], r['tend'] = int(m[-1][0]), float(m[-1][1])
    r['total'] = float(re.search(r'total elapsed time\s*:\s*(\S+)', s).group(1))
    r['perstep'] = float(re.search(r'time/timestep\s*:\s*(\S+)', s).group(1))
    rs = s[s.rfind('runtime statistics:'):]
    for line in rs.splitlines()[1:]:
        w = line.split()
        if len(w) >= 4 and w[1] == 'time':
            try: nums = [float(x) for x in w[2:]]
            except ValueError: continue
            r[w[0]] = nums                  # [count,] total, fraction
    vel = [int(x) for x in re.findall(r'^\s*\d+\s+(?:Helmh3 fluid|Hmholtz VEL[XY]:?|Helmholtz VEL[XY]:?)\s+(\d+)', s, re.M)]
    pre = [int(x) for x in re.findall(r'U-PRES gmres\s+(\d+)', s)]
    r['vit'] = sum(vel)/max(r['steps'], 1)
    r['pit'] = sum(pre)/max(r['steps'], 1)
    d = os.path.join(os.path.dirname(log), 'cyl_drag.dat')
    last = [l for l in open(d) if not l.startswith('#')][-1].split()
    r['cd'] = float(last[2])
    return r

def tot(r, k): return r[k][-2] if k in r else float('nan')

rows = []
for log in sorted(glob.glob(os.path.join(T, 'r*_?', 'run.log'))):
    rc = os.path.basename(os.path.dirname(log)); rnd, c = rc.split('_')
    try: rows.append((c, rnd, parse(log)))
    except Exception as e: print('skip', rc, e)
rows.sort()
print('%-2s %-3s %6s %6s %8s %8s %7s %7s %7s %7s %7s %7s %6s %6s %9s' % (
    'c', 'rnd', 'steps', 't_end', 'total s', 'ms/step', 'makf', 'hmhz', 'pres',
    'axhm', 'usbc', 'uchk', 'vit', 'pit', 'Cd'))
for c, rnd, r in rows:
    print('%-2s %-3s %6d %6.1f %8.1f %8.3f %7.1f %7.1f %7.1f %7.1f %7.2f %7.2f %6.2f %6.2f %9.5f' % (
        c, rnd, r['steps'], r['tend'], r['total'], 1e3*r['perstep'], tot(r, 'makf'),
        tot(r, 'hmhz'), tot(r, 'pres'), tot(r, 'axhm'), tot(r, 'usbc'), tot(r, 'uchk'),
        r['vit'], r['pit'], r['cd']))
print('\nper call: axhm matvec [us], uchk per step [us]')
for c, rnd, r in rows:
    ax = r.get('axhm', [0, 0, 0])
    print('%-2s %-3s  axhm %8.1f us x %7d   uchk %8.1f us/step   %s' % (
        c, rnd, 1e6*ax[1]/max(ax[0], 1), ax[0], 1e6*tot(r, 'uchk')/r['steps'], CASE[c]))
