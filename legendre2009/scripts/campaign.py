#!/usr/bin/env python3
"""campaign.py -- run matrix for Legendre, Lauga & Magnaudet (2009) figures 1
and 2, set-up of the run directories, and job lists split over the hosts.

   python3 scripts/campaign.py setup     (on a compute host, project root)
   xargs -P <slots> -I{} bash -c '{}' < jobs_<host>.txt     (on each host)

Families (runs/L<Re>_kn<Kn>, Kn = -1 is the shear-free wall):
  fig 2   Re = 50, 100, 200, 500, 800 at Legendre's Kn values, plus 0 and inf
  fig 1   separation boundary: two Kn around Legendre's Kn_sep(Re) for
          12 Re, and Re = 5.5, 7 at Kn = 0 (onset of separation, no slip);
          the shedding boundary comes from the fig 2 runs (lift growth rate)
  check   Re = 800, Kn = 0 and 0.2 on the refined mesh (Lfine800_kn*)
All on the production mesh (1040 elements, N = 7) unless stated.
"""
import os, re, sys, subprocess
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HOSTS = {'host1': 8, 'host2': 11}                # 4-rank job slots per host
KNC = {50: 0.08, 100: 0.36, 200: 0.47, 500: 0.50, 800: 0.52}   # Legendre fig 1
FIG2 = {50:  [0.01, 0.02, 0.05, 0.07, 0.08, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20],
        100: [0.01, 0.02, 0.2, 1, 2, 2.4, 5, 10, 20],   # 0-0.5, inf in runs/re100_kn*; 0.2 again for the fig 1 inset
        200: [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 1, 2, 5, 10],
        500: [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 1, 2, 5, 10, 20],
        800: [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 1, 2, 5]}
SEP = {10: 0.13, 15: 0.28, 20: 0.46, 30: 0.80, 40: 1.10, 65: 1.60, 100: 2.10,
       150: 2.41, 200: 2.45, 300: 2.43, 500: 2.25, 800: 2.05}         # Legendre Kn_sep(Re)
FINE = '--ntheta 56 --nr 34 --dr1 0.05'

def runs():
    R = []
    for Re, kns in FIG2.items():
        extra = [] if Re == 100 else [0.0, -1.0]
        for kn in extra + kns:
            R.append((Re, kn, ''))
    for Re, ks in SEP.items():
        for f in (0.85, 1.15):
            R.append((Re, round(ks*f, 3), ''))
        if Re in (10, 15, 20, 30, 40, 65, 150, 300):      # added: a second separated point
            R.append((Re, round(ks*0.7, 3), ''))
    for Re in (5.5, 7, 45):          # no slip: onset of separation (5.5, 7) and of shedding (45)
        R.append((Re, 0.0, ''))
    for kn in (0.0, 0.2):
        R.append((800, kn, FINE))
    seen, out = set(), []
    for r in R:
        if (r[0], r[1], r[2]) not in seen:
            seen.add((r[0], r[1], r[2])); out.append(r)
    return out

def settings(Re, kn):
    """dt, endTime, steady tolerance, diagnostics interval"""
    dt = 5e-3 if Re <= 200 else 2.5e-3
    knc = KNC.get(Re, 0.0 if Re < 47 else 0.52)
    unsteady = kn >= 0 and Re >= 47 and kn < 1.3*knc
    if not unsteady:
        return dt, 400.0, (1e-6 if Re <= 200 else 1e-5), 100
    near = kn > 0.7*knc
    T = {50: 300, 100: 300, 200: 250}.get(Re, 200)
    if near: T = {50: 600}.get(Re, 400)
    return dt, float(T), 0.0, 10

def name(Re, kn, mesh):
    return '%s%g_kn%g' % ('Lfine' if mesh else 'L', Re, kn)

def setup():
    jobs = []
    for Re, kn, mesh in runs():
        d = os.path.join(ROOT, 'runs', name(Re, kn, mesh))
        if os.path.exists(os.path.join(d, 'run.log')) and 'run successful' in open(os.path.join(d, 'run.log')).read():
            continue                                             # already done
        subprocess.run(['bash', os.path.join(ROOT, 'scripts', 'setup_run.sh'), d, repr(kn)] + mesh.split(),
                       check=True, stdout=subprocess.DEVNULL)
        dt, T, tol, npr = settings(Re, kn)
        p = os.path.join(d, 'cyl.par'); s = open(p).read()
        for key, val in (('viscosity', '%g' % -Re), ('endTime', '%g' % T), ('dt', '%g' % dt),
                         ('userParam02', '%g' % tol), ('userParam03', '%d' % npr)):
            s = re.sub(r'(?m)^%s *=.*$' % key, '%s = %s' % (key, val), s)
        s = re.sub(r'(?m)^(userParam03 *=.*)$', r'\1\nuserParam06 = 0.1\nuserParam07 = 1', s)
        open(p, 'w').write(s)
        for f in ('cyl_drag.dat', 'run.log'):
            if os.path.exists(os.path.join(d, f)): os.remove(os.path.join(d, f))
        nel = 1904 if mesh else 1040
        cost = T/dt*nel/1040.0*(0.5 if tol > 0 else 1.0)        # steady runs usually stop early
        jobs.append((cost, 'cd %s && OMP_NUM_THREADS=1 HWLOC_COMPONENTS=-gl mpirun --bind-to none '
                           '-np 4 ./nek5000 > run.log 2>&1 < /dev/null' % d))
    jobs.sort(reverse=True)                                      # longest first
    load = {h: 0.0 for h in HOSTS}; lists = {h: [] for h in HOSTS}
    for c, j in jobs:
        h = min(HOSTS, key=lambda h: load[h]/HOSTS[h]); load[h] += c; lists[h].append(j)
    for h in HOSTS:
        open(os.path.join(ROOT, 'jobs_%s.txt' % h), 'w').write('\n'.join(lists[h]) + '\n')
        print('%-6s %3d jobs, %2d slots, estimated %.0f min' % (h, len(lists[h]), HOSTS[h], load[h]/HOSTS[h]/40/60))
    print(len(jobs), 'runs set up')

if __name__ == '__main__':
    if sys.argv[1:] == ['setup']:
        setup()
    else:
        for r in runs():
            print(name(*r), settings(r[0], r[1]))
