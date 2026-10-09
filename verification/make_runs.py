#!/usr/bin/env python3
"""make_runs.py [NAME_PREFIX ...] -- set up the verification runs in
verification/runs/<name> and write verification/jobs.txt (one shell command
per run).  Needs NEK5000=/path/to/Nek5000 (gencon, genmap) and the executables
from build_verif.sh (build_N<N>, build_N<N>_3d).  Then:
    xargs -P 16 -I{} bash -c '{}' < verification/jobs.txt
Families (see verif.usr for the three exact solutions):
  p1_kn<Kn>_N<N>     case 1 (Stokes, translating outer cylinder), p-ref.
  p3_kn1_N<N>        case 3 (circular Couette, steady NS), p-refinement
  c3_kn<Kn>_N7       case 3 at N = 7 over a Kn sweep (velocity profiles)
  h1_kn<Kn>_nt<n>    case 1, h-refinement at N = 4
  t2_<scheme>_dt<m>  case 2 (decaying mode, NS), dt = 0.05/2^m at N = 9
  v3c3_kn<Kn>_N7     3D (4 periodic layers, Lz = 1): case 3, circular Couette
  v3c4_kn<Kn>_N7     3D: case 4, axial Couette (Robin term on w), Kn sweep
  v3p4_kn1_N<N>      3D: case 4 at Kn = 1, p-refinement
ALWAYS pass a name filter (e.g. v3): without one every run is set up again
and its old results are deleted.
Annulus a = 0.5, R = 1.5; coarse mesh 8 x 2 elements for p-refinement.
"""
import os, sys, shutil, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NEK  = os.environ.get('NEK5000') or sys.exit('set NEK5000=/path/to/Nek5000')
MESH = os.path.join(ROOT, 'cylinder2d', 'gen_cyl_re2.py')
A, R = 0.5, 1.5
only = sys.argv[1:]            # optional: only set up runs whose name starts with these

def par(case, kn, dt, nsteps, nu, ts='bdf2', ext='standard', nprint=None,
        vtol=1e-12, ptol=1e-11):
    L = ['[GENERAL]', 'stopAt = numSteps', 'numSteps = %d' % nsteps,
         'dt = %.17g' % dt, 'timeStepper = %s' % ts]
    if case != 1:
        L.append('extrapolation = %s' % ext)
        if ext.upper() == 'OIFS':
            L.append('targetCFL = 2.0')
    L += ['writeControl = timeStep', 'writeInterval = 100000000',
          'dealiasing = %s' % ('no' if case == 1 else 'yes'),
          'userParam01 = %g' % kn,
          'userParam03 = %d' % (nprint or max(nsteps // 10, 1)),
          'userParam04 = %d' % case,
          'userParam05 = %g' % (0.1 if case == 2 else 1.0),
          'userParam06 = 0.5', 'userParam07 = 1.0', '',
          '[PROBLEMTYPE]',
          'equation = %s' % ('stokes' if case == 1 else 'incompNS'),
          'stressFormulation = yes', '',
          '[PRESSURE]', 'residualTol = %.1e' % ptol, 'residualProj = no', '',
          '[VELOCITY]', 'residualTol = %.1e' % vtol, 'residualProj = no',
          'density = 1.0', 'viscosity = %g' % nu, '']
    return '\n'.join(L)

jobs = []
def make(name, N, parstr, ntheta=8, nr=2, nproc=1, nz=0, lz=1.0):
    if only and not any(name.startswith(o) for o in only):
        return
    d = os.path.join(HERE, 'runs', name)
    os.makedirs(d, exist_ok=True)
    for f in ('verif_summary.dat', 'verif.dat', 'run.log'):
        if os.path.exists(os.path.join(d, f)):
            os.remove(os.path.join(d, f))
    with open(os.path.join(d, 'mesh.log'), 'w') as log:
        subprocess.run([sys.executable, MESH,
                        '-o', 'verif', '--rinf', str(R), '--ntheta', str(ntheta),
                        '--nr', str(nr), '--dr1', repr((R - A)/nr)]
                       + (['--nz', str(nz), '--lz', repr(lz)] if nz else []),
                       cwd=d, check=True, stdout=log)
    for tool in ('gencon', 'genmap'):
        subprocess.run([os.path.join(NEK, 'bin', tool)], input='verif\n0.2\n', text=True,
                       cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    assert os.path.exists(os.path.join(d, 'verif.ma2')), name + ': genmap failed'
    shutil.copy(os.path.join(HERE, 'build_N%d%s' % (N, '_3d' if nz else ''), 'nek5000'), d)
    open(os.path.join(d, 'verif.par'), 'w').write(parstr)
    open(os.path.join(d, 'SESSION.NAME'), 'w').write('verif\n%s/\n' % d)
    jobs.append('cd %s && OMP_NUM_THREADS=1 HWLOC_COMPONENTS=-gl mpirun '
                '--bind-to none -np %d ./nek5000 > run.log 2>&1 < /dev/null'
                % (d, nproc))

# Case 1 is marched in time to its steady state.  The PN/PN-2 time stepping has
# a slow numerical transient in this closed annulus (error halving about every
# time unit at dt = 1e-2, faster at smaller dt): integrate to t = 20 at dt = 1e-3.
DT1, NS1 = 1e-3, 20000
for N in range(3, 14):                                   # p-refinement
    for kn in (0.0, 0.1, 1.0, 10.0, -1.0):
        make('p1_kn%g_N%d' % (kn, N), N, par(1, kn, DT1, NS1, 1.0))
    make('p3_kn1_N%d' % N, N, par(3, 1.0, 2e-3, 2500, 1.0))
for kn in (0.0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, -1.0):
    make('c3_kn%g_N7' % kn, 7, par(3, kn, 2e-3, 2500, 1.0))   # Couette Kn sweep
for nt, nr in ((8, 2), (16, 4), (32, 8), (64, 16)):     # h-refinement, N = 4
    for kn in (0.0, 0.1, 1.0, -1.0):
        make('h1_kn%g_nt%d' % (kn, nt), 4, par(1, kn, DT1, NS1, 1.0),
             nt, nr, nproc=1 if nt < 64 else 4)
for tag, ts, ext in (('bdf2oifs', 'bdf2', 'OIFS'),       # dt-refinement, N = 9
                     ('bdf2', 'bdf2', 'standard'), ('bdf3', 'bdf3', 'standard')):
    for m in range(6):
        dt = 0.05/2**m
        ns = int(round(1.0/dt))
        make('t2_%s_dt%d' % (tag, m), 9, par(2, 1.0, dt, ns, 0.1, ts, ext, nprint=ns))
# 3D: the same annulus extruded in 4 periodic layers (Lz = 1); the solutions are
# z-invariant, case 4 exercises the Robin term on the spanwise velocity
for kn in (0.0, 0.1, 1.0, 10.0, -1.0):
    make('v3c3_kn%g_N7' % kn, 7, par(3, kn, 2e-3, 2500, 1.0), nz=4, lz=1.0, nproc=4)
for kn in (0.0, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, -1.0):
    make('v3c4_kn%g_N7' % kn, 7, par(4, kn, 2e-3, 2500, 1.0), nz=4, lz=1.0, nproc=4)
for N in (3, 5, 7, 9):
    make('v3p4_kn1_N%d' % N, N, par(4, 1.0, 2e-3, 2500, 1.0), nz=4, lz=1.0, nproc=4)
open(os.path.join(HERE, 'jobs.txt'), 'w').write('\n'.join(jobs) + '\n')
print(len(jobs), 'runs set up; job list in', os.path.join(HERE, 'jobs.txt'))
