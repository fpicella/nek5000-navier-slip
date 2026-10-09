#!/usr/bin/env python3
"""plot_verif.py [RUNS_DIR] -- tabulate and plot the verification suite
(RUNS_DIR/*/verif_summary.dat, default data/ = the results shipped with the
repository; pass runs/ after running make_runs.py; see make_runs.py for the
run families).  Writes verif_results.txt and verif_convergence.png."""
import glob, os, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
runs = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'data')
S = {}
for f in glob.glob(os.path.join(runs, '*', 'verif_summary.dat')):
    S[os.path.basename(os.path.dirname(f))] = np.atleast_2d(np.loadtxt(f))[0]
# columns: 0 case 1 N 2 nel 3 Kn 4 dt 5 nsteps 6 err_inf 7 err_L2 8 err_wall
#          9 relQ 10 res_navier 11 sigma 12 sigma_exact
out = []
def P(s=''):
    print(s); out.append(s)
lab = {0.0: 'no-slip', 0.1: 'Kn = 0.1', 1.0: 'Kn = 1', 10.0: 'Kn = 10', -1.0: 'shear-free'}
fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))

# (a,b) p-refinement on 16 elements
P('p-refinement, 8 x 2 elements, annulus a = 0.5, R = 1.5')
P('case 1 = Stokes, outer cylinder translating (slip varies along the wall)')
P('   Kn       N   err_inf     err_L2      err_wall    relF        res_navier')
for i, kn in enumerate((0.0, 0.1, 1.0, 10.0, -1.0)):
    Ns, e, q, r = [], [], [], []
    for N in range(3, 14):
        k = 'p1_kn%g_N%d' % (kn, N)
        if k in S:
            s = S[k]; Ns.append(N); e.append(s[6]); q.append(s[9]); r.append(s[10])
            P('%6g  %4d  %.3e  %.3e  %.3e  %.3e  %.3e' % (kn, N, s[6], s[7], s[8], s[9], s[10]))
    if Ns:
        ax[0, 0].semilogy(Ns, e, 'o-', color='C%d' % i, label='case 1, ' + lab[kn])
        ax[0, 1].semilogy(Ns, q, 'o-', color='C%d' % i, label='force, ' + lab[kn])
Ns, e, q = [], [], []
P('case 3 = circular Couette, Kn = 1 (steady Navier-Stokes)')
for N in range(3, 14):
    k = 'p3_kn1_N%d' % N
    if k in S:
        s = S[k]; Ns.append(N); e.append(s[6]); q.append(s[9])
        P('%6g  %4d  %.3e  %.3e  %.3e  %.3e  %.3e' % (1, N, s[6], s[7], s[8], s[9], s[10]))
if Ns:
    ax[0, 0].semilogy(Ns, e, 'ks--', mfc='none', label='case 3 (Couette), Kn = 1')
    ax[0, 1].semilogy(Ns, q, 'ks--', mfc='none', label='torque, Couette, Kn = 1')
tt = [(k, S[k]) for k in sorted(S) if k.startswith('p1tight_')]
if tt:
    P('case 1 with solver tolerances 1e-14 (instead of 1e-12 velocity, 1e-11 pressure)')
    for k, s in tt:
        P('%6g  %4d  %.3e  %.3e  %.3e  %.3e  %.3e' % (s[3], s[1], s[6], s[7], s[8], s[9], s[10]))
    for kn, mk in ((0.0, 'x'), (1.0, '+')):
        z = [(s[1], s[6], s[9]) for k, s in tt if s[3] == kn]
        if z:
            z = np.array(z)
            ax[0, 0].semilogy(z[:, 0], z[:, 1], 'k' + mk, ms=10, label='tolerance 1e-14, ' + lab[kn])
            ax[0, 1].semilogy(z[:, 0], z[:, 2], 'k' + mk, ms=10, label='tolerance 1e-14, ' + lab[kn])
ax[0, 0].set(xlabel='polynomial order N', ylabel='max |u - u_exact|', title='(a) p-refinement: velocity error')
ax[0, 1].set(xlabel='polynomial order N', ylabel='relative error', title='(b) p-refinement: force / torque on the cylinder')
for a in ax[0]: a.grid(alpha=.3); a.legend(fontsize=7)

# (c) h-refinement at N = 4
P(''); P('h-refinement at N = 4, case 1 (element size h = (R-a)/nr, ntheta = 4 nr)')
P('   Kn   ntheta  err_inf     order   err_L2      order   relF        order')
for i, kn in enumerate((0.0, 0.1, 1.0, -1.0)):
    h, e, e2, q = [], [], [], []
    for nt in (8, 16, 32, 64):
        k = 'h1_kn%g_nt%d' % (kn, nt)
        if k in S:
            h.append(4.0/nt); e.append(S[k][6]); e2.append(S[k][7]); q.append(S[k][9])
    for j in range(len(h)):
        o = lambda v: ('%6.2f' % (np.log(v[j-1]/v[j])/np.log(h[j-1]/h[j]))) if j else '     -'
        P('%6g  %5d  %.3e  %s   %.3e  %s   %.3e  %s' % (kn, int(round(4/h[j])), e[j], o(e), e2[j], o(e2), q[j], o(q)))
    if h:
        c = 'C%d' % {0.0: 0, 0.1: 1, 1.0: 2, -1.0: 4}[kn]
        ax[1, 0].loglog(h, e, 'o-', color=c, label='max error, ' + lab[kn])
        ax[1, 0].loglog(h, e2, 's--', color=c, mfc='none', label='L2 error, ' + lab[kn])
if S.get('h1_kn1_nt8') is not None:
    hh = np.array([0.5, 0.0625])
    ax[1, 0].loglog(hh, 0.5*S['h1_kn1_nt8'][7]*(hh/0.5)**4, 'k:', label='slope 4 = N')
    ax[1, 0].loglog(hh, 2.0*S['h1_kn1_nt8'][6]*(hh/0.5)**3, 'k-.', label='slope 3 = N-1')
ax[1, 0].set(xlabel='element size h', ylabel='velocity error', title='(c) h-refinement at N = 4 (case 1)')
ax[1, 0].grid(alpha=.3, which='both'); ax[1, 0].legend(fontsize=7)

# (d) dt-refinement of the decay rate, case 2
P(''); P('dt-refinement, case 2 (decaying azimuthal mode, Navier-Stokes), Kn = 1, N = 9')
P('   scheme       dt          sigma_num             rel.err     order   err_inf(T)')
for i, (tag, name) in enumerate((('bdf2oifs', 'BDF2 + OIFS (production)'), ('bdf2', 'BDF2 / EXT2'), ('bdf3', 'BDF3 / EXT3'))):
    dts, er = [], []
    for m in range(6):
        k = 't2_%s_dt%d' % (tag, m)
        if k in S:
            s = S[k]; dts.append(s[4]); er.append(abs(s[11] - s[12])/s[12])
            o = ('%6.2f' % (np.log(er[-2]/er[-1])/np.log(dts[-2]/dts[-1]))) if len(er) > 1 else '     -'
            P('   %-10s %.5f  %.14f  %.3e  %s   %.3e' % (tag, s[4], s[11], er[-1], o, s[6]))
    if dts:
        ax[1, 1].loglog(dts, er, 'o-', color='C%d' % i, label=name)
k0 = S.get('t2_bdf2_dt0')
if k0 is not None:
    P('   exact sigma = nu k^2 = %.14f' % k0[12])
    dd = np.array([0.05, 0.0015625]); e0 = abs(k0[11] - k0[12])/k0[12]
    ax[1, 1].loglog(dd, e0*(dd/0.05)**2, 'k:', label='slope 2')
k3 = S.get('t2_bdf3_dt0')
if k3 is not None:
    dd = np.array([0.05, 0.0015625]); e3 = abs(k3[11] - k3[12])/k3[12]
    ax[1, 1].loglog(dd, e3*(dd/0.05)**3, 'k--', label='slope 3')
ax[1, 1].set(xlabel='time step dt', ylabel='|sigma - sigma_exact| / sigma_exact', title='(d) dt-refinement: decay rate (case 2)')
ax[1, 1].grid(alpha=.3, which='both'); ax[1, 1].legend(fontsize=7)
fig.suptitle('Navier slip on a curved wall in Nek5000: convergence to exact solutions', fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'verif_convergence.png'), dpi=140)
open(os.path.join(HERE, 'verif_results.txt'), 'w').write('\n'.join(out) + '\n')
