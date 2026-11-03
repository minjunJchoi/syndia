#!/usr/bin/env python3
"""
Regression test for the compatibility changes (compat.py, NumPy >= 2 / SciPy >= 1.14 / Python >= 3.12, lazy MDSplus,
offline ECEI parameters).  Runs OFFLINE (no MDSplus needed) on the shipped example: KSTAR #37389, 4150 ms, ECEI GT.

Checks
  1. the result reproduces data/..._REF.pkl (the original code's result)
  2. parallel run == sequential run
Results are saved to regression_np<ver>_sp<ver>.npz so runs on different NumPy/SciPy versions can be compared:
    python regression_test.py --compare a.npz b.npz

usage:  python regression_test.py [--quick]      (--quick: every 8th ECEI channel)
"""
import os, sys, io, time, pickle, contextlib
sys.path.insert(0, os.pardir)
import numpy as np
import scipy

HERE = os.path.dirname(os.path.abspath(__file__))
ECEI_PARAMS = {'hn': 1, 'lo': 77.0, 'sf': 360.0, 'sz': 150.0}   # MDSplus values for #37389 (GT)
RUN = dict(fstart=0, fend=0, Nf=1, zstart=0, zend=0, Nz=1, pstart=7.8, pend=-2, pint=0.1,
           Rinit=2.39, torbeam=0, select='max')


def run(efm, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        out = efm.run(**RUN, **kw)
    return [np.asarray(v, dtype=float) for v in out]


# Run only as a script: parallel runs use multiprocessing, and with the 'spawn' start method
# (Windows, macOS) the module is re-imported in each worker.
if __name__ == '__main__':
    os.chdir(HERE)
    if len(sys.argv) > 1 and sys.argv[1] == '--compare':
        a, b = np.load(sys.argv[2]), np.load(sys.argv[3])
        worst = max(np.nanmax(np.abs(a[k] - b[k])/np.maximum(np.abs(a[k]), 1e-30)) for k in a.files)
        print('RESULT:', 'IDENTICAL within 1e-6' if worst < 1e-6 else f'max relative difference {worst:.2e}')
        sys.exit(0)
    from efm import EceFwdMod
    ok_all = True
    print(f'Python {sys.version.split()[0]}, NumPy {np.__version__}, SciPy {scipy.__version__}')
    t0 = time.time()
    efm = EceFwdMod()
    with contextlib.redirect_stdout(io.StringIO()):
        efm.set_profile('data/KSTAR_EFIT01_g037389.004150', 'data/ece_Te_37389_4150ms.dat',
                        'data/tci_ne_37389_4150ms.dat', bfactor=1.01)
        efm.set_channel(37389, ['ECEI_GT0101-2408'], ecei_params=ECEI_PARAMS)
    ref = pickle.load(open('data/ecei_pos_37389_ECEI_GT0101-2408_4150ms_b1.01_max_REF.pkl', 'rb'))
    Rr, zr, taur, Trr = (np.asarray(x, float) for x in ref[1:5])
    if '--quick' in sys.argv:
        sel = np.arange(0, len(Rr), 8)
        with contextlib.redirect_stdout(io.StringIO()):
            efm.diag.clist = [efm.diag.clist[i] for i in sel]
            efm.diag.vn_list, efm.diag.fn_list, efm.diag.ece_freq = efm.diag.vn_list[sel], efm.diag.fn_list[sel], efm.diag.ece_freq[sel]
        Rr, zr, taur, Trr = Rr[sel], zr[sel], taur[sel], Trr[sel]
    S = run(efm, parallel=False)
    dR, dz = np.max(np.abs(S[0] - Rr)), np.max(np.abs(S[1] - zr))
    dtau, dT = np.max(np.abs(S[3] - taur)/taur), np.max(np.abs(S[4] - Trr)/Trr)
    ok1 = dR < 1e-5 and dz < 1e-5 and dtau < 1e-3 and dT < 1e-4
    print(f"[{'PASS' if ok1 else 'FAIL'}] reproduces _REF.pkl: max|dR| {dR:.1e} m, max|dz| {dz:.1e} m, "
          f"max rel dtau {dtau:.1e}, max rel dTrad {dT:.1e}")
    P = run(efm, parallel=True, n_processes=2)
    dpar = max(np.max(np.abs(P[k] - S[k])) for k in (0, 1, 4))
    ok2 = dpar < 1e-12
    print(f"[{'PASS' if ok2 else 'FAIL'}] parallel == sequential: max difference {dpar:.1e}")
    ok_all = ok1 and ok2
    fn = f'regression_np{np.__version__}_sp{scipy.__version__}.npz'
    np.savez(fn, R=S[0], z=S[1], tau=S[3], Trad=S[4])
    print(f'saved {fn}   ({time.time()-t0:.0f} s)')
    print('OVERALL:', 'PASS' if ok_all else 'FAIL')
    sys.exit(0 if ok_all else 1)
