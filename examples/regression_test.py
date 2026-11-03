#!/usr/bin/env python3
"""
Regression test for the compatibility layer and the line-shape / absorption-prefactor options
(shape_model, amp_model).  Runs OFFLINE (no MDSplus needed) on the shipped example: KSTAR #37389, 4150 ms, ECEI GT.

Checks
  1. shape_model='legacy', amp_model='legacy' reproduces data/..._REF.pkl (the original code's result)
  2. parallel run == sequential run
  3. 'general' line shape is normalised, Int phi dw = 2 pi ('legacy' is reported)
  4. ECE radiometer (X-mode, m = 2): 'general' == 'legacy' line shape
  5. O1 absolute absorption (100 eV, 90 deg, low density): Int alpha dw = (pi/2) wp^2 Theta / c with amp_model='corrected'
     (Bekefi single-particle emission + Kirchhoff); 'legacy' gives 2 pi times this
  6. X-mode m = 2 at 90 deg: 'corrected' == 'legacy' (the Appleton-Hartree change only matters at oblique angles)
and reports the change of the ECEI O-mode result for the new defaults.
Results are saved to regression_np<ver>_sp<ver>.npz (compare runs with  python regression_test.py --compare a.npz b.npz).

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
ok_all = True


def check(name, cond, detail):
    global ok_all
    ok_all &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")


def run(efm, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        out = efm.run(**RUN, **kw)
    return [np.asarray(v, dtype=float) for v in out]


def int_alpha(eceint, m, TeV, ne, thd, amp, B=2.5):
    """Int alpha domega at one point (constant B, T_e, n_e) from eceint.ece_intensity."""
    e_, me_ = 1.602e-19, 9.109e-31
    wce = e_*B/me_
    s = np.array([0.0, 0.01, 0.02]); th = np.full(3, np.radians(thd))
    T = TeV*e_
    xs = np.linspace(m*(1 - 0.03), m*(1 + 0.005), 1201)
    a = [eceint.ece_intensity(s, s, 0*s, th, x*wce, m, lambda q: B + 0*q, lambda q: T + 0*q, lambda q: ne + 0*q,
                              select='max', use_cache=False, amp_model=amp)[-1][1] for x in xs]
    return np.trapezoid(a, xs*wce) if hasattr(np, 'trapezoid') else np.trapz(a, xs*wce)


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
    import eceint
    from compat import trapezoid
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
        sel = np.arange(0, len(efm.diag.clist), 8)
        efm.diag.clist = [efm.diag.clist[i] for i in sel]
        efm.diag.vn_list, efm.diag.fn_list, efm.diag.ece_freq = efm.diag.vn_list[sel], efm.diag.fn_list[sel], efm.diag.ece_freq[sel]
        Rr, zr, taur, Trr = Rr[sel], zr[sel], taur[sel], Trr[sel]

    # 1. legacy reproduces the original result
    L = run(efm, parallel=False, shape_model='legacy', amp_model='legacy')
    dR, dz = np.max(np.abs(L[0] - Rr)), np.max(np.abs(L[1] - zr))
    dtau, dT = np.max(np.abs(L[3] - taur)/taur), np.max(np.abs(L[4] - Trr)/Trr)
    check('legacy reproduces _REF.pkl', dR < 1e-5 and dz < 1e-5 and dtau < 1e-3 and dT < 1e-4,
          f'max|dR| {dR:.1e} m, max|dz| {dz:.1e} m, max rel dtau {dtau:.1e}, max rel dTrad {dT:.1e}')

    # 2. parallel == sequential (defaults)
    G = run(efm, parallel=False)
    P = run(efm, parallel=True, n_processes=2)
    dpar = max(np.max(np.abs(P[k] - G[k])) for k in (0, 1, 4))
    check('parallel == sequential', dpar < 1e-12, f'max difference {dpar:.1e}')

    # 3. line-shape normalisation
    e_, me_ = 1.602e-19, 9.109e-31
    wce = e_*2.5/me_
    norms = {}
    for mdl in ('general', 'legacy'):
        worst = 0.0
        for m, TeV, thd in [(1, 100, 60), (1, 300, 90), (1, 2000, 80), (2, 50, 60), (2, 100, 60), (2, 300, 60), (2, 2000, 80), (2, 1500, 88)]:
            zeta = eceint.mc2/(TeV*e_)/2
            th = np.radians(thd)
            vth = 1/np.sqrt(zeta)
            xs = np.linspace(1 - 10*vth*abs(np.cos(th)) - 8/zeta, min(1 + 10*vth*abs(np.cos(th)) + 1e-3, 1/np.sin(th) - 1e-9), 2001)
            ph = [eceint.shape_fn(m, zeta, wce, x*m*wce, th, mdl) for x in xs]
            worst = max(worst, abs(trapezoid(ph, x=xs*m*wce)/(2*np.pi) - 1))
        norms[mdl] = worst
    check("'general' line shape normalised (Int phi dw = 2 pi)", norms['general'] < 1e-3,
          f"max |norm-1| = {norms['general']:.1e}   (legacy: {norms['legacy']:.2f}; scipy quad defaults at <= 300 eV, 60 deg)")

    # 4. ECE radiometer, m = 2: general == legacy line shape
    ece = EceFwdMod()
    with contextlib.redirect_stdout(io.StringIO()):
        ece.set_profile('data/KSTAR_EFIT01_g037389.004150', 'data/ece_Te_37389_4150ms.dat',
                        'data/tci_ne_37389_4150ms.dat', bfactor=1.01)
        ece.set_channel(37389, ['ECE{:02d}'.format(i) for i in range(25, 49)])
    EL = run(ece, parallel=False, shape_model='legacy', amp_model='legacy')
    EG = run(ece, parallel=False, shape_model='general', amp_model='legacy')
    ok = np.isfinite(EL[0]) & (EL[2] > 0)
    d4 = max(np.nanmax(np.abs(EG[0][ok] - EL[0][ok])), np.nanmax(np.abs(EG[4][ok] - EL[4][ok])/EL[4][ok]))
    check("ECE X-mode m=2: 'general' == 'legacy' line shape", d4 < 1e-6, f'{ok.sum()} channels, max difference {d4:.1e}')

    # 5. O1 absolute absorption
    ref5 = 0.5*np.pi*1e17*e_**2/(eceint.eps*me_)/eceint.c*(100*e_/eceint.mc2)
    r5c, r5l = int_alpha(eceint, 1, 100, 1e17, 90, 'corrected')/ref5, int_alpha(eceint, 1, 100, 1e17, 90, 'legacy')/ref5
    check('O1 Int alpha dw = (pi/2) wp^2 Theta/c (100 eV, 90 deg)', abs(r5c - 1) < 0.01,
          f"corrected {r5c:.4f}, legacy {r5l:.3f} (= 2 pi x {r5l/(2*np.pi):.4f})")

    # 6. X2 at 90 deg: corrected == legacy; oblique ratio reported
    x90 = int_alpha(eceint, 2, 1000, 1e19, 90, 'corrected', B=2.0)/int_alpha(eceint, 2, 1000, 1e19, 90, 'legacy', B=2.0)
    x80 = int_alpha(eceint, 2, 1000, 1e19, 80, 'corrected', B=2.0)/int_alpha(eceint, 2, 1000, 1e19, 80, 'legacy', B=2.0)
    check("X2 at 90 deg: amp 'corrected' == 'legacy'", abs(x90 - 1) < 1e-9, f'ratio {x90:.12f}; at 80 deg corrected/legacy = {x80:.3f}')

    # report: ECEI O-mode change for the defaults
    dRg = 1e3*(G[0] - L[0])
    print(f"\nECEI O-mode, defaults - legacy: tau {L[3].min():.1f}..{L[3].max():.1f} -> {G[3].min():.1f}..{G[3].max():.1f}, "
          f"dR = {dRg.min():+.2f} .. {dRg.max():+.2f} mm, dTrad = {100*np.min(G[4]/L[4]-1):+.2f} .. {100*np.max(G[4]/L[4]-1):+.2f} %")

    fn = f'regression_np{np.__version__}_sp{scipy.__version__}.npz'
    np.savez(fn, ecei_legacy_R=L[0], ecei_legacy_tau=L[3], ecei_legacy_Trad=L[4],
             ecei_default_R=G[0], ecei_default_tau=G[3], ecei_default_Trad=G[4], ece_legacy_R=EL[0], ece_legacy_Trad=EL[4])
    print(f'\nsaved {fn}   ({time.time()-t0:.0f} s)')
    print('OVERALL:', 'PASS' if ok_all else 'FAIL')
    sys.exit(0 if ok_all else 1)
