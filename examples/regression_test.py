#!/usr/bin/env python3
"""
Regression test for SYNDIA (compatibility layer, shape_model, amp_model and dist options).
Runs OFFLINE (no MDSplus needed) on the shipped example: KSTAR #37389, 4150 ms, ECEI GT.

Checks
  1. shape_model='legacy', amp_model='legacy' reproduces data/..._REF.pkl  (the original code's result)
  2. parallel run == sequential run
  3. line-shape normalisation  Int phi dw = 2 pi  ('general' must pass; 'legacy' is reported)
  4. ECE radiometer (X-mode, m = 2): 'general' == 'legacy' (the m = 2 integrand is unchanged)
  5. dist='bimaxwell' with t_r = 1 == dist='maxwell' (ECEI O-mode and ECE X-mode)
  6. dist='bimaxwell', t_r = 2: reported T_rad / tau changes; the ECEI O1 result agrees with the kinetic
     kernel for the same bi-Maxwellian (tau ratio within 3 %, T_rad within 1 %), see item 7
  7. dist='general' (kinetic kernel) with a Maxwellian agrees with the corrected Hutchinson model
     (also the ECEI O1 emission position on the midplane row and on an oblique row, and a bi-Maxwellian check)
  8. dist='general': parallel == sequential with a lab-frame Grad13 model (pickling)
and reports the 'general' - 'legacy' change for the ECEI (O-mode, m = 1) channels.
Results are also saved to regression_np<ver>_sp<ver>.npz so runs on different NumPy/SciPy
versions can be compared with   python regression_test.py --compare a.npz b.npz

usage:  python regression_test.py [--quick]      (--quick: every 8th ECEI channel)
"""
import os, sys, io, time, pickle, contextlib
sys.path.insert(0, os.pardir)
import numpy as np
import scipy

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

# Run only as a script: parallel runs use multiprocessing, and with the 'spawn' start method
# (Windows, macOS) the module is re-imported in each worker.
if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--compare':
        a, b = np.load(sys.argv[2]), np.load(sys.argv[3])
        worst = 0.0
        for k in a.files:
            if a[k].dtype.kind == 'f':
                d = np.nanmax(np.abs(a[k] - b[k]) / np.maximum(np.abs(a[k]), 1e-30))
                worst = max(worst, d)
                print(f'{k:28s} max relative difference {d:.2e}')
        print('RESULT:', 'IDENTICAL within 1e-6' if worst < 1e-6 else f'max relative difference {worst:.2e}')
        sys.exit(0)

    from efm import EceFwdMod
    import eceint
    from compat import trapezoid

    QUICK = '--quick' in sys.argv
    ECEI_PARAMS = {'hn': 1, 'lo': 77.0, 'sf': 360.0, 'sz': 150.0}   # MDSplus values for #37389 (GT)
    RUN = dict(fstart=0, fend=0, Nf=1, zstart=0, zend=0, Nz=1, pstart=7.8, pend=-2, pint=0.1,
               Rinit=2.39, torbeam=0, select='max')
    results, ok_all = {}, True


    def check(name, cond, detail):
        global ok_all
        ok_all &= bool(cond)
        print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")


    def run(efm, **kw):
        with contextlib.redirect_stdout(io.StringIO()):
            out = efm.run(**RUN, **kw)
        return [np.asarray(x, dtype=float) for x in out]   # Rch, zch, int_meas, tau, rad_temp, abs_temp


    print(f'Python {sys.version.split()[0]}, NumPy {np.__version__}, SciPy {scipy.__version__}')
    t0 = time.time()
    efm = EceFwdMod()
    with contextlib.redirect_stdout(io.StringIO()):
        efm.set_profile('data/KSTAR_EFIT01_g037389.004150', 'data/ece_Te_37389_4150ms.dat',
                        'data/tci_ne_37389_4150ms.dat', bfactor=1.01)
        efm.set_channel(37389, ['ECEI_GT0101-2408'], ecei_params=ECEI_PARAMS)
    ref = pickle.load(open('data/ecei_pos_37389_ECEI_GT0101-2408_4150ms_b1.01_max_REF.pkl', 'rb'))
    Rr, zr, taur, Trr = (np.asarray(x, float) for x in ref[1:5])
    if QUICK:
        sel = np.arange(0, len(efm.diag.clist), 8)
        efm.diag.clist = [efm.diag.clist[i] for i in sel]
        efm.diag.vn_list, efm.diag.fn_list, efm.diag.ece_freq = efm.diag.vn_list[sel], efm.diag.fn_list[sel], efm.diag.ece_freq[sel]
        Rr, zr, taur, Trr = Rr[sel], zr[sel], taur[sel], Trr[sel]

    # 1. legacy vs REF
    L = run(efm, parallel=False, shape_model='legacy', amp_model='legacy')
    dR, dz = np.max(np.abs(L[0] - Rr)), np.max(np.abs(L[1] - zr))
    dtau, dT = np.max(np.abs(L[3] - taur) / taur), np.max(np.abs(L[4] - Trr) / Trr)
    check('legacy reproduces _REF.pkl', dR < 1e-5 and dz < 1e-5 and dtau < 1e-3 and dT < 1e-4,
          f'max|dR| {dR:.1e} m, max|dz| {dz:.1e} m, max rel dtau {dtau:.1e}, max rel dTrad {dT:.1e}')

    # 2. parallel == sequential
    P = run(efm, parallel=True, n_processes=2, shape_model='legacy', amp_model='legacy')
    dpar = max(np.max(np.abs(P[k] - L[k])) for k in (0, 1, 4))
    check('parallel == sequential', dpar < 1e-12, f'max difference {dpar:.1e}')

    # 3. line-shape normalisation
    e_, me_ = 1.602e-19, 9.109e-31
    wce = e_*2.5/me_
    norms = {}
    for mdl in ('general', 'legacy'):
        worst = 0.0
        for m, TeV, thd in [(1, 100, 60), (1, 2000, 80), (2, 50, 60), (2, 100, 60), (2, 300, 60), (2, 2000, 80), (2, 1500, 88)]:
            zeta = eceint.mc2/(TeV*e_)/2
            th = np.radians(thd)
            vth = 1/np.sqrt(zeta)
            xs = np.linspace(1 - 10*vth*abs(np.cos(th)) - 8/zeta, min(1 + 10*vth*abs(np.cos(th)) + 1e-3, 1/np.sin(th) - 1e-9), 2001)
            ph = [eceint.shape_fn(m, zeta, wce, x*m*wce, th, mdl) for x in xs]
            n = trapezoid(ph, x=xs*m*wce) / (2*np.pi)
            worst = max(worst, abs(n - 1))
        norms[mdl] = worst
    check("'general' line shape normalised (Int phi dw = 2 pi)", norms['general'] < 1e-3,
          f"max |norm-1| = {norms['general']:.1e}   (legacy: {norms['legacy']:.2f}, i.e. the known quad problem at <=300 eV, 60 deg)")

    # 4. ECE radiometer, m = 2: general == legacy
    ece = EceFwdMod()
    with contextlib.redirect_stdout(io.StringIO()):
        ece.set_profile('data/KSTAR_EFIT01_g037389.004150', 'data/ece_Te_37389_4150ms.dat',
                        'data/tci_ne_37389_4150ms.dat', bfactor=1.01)
        ece.set_channel(37389, ['ECE{:02d}'.format(i) for i in range(25, 49)])   # 126-162 GHz, inside the plasma
    EL = run(ece, parallel=False, shape_model='legacy', amp_model='legacy')
    EG = run(ece, parallel=False, shape_model='general', amp_model='legacy')
    ok = np.isfinite(EL[0]) & (EL[2] > 0)
    d2 = max(np.nanmax(np.abs(EG[0][ok] - EL[0][ok])), np.nanmax(np.abs(EG[4][ok] - EL[4][ok]) / EL[4][ok]))
    check("ECE X-mode m=2: 'general' == 'legacy'", d2 < 1e-6, f'{ok.sum()} channels, max difference {d2:.1e}')

    # report: general - legacy for ECEI O-mode (m = 1)
    Gs = run(efm, parallel=False, shape_model='general', amp_model='legacy')   # shape fix only
    G = run(efm, parallel=False)                                               # defaults: shape + amplitude fixes
    EG = run(ece, parallel=False)

    # 5. bi-Maxwellian, t_r = 1 == Maxwellian
    with contextlib.redirect_stdout(io.StringIO()):
        efm.set_anisotropy(1.0); ece.set_anisotropy(1.0)
    B1 = run(efm, parallel=False, dist='bimaxwell')
    EB1 = run(ece, parallel=False, dist='bimaxwell')
    d5 = max(max(np.max(np.abs(B1[k] - G[k])/np.abs(G[k])) for k in (0, 1, 3, 4)),
             max(np.nanmax(np.abs(EB1[k][ok] - EG[k][ok])/np.abs(EG[k][ok])) for k in (0, 4)))
    check("dist='bimaxwell', t_r = 1 == 'maxwell'", d5 < 1e-9, f'max relative difference {d5:.1e}')

    # 6. bi-Maxwellian, t_r = 2
    with contextlib.redirect_stdout(io.StringIO()):
        efm.set_anisotropy(2.0); ece.set_anisotropy(2.0)
    B2 = run(efm, parallel=False, dist='bimaxwell')
    EB2 = run(ece, parallel=False, dist='bimaxwell')
    rt = B2[3]/G[3]
    print(f'    bi-Maxwellian t_r = 2: ECEI O-mode tau ratio {rt.min():.3f} .. {rt.max():.3f} '
          '(O1 strength ~ T_par; below 2 where the line is Doppler broadened, checked against the kinetic kernel in 7)')
    print(f"    ECEI O-mode   dR = {1e3*np.min(B2[0]-G[0]):+.2f} .. {1e3*np.max(B2[0]-G[0]):+.2f} mm, "
          f"dTrad = {100*np.min(B2[4]/G[4]-1):+.2f} .. {100*np.max(B2[4]/G[4]-1):+.2f} %")
    print(f"    ECE X-mode    dR = {1e3*np.nanmin(EB2[0][ok]-EG[0][ok]):+.2f} .. {1e3*np.nanmax(EB2[0][ok]-EG[0][ok]):+.2f} mm, "
          f"dTrad = {100*np.nanmin(EB2[4][ok]/EG[4][ok]-1):+.2f} .. {100*np.nanmax(EB2[4][ok]/EG[4][ok]-1):+.2f} %  (largest at the optically thin edge)")
    dRg = 1e3*(Gs[0] - L[0])
    print(f"\nECEI O-mode (m=1), shape_model 'general' - 'legacy': dR = {dRg.min():+.2f} .. {dRg.max():+.2f} mm (mean {dRg.mean():+.2f}), "
          f"dTrad = {100*np.min((Gs[4]-L[4])/L[4]):+.2f} .. {100*np.max((Gs[4]-L[4])/L[4]):+.2f} %, "
          f"dtau = {100*np.min((Gs[3]-L[3])/L[3]):+.1f} .. {100*np.max((Gs[3]-L[3])/L[3]):+.1f} %")
    dRa = 1e3*(G[0] - Gs[0])
    print(f"ECEI O-mode (m=1), amp_model 'corrected' - 'legacy' (O1 prefactor / 2 pi): tau {L[3].min():.1f}..{L[3].max():.1f} -> "
          f"{G[3].min():.1f}..{G[3].max():.1f}, dR = {dRa.min():+.2f} .. {dRa.max():+.2f} mm, "
          f"dTrad = {100*np.min(G[4]/Gs[4]-1):+.2f} .. {100*np.max(G[4]/Gs[4]-1):+.2f} %")

    # 7. kinetic kernel (dist='general') with a Maxwellian vs the corrected Hutchinson model
    #    (ECEI: one vertical row, vn = 12, 8 frequencies; ECE: the same 24 channels)
    from distributions import MaxwellModel, Grad13Model
    efk = EceFwdMod()
    with contextlib.redirect_stdout(io.StringIO()):
        efk.set_profile('data/KSTAR_EFIT01_g037389.004150', 'data/ece_Te_37389_4150ms.dat',
                        'data/tci_ne_37389_4150ms.dat', bfactor=1.01)
        efk.set_channel(37389, ['ECEI_GT1201-1208'], ecei_params=ECEI_PARAMS)
    efk.set_distribution(MaxwellModel()); ece.set_distribution(MaxwellModel())
    Gk = run(efk, parallel=False)
    K = run(efk, parallel=False, dist='general')
    EK = run(ece, parallel=False, dist='general')
    dKR, dKe = 1e3*(K[0] - Gk[0]), 1e3*(EK[0][ok] - EG[0][ok])
    dTi, dTe = K[4]/Gk[4] - 1, EK[4][ok]/EG[4][ok] - 1
    # an off-midplane row (oblique theta), where the O1 coupling weight matters; N = 1 in the resonance as in eceint
    efo = EceFwdMod()
    with contextlib.redirect_stdout(io.StringIO()):
        efo.set_profile('data/KSTAR_EFIT01_g037389.004150', 'data/ece_Te_37389_4150ms.dat',
                        'data/tci_ne_37389_4150ms.dat', bfactor=1.01)
        efo.set_channel(37389, ['ECEI_GT0201-0208'], ecei_params=ECEI_PARAMS)
    efo.set_distribution(MaxwellModel())
    Go, Lo = run(efo, parallel=False), run(efo, parallel=False, shape_model='legacy')
    Ko = run(efo, parallel=False, dist='general', kin_opts={'N_res': 'vacuum'})
    dRo, dRl = 1e3*(Go[0] - Ko[0]), 1e3*(Lo[0] - Ko[0])
    print(f"    ECEI row vn=2 (oblique) vs kinetic (N_res='vacuum'): dR general {dRo.min():+.2f}..{dRo.max():+.2f} mm, "
          f"legacy shape {dRl.min():+.2f}..{dRl.max():+.2f} mm")
    # bi-Maxwellian t_r = 2 on the midplane row: eceint vs kinetic kernel
    from distributions import BiMaxwellModel
    with contextlib.redirect_stdout(io.StringIO()):
        efk.set_anisotropy(2.0)
    Bk = run(efk, parallel=False, dist='bimaxwell')
    Kv = run(efk, parallel=False, dist='general', kin_opts={'N_res': 'vacuum'})
    efk.set_distribution(BiMaxwellModel(2.0))
    KB = run(efk, parallel=False, dist='general', kin_opts={'N_res': 'vacuum'})
    efk.set_distribution(MaxwellModel())
    dtr = (Bk[3]/Gk[3])/(KB[3]/Kv[3]) - 1
    dTb = Bk[4]/KB[4] - 1
    print(f"    bi-Maxwellian t_r = 2, midplane row: tau ratio SYNDIA {np.min(Bk[3]/Gk[3]):.3f}..{np.max(Bk[3]/Gk[3]):.3f}, "
          f"kinetic {np.min(KB[3]/Kv[3]):.3f}..{np.max(KB[3]/Kv[3]):.3f}; T_rad SYNDIA/kinetic - 1 = {100*dTb.min():+.2f}..{100*dTb.max():+.2f} %")
    check("dist='general' (Maxwellian) vs corrected Hutchinson",
          np.max(np.abs(dTe)) < 0.02 and np.max(np.abs(dTi)) < 0.05 and np.max(np.abs(dKR)) < 1.5 and np.max(np.abs(dRo)) < 1.0
          and np.max(np.abs(dtr)) < 0.03 and np.max(np.abs(dTb)) < 0.01,
          f"ECE X2: dTrad {100*dTe.min():+.2f}..{100*dTe.max():+.2f} %, dR {dKe.min():+.1f}..{dKe.max():+.1f} mm, "
          f"tau ratio {np.nanmin(EK[3][ok]/EG[3][ok]):.3f}..{np.nanmax(EK[3][ok]/EG[3][ok]):.3f}; "
          f"ECEI O1: dTrad {100*dTi.min():+.2f}..{100*dTi.max():+.2f} %, dR {dKR.min():+.1f}..{dKR.max():+.1f} mm, "
          f"tau ratio {np.min(K[3]/Gk[3]):.3f}..{np.max(K[3]/Gk[3]):.3f}")

    # 8. general distribution in parallel mode (pickling of a lab-frame Grad13 model)
    efk.set_distribution(Grad13Model(pi_hat=[[0.1, 0.05, 0], [0.05, -0.2, 0.1], [0, 0.1, 0.1]], q_hat=[0.02, 0.03, -0.02]))
    with contextlib.redirect_stdout(io.StringIO()):
        efk.diag.clist = efk.diag.clist[:4]
        efk.diag.vn_list, efk.diag.fn_list, efk.diag.ece_freq = efk.diag.vn_list[:4], efk.diag.fn_list[:4], efk.diag.ece_freq[:4]
    KP = run(efk, parallel=True, n_processes=2, dist='general')
    KS = run(efk, parallel=False, dist='general')
    dpk = max(np.max(np.abs(KP[k] - KS[k])) for k in (0, 1, 3, 4))
    check("dist='general': parallel == sequential (Grad13Model, lab frame)", dpk < 1e-12, f'max difference {dpk:.1e}')

    fn = f'regression_np{np.__version__}_sp{scipy.__version__}.npz'
    np.savez(fn, ecei_legacy_R=L[0], ecei_legacy_z=L[1], ecei_legacy_tau=L[3], ecei_legacy_Trad=L[4],
             ecei_general_R=G[0], ecei_general_z=G[1], ecei_general_tau=G[3], ecei_general_Trad=G[4],
             ece_legacy_R=EL[0], ece_legacy_Trad=EL[4], ece_general_R=EG[0], ece_general_Trad=EG[4],
             ecei_biM2_R=B2[0], ecei_biM2_tau=B2[3], ecei_biM2_Trad=B2[4], ece_biM2_R=EB2[0], ece_biM2_Trad=EB2[4],
             ecei_kin_R=K[0], ecei_kin_tau=K[3], ecei_kin_Trad=K[4], ece_kin_R=EK[0], ece_kin_Trad=EK[4],
             ecei_grad_R=KS[0], ecei_grad_Trad=KS[4])
    print(f'\nsaved {fn}   ({time.time()-t0:.0f} s)')
    print('OVERALL:', 'PASS' if ok_all else 'FAIL')
    sys.exit(0 if ok_all else 1)
