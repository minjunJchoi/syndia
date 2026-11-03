#!/usr/bin/env python3
"""
Bi-Maxwellian (T_par != T_perp) forward modelling example, KSTAR #37389 4150 ms (offline, no MDSplus).

    t_r = T_par / T_perp,  the Te profile given to set_profile is T_perp.

Emission and absorption are computed separately (see eceint.py, dist='bimaxwell'):
    j     = a_m(T_perp) phi_j     I_BB(T_perp)
    alpha = a_m(T_perp) phi_alpha,     phi_alpha = phi_j with weight W = 1 - beta_par cos(theta) (1 - 1/t_r)
t_r = 1 gives exactly the Maxwellian ('general') result.

usage:  python biM_example.py            (ECE radiometer, 24 channels; about 1 min)
"""
import os, sys
sys.path.insert(0, os.pardir)
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from efm import EceFwdMod

# Run only as a script: parallel runs use multiprocessing, and with the 'spawn' start method
# (Windows, macOS) the module is re-imported in each worker.
if __name__ == '__main__':
    RUN = dict(fstart=0, fend=0, Nf=1, zstart=0, zend=0, Nz=1, pstart=7.8, pend=-2, pint=0.1,
               Rinit=2.39, torbeam=0, select='max', parallel=True)

    efm = EceFwdMod()
    efm.set_profile('data/KSTAR_EFIT01_g037389.004150', 'data/ece_Te_37389_4150ms.dat',
                    'data/tci_ne_37389_4150ms.dat', bfactor=1.01)          # Te file = T_perp
    efm.set_channel(37389, ['ECE{:02d}'.format(i) for i in range(25, 49)])

    # Maxwellian reference
    R0, z0, I0, tau0, Trad0, Tabs0 = efm.run(**RUN)

    # 1) uniform anisotropy
    efm.set_anisotropy(2.0)                                                # T_par = 2 T_perp everywhere
    R1, z1, I1, tau1, Trad1, Tabs1 = efm.run(**RUN, dist='bimaxwell')

    # 2) anisotropy profile t_r(psi_n): a file with two columns [norm. psi, t_r] works the same way,
    #    efm.set_anisotropy('data/tr_profile.dat')
    psin = np.array([0.0, 0.3, 0.6, 1.0])
    tr = np.array([0.6, 0.7, 0.9, 1.0])                                    # e.g. T_perp > T_par in the core (ECRH)
    efm.set_anisotropy((psin, tr))
    R2, z2, I2, tau2, Trad2, Tabs2 = efm.run(**RUN, dist='bimaxwell')

    print('\n  f[GHz]   R[m]   tau   Trad(Maxw)[keV]   dTrad t_r=2 [%]   dTrad profile [%]   dR t_r=2 [mm]')
    for k in range(len(R0)):
        print('  {:6.1f} {:6.3f} {:6.2f} {:12.3f} {:17.2f} {:18.2f} {:16.2f}'.format(
            efm.diag.ece_freq[k], R0[k], tau0[k], Trad0[k], 100*(Trad1[k]/Trad0[k] - 1),
            100*(Trad2[k]/Trad0[k] - 1), 1e3*(R1[k] - R0[k])))
