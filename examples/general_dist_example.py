#!/usr/bin/env python3
"""
General electron distributions in SYNDIA (dist='general'), KSTAR #37389 4150 ms example (offline).

The whole chain is:   distribution model -> local f at every path point -> gyro-average F(u_perp, u_par)
  -> emission j and absorption alpha (kinetic_ece: relativistic, exact Bessel, cold-plasma polarisation)
  -> SYNDIA ray tracing (cold-plasma refraction) and radiative transfer -> T_rad, emission position, tau.

Four ways to give a distribution:
  1. a built-in model: MaxwellModel, KappaModel, BiMaxwellModel, Grad13Model
  2. a function of the local plasma that returns a Distribution (FunctionModel) -- shown with a
     NON-SEPARABLE, NON-GYROTROPIC f(u_R, u_phi, u_Z) in the lab frame (an off-diagonal 'viscous stress' term)
  3. a gyrotropic F(u_perp, u_par) (distributions.Gyrotropic) -- shown with a Maxwellian + a parallel beam
  4. a table F(u_perp, u_par), e.g. from a Fokker-Planck code (distributions.Tabulated)
Functions used by a model must be defined at module level (not lambdas) so that parallel runs can pickle them.

usage:  python general_dist_example.py        (8 ECE channels; a few minutes)
"""
import os, sys
sys.path.insert(0, os.pardir)
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from efm import EceFwdMod
from distributions import (MaxwellModel, KappaModel, Grad13Model, FunctionModel, Arbitrary3D, Gyrotropic,
                           Juttner, MC2_J)

# ---- 2. a non-separable lab-frame distribution:  f = f_J(u) [1 + a (u_R u_phi + u_phi u_Z)/Theta]
#         (u1, u2, u3) = (u_R, u_phi, u_Z); only the part of the stress that projects on b survives the gyro-average.
#         The correction integrates to zero, so f stays normalised by the Juttner constant Z.
def f_stress(u1, u2, u3, Theta=0.004, Z=1.0, a=0.1):
    g = np.sqrt(1 + u1*u1 + u2*u2 + u3*u3)
    return np.exp(-(g - 1)/Theta)*(1 + a*(u1*u2 + u2*u3)/Theta)/Z


def stress_model(R, z, Te, ne, Bvec, psin, a=0.1):
    Th = Te/MC2_J
    b = np.asarray(Bvec)/np.linalg.norm(Bvec)                  # field direction in the (R, phi, Z) frame
    return Arbitrary3D(f_stress, b, scale=np.sqrt(Th), fargs={'Theta': Th, 'Z': Juttner(Th).Z, 'a': a},
                       nphi=16, normalize=False)


# ---- 3. a gyrotropic distribution: Maxwellian core + 1 % parallel beam at u_par = 0.3 (width 0.05)
class BeamF(object):
    """Picklable callable F(uperp, upar), normalised analytically (core 99 %, beam 1 %)."""
    def __init__(self, Theta, frac=0.01, u0=0.3, w=0.05):
        self.Theta, self.frac, self.u0, self.w = Theta, frac, u0, w
        self.ZJ = Juttner(Theta).Z

    def __call__(self, uperp, upar):
        g = np.sqrt(1 + uperp**2 + upar**2)
        core = np.exp(-(g - 1)/self.Theta)/self.ZJ
        beam = np.exp(-(uperp**2 + (upar - self.u0)**2)/(2*self.w**2))/((2*np.pi)**1.5*self.w**3)
        return (1 - self.frac)*core + self.frac*beam


def beam_model(R, z, Te, ne, Bvec, psin):
    Th = Te/MC2_J
    return Gyrotropic(BeamF(Th), scale=np.sqrt(Th), normalize=False)


# Run only as a script: parallel runs use multiprocessing, and with the 'spawn' start method
# (Windows, macOS) the module is re-imported in each worker.
if __name__ == '__main__':
    RUN = dict(fstart=0, fend=0, Nf=1, zstart=0, zend=0, Nz=1, pstart=7.8, pend=-2, pint=0.1,
               Rinit=2.39, torbeam=0, select='max', parallel=True)

    efm = EceFwdMod()
    efm.set_profile('data/KSTAR_EFIT01_g037389.004150', 'data/ece_Te_37389_4150ms.dat',
                    'data/tci_ne_37389_4150ms.dat', bfactor=1.01)
    efm.set_channel(37389, ['ECE{:02d}'.format(i) for i in range(27, 49, 3)])

    cases = [('Maxwell (kinetic)', MaxwellModel()),
             ('kappa = 5', KappaModel(5.0)),
             ('Grad13 lab-frame stress', Grad13Model(pi_hat=[[0, 0.1, 0], [0.1, 0, 0.1], [0, 0.1, 0]])),
             ('non-separable f (FunctionModel)', FunctionModel(stress_model, a=0.1)),
             ('Maxwell + 1% beam', FunctionModel(beam_model))]
    res = {}
    for lab, model in cases:
        efm.set_distribution(model)
        res[lab] = efm.run(**RUN, dist='general')

    R0, T0 = np.asarray(res['Maxwell (kinetic)'][0]), np.asarray(res['Maxwell (kinetic)'][4])
    print('\nf [GHz]  R [m]  Trad(Maxwell) [keV] | T_rad change [%]: ' + ' | '.join(l for l, _ in cases[1:]))
    for k in range(len(R0)):
        print('{:6.1f} {:6.3f} {:10.3f}  '.format(efm.diag.ece_freq[k], R0[k], T0[k]) +
              ''.join('{:14.2f}'.format(100*(res[l][4][k]/T0[k] - 1)) for l, _ in cases[1:]))
