#!/usr/bin/env python3
"""
Unit tests of the general-distribution ECE kernel (kinetic_ece.py, distributions.py);
items 9-11: bi-Maxwellian normalisation at high T, O1 line shape and bi-Maxwellian O1 strength against the kernel.
No SYNDIA data needed; about 1-2 minutes.   usage:  python kinetic_tests.py
"""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
import numpy as np
import scipy
import kinetic_ece as K
from distributions import (Juttner, Kappa, BiMaxwell, Gyrotropic, Tabulated, Arbitrary3D, Grad13)
from compat import trapezoid

ok_all = True


def check(name, cond, detail):
    global ok_all
    ok_all &= bool(cond)
    print("[{}] {}: {}".format('PASS' if cond else 'FAIL', name, detail))


print('Python {}, NumPy {}, SciPy {}'.format(sys.version.split()[0], np.__version__, scipy.__version__))
t0 = time.time()
B = 2.5
wce = K.e*B/K.me
th_list = np.radians([60, 80, 88, 90, 110])

# 1. normalisation
norms = [Juttner(0.004).norm(), Kappa(0.004, 3.0).norm(umax=60*0.063, n=1201), Kappa(0.004, 3.0, 'classic').norm(),
         BiMaxwell(0.004, 2.0).norm(), Grad13(0.004, [0, 0, 1], q_hat=[0, 0, 0.1]).norm()]
check('distributions normalised', max(abs(np.array(norms) - 1)) < 2e-4,
      'Juttner, kappa(rel), kappa(classic), bi-Maxwell, Grad13: ' + ', '.join('{:.6f}'.format(v) for v in norms))

# 2. Kirchhoff for Juttner
worst = 0.0
for TkeV in (0.3, 3.0, 10.0):
    d = Juttner(TkeV*1e3*K.e/K.mc2)
    for md, m in (('X', 2), ('O', 1), ('X', 3)):
        for th in th_list:
            for ne in (1e17, 3e19):
                for x in (0.96, 0.995, 1.003):
                    w = x*m*wce
                    j, a = K.emission_absorption(d, w, B, ne, th, md)
                    if j > 0:
                        worst = max(worst, abs(j/(a*K.I_BB(w, TkeV*1e3*K.e)) - 1))
check('Kirchhoff j = alpha I_BB (Juttner)', worst < 1e-12, 'max |j/(alpha I_BB) - 1| = {:.1e}'.format(worst))

# 3. absolute O1 absorption vs the hot-plasma result Int alpha domega = (pi/2) wp^2 Theta / c
ne, Th = 1e17, 100*K.e/K.mc2
ws = wce*np.linspace(0.985, 1.01, 1501)
a = np.array([K.emission_absorption(Juttner(Th), w, B, ne, np.pi/2, 'O', harmonics=[1])[1] for w in ws])
r = trapezoid(a, ws)/(0.5*np.pi*ne*K.e**2/(K.eps*K.me)/K.c*Th)
check('O1 absolute absorption (100 eV, 90 deg)', abs(r - 1) < 0.01, 'ratio to (pi/2) wp^2 Theta/c = {:.4f}'.format(r))

# 4. polarisation: X-mode at small theta is the R wave (rotates with the electrons)
N, ev, et2 = K.mode_props(0.9*wce, B, 1e18, 1e-3, 'X')
check('X-mode -> R wave at theta -> 0', abs(ev[1]/ev[0] - 1j) < 1e-3, 'e_y/e_x = {:.4f}'.format(ev[1]/ev[0]))

# 5. bi-Maxwellian: t_r = 1 is exactly Juttner; at 90 deg the local T_rad = j/(alpha I_BB) is T_perp to O(Theta)
#    (the relativistic bi-Maxwellian is exp(-(sqrt(1+u_perp^2+u_par^2/t_r)-1)/Th_perp))
d1, dJ0 = BiMaxwell(0.004, 1.0), Juttner(0.004)
ug = np.linspace(0, 0.6, 7)
dev_id = max(abs(d1.F(ug, 0.3*ug)/dJ0.F(ug, 0.3*ug) - 1))
worst = 0.0
for tr in (0.4, 2.5):
    d = BiMaxwell(0.004, tr)
    for x in (0.98, 0.995):
        j, a = K.emission_absorption(d, x*2*wce, B, 2e19, np.pi/2, 'X')
        worst = max(worst, abs(j/(a*K.I_BB(x*2*wce, 0.004*K.mc2)) - 1))
check('bi-Maxwellian: t_r = 1 is Juttner; 90 deg T_rad = T_perp to O(Theta)', dev_id < 1e-12 and worst < 2e-2,
      't_r = 1 deviation {:.1e}; T_rad/T_perp - 1 up to {:.1e} at Theta = 0.004'.format(dev_id, worst))

# 6. kappa -> Juttner for large kappa
th = np.radians(80)
dJ = Juttner(0.004)
dev = []
for kap in (20.0, 200.0):
    dk = Kappa(0.004, kap)
    rr = [K.emission_absorption(dk, x*2*wce, B, 2e19, th, 'X')[0]/K.emission_absorption(dJ, x*2*wce, B, 2e19, th, 'X')[0] - 1
          for x in (0.99, 1.0)]
    dev.append(max(abs(np.array(rr))))
check('kappa -> Juttner', dev[1] < 0.2*dev[0] and dev[1] < 0.02, 'max |j_kappa/j_J - 1|: kappa=20 {:.3f}, kappa=200 {:.4f}'.format(*dev))

# 7. non-gyrotropic parts drop out; lab-frame stress = its field-aligned projection
Th = 0.004
pi_perp_offdiag = [[0.3, 0.2, 0.15], [0.2, -0.3, 0.1], [0.15, 0.1, 0.0]]     # traceless, b = z: no b.pi.b part
g1 = Grad13(Th, [0, 0, 1], pi_hat=pi_perp_offdiag, q_hat=[0.2, -0.1, 0.0])
dev1 = 0.0
for x, md, m in ((0.995, 'X', 2), (0.99, 'O', 1)):
    jg, ag = K.emission_absorption(g1, x*m*wce, B, 2e19, th, md)
    jj, aj = K.emission_absorption(Juttner(Th), x*m*wce, B, 2e19, th, md)
    dev1 = max(dev1, abs(jg/jj - 1), abs(ag/aj - 1))
check('Grad13: perpendicular/off-diagonal stress and q_perp give the Maxwellian', dev1 < 1e-7, 'max deviation {:.1e} (numerical gyro-average and derivatives)'.format(dev1))

rng = np.random.RandomState(1)
A = rng.normal(size=(3, 3)); A = 0.15*(A + A.T); A -= np.trace(A)/3*np.eye(3)
q = np.array([0.05, -0.08, 0.06])
b = np.array([0.3, 0.9, -0.2]); b /= np.linalg.norm(b)
A = A + 0.3*(np.outer(b, b) - np.eye(3)/3)
g_lab = Grad13(Th, b, pi_hat=A, q_hat=q)
pp, qp = g_lab.projected()
g_fld = Grad13(Th, [0, 0, 1], pi_hat=np.diag([-pp/2, -pp/2, pp]), q_hat=[0, 0, qp])
dev2 = 0.0
for x, md, m in ((0.995, 'X', 2), (1.004, 'X', 2), (0.99, 'O', 1)):
    j1, a1 = K.emission_absorption(g_lab, x*m*wce, B, 2e19, th, md)
    j2, a2 = K.emission_absorption(g_fld, x*m*wce, B, 2e19, th, md)
    dev2 = max(dev2, abs(j1/j2 - 1), abs(a1/a2 - 1))
check('Grad13 in a tilted lab frame = its field-aligned projection', dev2 < 1e-9,
      'pi_par = {:.3f}, q_par = {:.3f}, max deviation {:.1e}'.format(pp, qp, dev2))

# 8. generic Arbitrary3D (non-separable f) and user gyrotropic function / table reproduce Juttner
def f_xy(u1, u2, u3, Theta=0.004):
    g = np.sqrt(1 + u1*u1 + u2*u2 + u3*u3)
    return np.exp(-(g - 1)/Theta)*(1 + 3.0*u1*u2/Theta)


def F_J(uperp, upar):
    return np.exp(-(np.sqrt(1 + uperp**2 + upar**2) - 1)/0.004)


a3 = Arbitrary3D(f_xy, [0, 0, 1], scale=np.sqrt(0.004))
gy = Gyrotropic(F_J, scale=np.sqrt(0.004))
up = np.linspace(0, 1.0, 801); ul = np.linspace(-1.0, 1.0, 1601)
tb = Tabulated(up, ul, F_J(up[:, None], ul[None, :]))
dev3 = {}
for lab, d in (('Arbitrary3D(u_x u_y)', a3), ('Gyrotropic(user F)', gy), ('Tabulated', tb)):
    dv = 0.0
    for x, md, m in ((0.995, 'X', 2), (0.99, 'O', 1)):
        j1, a1 = K.emission_absorption(d, x*m*wce, B, 2e19, th, md)
        j2, a2 = K.emission_absorption(Juttner(0.004), x*m*wce, B, 2e19, th, md)
        dv = max(dv, abs(j1/j2 - 1), abs(a1/a2 - 1))
    dev3[lab] = dv
check('user-supplied distributions reproduce Juttner', dev3['Arbitrary3D(u_x u_y)'] < 1e-6 and dev3['Gyrotropic(user F)'] < 1e-6
      and dev3['Tabulated'] < 2e-3, ', '.join('{} {:.1e}'.format(k, v) for k, v in dev3.items()))

# 9. bi-Maxwellian normalisation at high temperature and t_r > 1
import warnings
worst, nwarn = 0.0, 0
for TkeV in (2.0, 10.0, 20.0):
    for tr in (0.5, 1.5, 2.0, 3.0):
        d = BiMaxwell(TkeV*1e3*K.e/K.mc2, tr)
        with warnings.catch_warnings(record=True) as wl:
            warnings.simplefilter('always')
            nv = d.norm(umax=30*d.scale, n=801)
        nwarn += len(wl)
        worst = max(worst, abs(nv - 1))
check('bi-Maxwellian normalised for T_perp = 2-20 keV, t_r = 0.5-3', worst < 1e-3 and nwarn == 0,
      'max |norm - 1| = {:.1e}, tail warnings {}'.format(worst, nwarn))

# 10. SYNDIA O1 line shape (eceint, shape_model='general') vs the kinetic kernel (m = 1, N = 1, low density)
import eceint as E
def _stats(x, y):
    y = np.clip(np.nan_to_num(y), 0, None); n = trapezoid(y, x); mu_ = trapezoid(x*y, x)/n
    return mu_ - 1, np.sqrt(trapezoid((x - mu_)**2*y, x)/n)
xs = np.linspace(0.97, 1.02, 801)
dev_sh, dev_w = 0.0, 0.0
for thd in (90, 85, 80):
    thr = np.radians(thd); T = 1e3*K.e; zeta = K.mc2/T/2
    ak = [K.emission_absorption(Juttner(T/K.mc2), x*wce, B, 1e17, thr, 'O', harmonics=[1], N_res='vacuum')[1] for x in xs]
    ph = [E.shape_fn(1, zeta, wce, x*wce, thr, 'general') for x in xs]
    (mk, sk), (mp, sp) = _stats(xs, ak), _stats(xs, ph)
    dev_sh, dev_w = max(dev_sh, abs(mp/mk - 1)), max(dev_w, abs(sp/sk - 1))
check("O1 line shape ('general') = kinetic kernel, 1 keV, 80-90 deg", dev_sh < 0.05 and dev_w < 0.05,
      'mean shift within {:.1%}, rms width within {:.1%}'.format(dev_sh, dev_w))

# 11. SYNDIA bi-Maxwellian O1 (eceint.shape_pair_biM) vs the kinetic kernel: t_r dependence of the integrated absorption
th85 = np.radians(85); T = 1e3*K.e; zeta = K.mc2/T/2
xs = np.linspace(0.97, 1.02, 601)
def _int_alpha_kin(tr):
    return trapezoid([K.emission_absorption(BiMaxwell(T/K.mc2, tr), x*wce, B, 1e17, th85, 'O', harmonics=[1], N_res='vacuum')[1]
                      for x in xs], xs)
def _int_phi(tr):
    return trapezoid([E.shape_pair_biM(1, zeta, tr, wce, x*wce, th85)[1] for x in xs], xs)
rk = [_int_alpha_kin(tr)/_int_alpha_kin(1.0) for tr in (0.5, 2.0)]
rs = [_int_phi(tr)/_int_phi(1.0) for tr in (0.5, 2.0)]
dev11 = max(abs(rs[i]/rk[i] - 1) for i in range(2))
check('bi-Maxwellian O1 strength vs t_r: SYNDIA = kinetic kernel (1 keV, 85 deg)', dev11 < 0.03,
      'Int alpha(t_r)/Int alpha(1): kinetic {:.3f}, {:.3f}; SYNDIA {:.3f}, {:.3f} (t_r = 0.5, 2)'.format(*rk, *rs))

print('\n({:.0f} s)'.format(time.time() - t0))
print('OVERALL:', 'PASS' if ok_all else 'FAIL')
sys.exit(0 if ok_all else 1)
