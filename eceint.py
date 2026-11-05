"""
ECE intensity calculation
Repo : https://github.com/minjunJchoi/syndia
Author : Minjun J. Choi (mjchoi@nfri.re.kr)
Collaborators : Jieun Lee, Yoonbum Nam
References : Hutchinson, Bornatici, Rathgeber
"""
import numpy as np
import math
import cmath
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from compat import trapezoid, cumulative_trapezoid, simpson, scalar

# all mks units except Temperature in [J]
# s : beam path [m] from hfs to lfs
# Rp : R coordinates on beam path [m]
# zp : z coordinate on beam path [m]
# th : angle between field and beam path [rad]
# omega : emission frequency [rad/s]
# m : harmonic number
# B : B field function in R-z coordinate [T]
# F_Te : 2d TriScatteredInterp Te function in R-z coordinates [J]
# F_ne : 2d TriScatteredInterp ne function in R-z coordinates [m^-3]
# @(x,y) : (s) coordinate
# select : 'mean' or 'max'

# it will calculate emission profile of frequency omega along s
# Rm : maximum emission position [m]
# zm : maximum emission position [m]
# ece_int : ece intensity [W/m^2 Hz]

e = 1.602*1e-19
me = 9.109*1e-31
eps = 8.854*1e-12
c = 299792458
mc2 = me*c**2

# ---------------------------------------------------------------------------------
# ECE line shape.  shape_model:
#   'general' (default): coupling-weighted integrand + peak-adaptive quadrature, Int phi dw = 2 pi.
#     X-mode, m >= 2 (beta_perp^(2m) weight, Rathgeber 2013 / Leem 2019 written in beta_par):
#       phi_m(w) = (4 sqrt(pi)/m!) zeta^((2m+3)/2) w/(m wce)^2
#                  * Int dba (1 - ba cos th) bperp2^m exp(-zeta (1 - (1 - ba cos th)^2 mu)),
#       bperp2 = 1 - ba^2 - (1 - ba cos th)^2 mu  (= beta_perp^2 on the resonance).  Same as the original for m = 2.
#     O-mode, m = 1 (SYNDIA's hn = 1): the electric field is along B, so the coupling is
#       beta_par^2 J_1^2 ~ beta_par^2 beta_perp^2 (Bekefi):
#       phi_1(w) = 8 sqrt(pi) zeta^(7/2) w/wce^2 * Int dba (1 - ba cos th) ba^2 bperp2 exp(-zeta (...)),
#       (<beta_par^2 beta_perp^2> = 1/(2 zeta^2) fixes the prefactor).  At theta = 90 deg the resonance depends only
#       on |beta|, the angle averages of beta_perp^4 and beta_par^2 beta_perp^2 are both ~ beta^4, and this is
#       identical to the original (m = 2 integrand); at oblique angles only this weight reproduces the width of a
#       fully relativistic kinetic calculation.
#   'legacy' : the original expression, unchanged (m = 2 integrand for every m, scipy quad defaults),
#       kept to reproduce earlier results bit-for-bit.
# The scipy quad defaults (an empirical choice that is fine in SYNDIA's usual range) can miss the narrow peak for
# T_e <~ 300 eV at |90 deg - theta| >~ 10 deg (line-shape normalisation off by up to +-60 %); 'general' finds the
# peak first and integrates it densely, so the line shape stays normalised there as well.
# ---------------------------------------------------------------------------------
SHAPE_MODELS = ('general', 'legacy')


def _robust_int(F, a, b, n=2001):
    """Peak-adaptive Simpson: find the support (>= 1e-13 of max) on a coarse grid, integrate densely there."""
    g = np.linspace(a, b, n)
    v = F(g)
    if not np.any(v > 0):
        return 0.0
    idx = np.where(v >= 1e-13*v.max())[0]
    lo, hi = g[max(idx[0]-1, 0)], g[min(idx[-1]+1, n-1)]
    g2 = np.linspace(lo, hi, n)
    return simpson(F(g2), x=g2)


def shape_fn(m, zeta_s, wce_s, w, th, shape_model='general'):
    """Normalized line shape phi_m(w) [s] at one point of the path (Int phi dw = 2 pi)."""
    mu_s = (w/(m*wce_s))**2
    ct = np.cos(th)
    ba1 = (mu_s*ct - np.lib.scimath.sqrt(1 - mu_s*np.sin(th)**2))/(1 + mu_s*ct**2)
    ba2 = (mu_s*ct + np.lib.scimath.sqrt(1 - mu_s*np.sin(th)**2))/(1 + mu_s*ct**2)
    if shape_model == 'legacy':
        return 2.0*np.pi*zeta_s**(3.5)*w/(np.sqrt(np.pi)*(m*wce_s)**2.0) \
            *integrate.quad(lambda ba: (1 - ba**2.0 - (1 - ba*ct)**2.0*mu_s)**2.0 \
            *(1 - ba*ct)*np.exp(-zeta_s*(1 - (1 - ba*ct)**2.0*mu_s)), ba1, ba2)[0]
    if shape_model != 'general':
        raise ValueError("shape_model must be one of {}".format(SHAPE_MODELS))
    if mu_s*np.sin(th)**2 >= 1:
        return 0.0
    ba1, ba2 = float(np.real(ba1)), float(np.real(ba2))
    if m == 1:   # O-mode fundamental: beta_par^2 beta_perp^2 coupling
        F = lambda ba: ba**2*np.clip(1 - ba**2 - (1 - ba*ct)**2*mu_s, 0, None) * (1 - ba*ct) \
            * np.exp(-zeta_s*(1 - (1 - ba*ct)**2*mu_s))
        pref = 8.0*np.sqrt(np.pi)*zeta_s**3.5*w/wce_s**2
    else:        # X-mode: beta_perp^(2m) coupling
        F = lambda ba: np.clip(1 - ba**2 - (1 - ba*ct)**2*mu_s, 0, None)**m * (1 - ba*ct) \
            * np.exp(-zeta_s*(1 - (1 - ba*ct)**2*mu_s))
        pref = 4.0*np.sqrt(np.pi)/math.factorial(m)*zeta_s**((2*m+3)/2.0)*w/(m*wce_s)**2
    return pref*_robust_int(F, ba1, ba2)


# ---------------------------------------------------------------------------------
# Electron distribution.  dist:
#   'maxwell'   (default): isotropic Maxwellian, j = alpha * I_BB(T_e)  (Kirchhoff), as before.
#   'bimaxwell' : bi-Maxwellian with t_r = T_par/T_perp (T_e profile = T_perp).  j and alpha are
#       computed separately:
#         j     = a_m(T_perp) * phi_j     * I_BB(T_perp)
#         alpha = a_m(T_perp) * phi_alpha
#       phi_j     = (4 sqrt(pi)/m!) zeta_perp^(m+1) zeta_par^(1/2) w/(m wce)^2
#                   * Int dba (1 - ba cos th) bperp2^m exp(-zeta_perp (1 - (1 - ba cos th)^2 mu) + (zeta_perp - zeta_par) ba^2)
#       phi_alpha = same integral with the weight W = 1 - ba cos th (1 - 1/t_r)
#       (detailed balance with the photon momentum == Stix resonant operator; N = 1 as in the resonance).
#       For m = 1 (O-mode) the integrand carries the beta_par^2 beta_perp^2 coupling and the prefactor is normalised
#       with zeta_perp, so Int phi dw = 2 pi t_r: the T_par dependence of the O1 strength comes out of the integral
#       t_r = 1 reproduces shape_model='general' (to rounding).  Requires shape_model='general'.
#       W < 0 (possible only for strong T_perp > T_par) means negative absorption (maser); a warning is issued.
# ---------------------------------------------------------------------------------
DIST_MODELS = ('maxwell', 'bimaxwell', 'general')
AMP_MODELS = ('corrected', 'legacy')
# amp_model: Hutchinson absorption prefactors a_m (used with dist='maxwell' and 'bimaxwell').
#   'corrected' (default):
#     (1) cold-plasma index by Appleton-Hartree, Delta = sqrt((Y^2 sin^2/2)^2 + (1-X)^2 Y^2 cos^2), Y = 1/m.
#         The original had sqrt((sin^2/2m^2)^2 - (1-X) cos^2/m^2) (sign error, (1-X) not squared): identical at
#         90 deg, but N^2 wrong/complex at oblique angles (X2 a_2 too small by 5-7 % at 85 deg and 20-40 % at 80 deg
#         for X = 0.3-0.01).
#     (2) O-mode fundamental: a_1 divided by 2 pi.  With SYNDIA's line shape (Int phi domega = 2 pi) the original
#         gave Int alpha domega = pi^2 wp^2 Theta / c, whereas the hot-plasma result (eps_zz, n = 1) and the
#         exact kernel (kinetic_ece) and Bekefi's single-particle emission with Kirchhoff's law give (pi/2) wp^2 Theta / c.  The X-mode prefactor was already consistent.
#   'legacy' : the original expressions (reproduces earlier results).
# dist='general': any electron distribution through kinetic_ece.emission_absorption
#   (fully relativistic, exact Bessel, cold-plasma polarisation, j and alpha from one kernel).  Needs
#   dists = list of distributions.Distribution, one per path point (built by efm from a DistributionModel).


def _support_grid(F, a, b, n=2001):
    """Same support search as _robust_int, returning the dense grid (None if F == 0 everywhere)."""
    g = np.linspace(a, b, n)
    v = F(g)
    if not np.any(v > 0):
        return None
    idx = np.where(v >= 1e-13*v.max())[0]
    lo, hi = g[max(idx[0]-1, 0)], g[min(idx[-1]+1, n-1)]
    return np.linspace(lo, hi, n)


def shape_pair_biM(m, zeta_perp, t_r, wce_s, w, th):
    """(phi_j, phi_alpha) [s] for a bi-Maxwellian; zeta_perp = m c^2/(2 T_perp), t_r = T_par/T_perp."""
    zeta_par = zeta_perp/t_r
    mu_s = (w/(m*wce_s))**2
    ct, st = np.cos(th), np.sin(th)
    if mu_s*st**2 >= 1:
        return 0.0, 0.0
    sq = np.sqrt(1 - mu_s*st**2)
    ba1 = (mu_s*ct - sq)/(1 + mu_s*ct**2)
    ba2 = (mu_s*ct + sq)/(1 + mu_s*ct**2)
    wpar = (lambda ba: ba**2) if m == 1 else (lambda ba: 1.0)      # O1: beta_par^2 beta_perp^2 coupling
    F = lambda ba: wpar(ba)*np.clip(1 - ba**2 - (1 - ba*ct)**2*mu_s, 0, None)**m * (1 - ba*ct) \
        * np.exp(-zeta_perp*(1 - (1 - ba*ct)**2*mu_s) + (zeta_perp - zeta_par)*ba**2)
    g = _support_grid(F, ba1, ba2)
    if g is None:
        return 0.0, 0.0
    Fg = F(g)
    pref = 4.0*np.sqrt(np.pi)/math.factorial(m)*zeta_perp**(m + 1)*np.sqrt(zeta_par)*w/(m*wce_s)**2
    if m == 1:
        pref *= 2.0*zeta_perp   # normalised with <beta_par^2> at T_perp -> Int phi dw = 2 pi t_r
    return pref*simpson(Fg, x=g), pref*simpson(Fg*(1 - g*ct*(1 - 1/t_r)), x=g)


def ece_intensity(s, Rp, zp, th, omega, m, F_B, F_Te, F_ne, select='mean', use_cache=True, shape_model='general',
                  dist='maxwell', F_tr=None, amp_model='corrected', dists=None, kin_opts=None): # [m], [m], [rad], [rad/s], harmonic number
    """dist='bimaxwell' needs F_tr(s) = T_par/T_perp along the path (F_Te is then T_perp)."""
    if dist not in DIST_MODELS:
        raise ValueError("dist must be one of {}".format(DIST_MODELS))
    if amp_model not in AMP_MODELS:
        raise ValueError("amp_model must be one of {}".format(AMP_MODELS))
    if dist == 'general' and (dists is None or len(dists) != len(s)):
        raise ValueError("dist='general' requires dists (one Distribution per path point)")
    if dist == 'bimaxwell':
        if shape_model != 'general':
            raise ValueError("dist='bimaxwell' requires shape_model='general'")
        if F_tr is None:
            raise ValueError("dist='bimaxwell' requires F_tr (T_par/T_perp along the path)")
    if use_cache:
        # Small cache for evaluations to reduce expensive F_B, F_ne calls
        _cache = {}
        def _q(v):
            return round(scalar(v), 5)

    # characteristic frequencies
    if use_cache:
        def wce(s):
            key = ('wce', _q(s))
            val = _cache.get(key)
            if val is None:
                val = e*F_B(s)/me # [rad/s]
                _cache[key] = val
            return val   
        def wpe(s):
            key = ('wpe', _q(s))
            val = _cache.get(key)
            if val is None:
                val = np.sqrt(F_ne(s)*e**2/(eps*me)) # [rad/s]
                _cache[key] = val
            return val
    else:
        wce = lambda s: e*F_B(s)/me # [rad/s]
        wpe = lambda s: np.sqrt(F_ne(s)*e**2/(eps*me)) # [rad/s]

    # characteristic functions
    if use_cache:
        def zeta(s):
            key = ('zeta', _q(s))
            val = _cache.get(key)
            if val is None:
                val = mc2/F_Te(s)/2
                _cache[key] = val
            return val
        def mu(s, w):
            key = ('mu', _q(s), _q(w))
            val = _cache.get(key)
            if val is None:
                val = (w/(m*wce(s)))**2
                _cache[key] = val
            return val
    else:
        zeta = lambda s: mc2/F_Te(s)/2
        mu = lambda s,w: (w/(m*wce(s)))**2

    # refractive index
    X = lambda s: ( wpe(s)/(m*wce(s)) )**2
    if m == 1: # O-mode
    #     Nperp2C = @(x,y,w) 1 - (wpe(x,y)./w).^2; % B(3.1.22)
        if amp_model == 'legacy':
            N1OCsq = lambda s,th: 1 - X(s)*(1-X(s))/(1 - X(s) - np.sin(th)**2/(2*m**2) + cmath.sqrt((np.sin(th)**2/(2*m**2))**2 - (1 - X(s))*np.cos(th)**2/m**2))
        else:
            N1OCsq = lambda s,th: 1 - X(s)*(1-X(s))/(1 - X(s) - np.sin(th)**2/(2*m**2) + cmath.sqrt((np.sin(th)**2/(2*m**2))**2 + (1 - X(s))**2*np.cos(th)**2/m**2))
        N1OCre = lambda s,th: (np.lib.scimath.sqrt(N1OCsq(s,th))).real
    elif m == 2: # X-mode
    #     Nperp2C = @(x,y,w) 1 - (wpe(x,y)./w).^2.*(w.^2 - wpe(x,y).^2)./(w.^2 - wce(x,y).^2 - wpe(x,y).^2); % B(3.1.12)
        if amp_model == 'legacy':
            N2XCsq = lambda s,th: 1 - X(s)*(1-X(s))/(1 - X(s) - np.sin(th)**2/(2*m**2) - cmath.sqrt((np.sin(th)**2/(2*m**2))**2 - (1 - X(s))*np.cos(th)**2/m**2)) # H(5.2.48)
        else:
            N2XCsq = lambda s,th: 1 - X(s)*(1-X(s))/(1 - X(s) - np.sin(th)**2/(2*m**2) - cmath.sqrt((np.sin(th)**2/(2*m**2))**2 + (1 - X(s))**2*np.cos(th)**2/m**2)) # Appleton-Hartree, Y = 1/m
        N2XCre = lambda s,th: (np.lib.scimath.sqrt(N2XCsq(s,th))).real

    # absorption coefficient
    if m == 1: # O-mode
        amO = lambda s,th: np.pi/(2*c)*wpe(s)**2*N1OCre(s,th)*(1 + 2*np.cos(th)**2)**2*np.sin(th)**4/(1 + np.cos(th)**2)**3*F_Te(s)/mc2 # H(5.2.52)
        o1fac = 1.0 if amp_model == 'legacy' else 1.0/(2*np.pi)  # see AMP_MODELS
        amO_ = amO
        amO = lambda s,th: o1fac*amO_(s,th)

        amRz = lambda s,th: amO(s,th)
    elif m == 2: # X-mode
        am = lambda s: e**2*F_ne(s)/(4*c*me*eps) * m**(2*m-1)/math.factorial(m-1) * (F_Te(s)/(2*mc2))**(m-1) # H(5.2.39)

        a2sq = lambda s,th: (1 + (1 - X(s))*N2XCre(s,th)**2*np.cos(th)**2/(1 - X(s) - N2XCre(s,th)**2*np.sin(th)**2)**2 \
        * m**2*(1 - (m**2 -1)/m**2/X(s)*(1 - N2XCre(s,th)**2))**2)**2*np.sin(th)**2

        b2sq = lambda s,th: (1 + (1 - X(s))/(1 - X(s) - N2XCre(s,th)**2*np.sin(th)**2) \
        * m**2*(1 - (m**2 -1)/m**2/X(s)*(1 - N2XCre(s,th)**2))**2)**2*np.cos(th)**2

        eta2X = lambda s,th: N2XCre(s,th)**(2*m-3)*(m - 1)**2*(1 - (m+1)/m/X(s)*(1 - N2XCre(s,th)**2))**2 \
        / ((1 + np.cos(th)**2)*(a2sq(s,th) + b2sq(s,th))**(1.0/2.0)) # H(5.2.47)

        amX = lambda s,th: am(s) * eta2X(s,th) # H(5.2.54)

        amRz = lambda s,th: amX(s,th)

    # shape Maxwellian (relativistic + Doppler)
    ba1 = lambda s,w,th: (mu(s,w)*np.cos(th) - np.lib.scimath.sqrt(1 - mu(s,w)*np.sin(th)**2))/(1 + mu(s,w)*np.cos(th)**2)
    ba2 = lambda s,w,th: (mu(s,w)*np.cos(th) + np.lib.scimath.sqrt(1 - mu(s,w)*np.sin(th)**2))/(1 + mu(s,w)*np.cos(th)**2)

    # line shape (see shape_fn above; shape_model='legacy' reproduces the original code)
    shape = lambda s,w,th: shape_fn(m, zeta(s), wce(s), w, th, shape_model)

    # perfect blackbody intensity
    Ibb = lambda s,w: (w/(2*np.pi*c))**2*F_Te(s) # H(5.2.37)

    ## start calculation
    # absorption coefficent along beam path s (from hfs to lfs)
    ams = np.zeros(s.size)
    if dist in ('bimaxwell', 'general'):
        ems = np.zeros(s.size)  # bimaxwell: j/I_BB(T_perp); general: j itself
    if dist == 'general':
        import kinetic_ece
        ko = dict(N_res='cold', n_harm=3)
        ko.update(kin_opts or {})
        mode = 'O' if m == 1 else 'X'
        for i in range(s.size):
            Bi, nei = scalar(F_B(s[i])), scalar(F_ne(s[i]))
            ji, ai = kinetic_ece.emission_absorption(dists[i], omega, Bi, nei, th[i], mode, **ko)
            ams[i] = ai
            ems[i] = ji  # the emission coefficient itself (not j/I_BB: I_BB -> 0 where the T_e spline <= 0)
        nneg = sum(getattr(d, 'n_negative', 0) for d in dists)
        if nneg:   # reported here so that it is also seen in parallel runs (each worker has its own copies)
            import warnings
            warnings.warn('dist=general: negative values of the gyro-averaged distribution (far tail of a Grad/Arbitrary3D '
                          'expansion) were clipped to 0; the count per distribution is in Distribution.n_negative',
                          RuntimeWarning)
    for i in (range(s.size) if dist != 'general' else []):
        if mu(s[i], omega)*np.sin(th[i])**2 < 1: # the integral blows up with the imaginary beta
            if dist == 'bimaxwell':
                t_r = scalar(F_tr(s[i]))
                a_m = (amRz(s[i], th[i])).real  # a_m(T_perp); for O1 the T_par dependence is in phi (see shape_pair_biM)
                phi_j, phi_a = shape_pair_biM(m, zeta(s[i]), t_r, wce(s[i]), omega, th[i])
                ems[i] = a_m*phi_j
                ams[i] = a_m*phi_a
            else:
                ams[i] = (amRz(s[i], th[i])).real*(shape(s[i], omega, th[i])).real
        else:
            ams[i] = 0
            if dist == 'bimaxwell':
                ems[i] = 0
    if dist in ('bimaxwell', 'general') and np.any(ams < 0):
        import warnings
        warnings.warn('{}: negative absorption (population inversion / maser regime) at {} path points'.format(dist, 
            int(np.sum(ams < 0))), RuntimeWarning)

    # # define path from the inside
    # s = np.zeros(Rp.size)
    ds = np.zeros(Rp.size)
    for i in range(1,Rp.size):
        ds[i] = np.sqrt((Rp[i] - Rp[i-1])**2 + (zp[i] - zp[i-1])**2)
    #     s[i] = s[i-1] + ds[i]

    # calculate differential optical depth dtau
    dtau = trapezoid(ams,x=s) - np.append(0,cumulative_trapezoid(ams,x=s)) # from zero; a single point matters

    # emissitivty after reabsorption
    jms = np.zeros(s.size)
    for i in range(s.size):
        if dist == 'general':
            jem = ems[i]  # kinetic kernel returns j directly
        else:
            emi = ems[i] if dist == 'bimaxwell' else ams[i] # Kirchhoff (j = alpha I_BB) only for a Maxwellian
            jem = emi*Ibb(s[i],omega)
        jms[i] = jem*np.exp(-dtau[i]) # emissivity after reabsorption. B(2.2.13), B(2.2.15)
   
    # max or mean emissivity position
    if select == 'max': 
        # maximum emissivity position
        midx = np.where(jms == jms.max())
        Rm = np.mean(Rp[midx]) # for multiple maximum cases
        zm = np.mean(zp[midx])
        thm = np.mean(th[midx])
    elif select == 'mean':
        # mean emissivity position
        njms = jms/np.sum(jms) 
        Rm = np.sum(Rp * njms)
        zm = np.sum(zp * njms)
        thm = np.sum(th * njms)

    # plt.plot(s,ams)
    # plt.plot(Rp,jms)
    # plt.plot(Rp,jms/jms.max())
    # plt.plot(Rp,F_Te(s)/(1000.0*e))
    # plt.plot(Rp,F_ne(s)/(1e19))
    # plt.axvline(x=Rm, color='k')
    # plt.show()

    # total intensity measured at outisde
    # ece_int = integrate.simps(jms,x=s)
    ece_int = np.sum(jms*ds)

    return ece_int, Rm, zm, thm, s, jms, ams


#ece_int, Rm, zm, thm, s, jms, ams = ece_intensity(np.array([1.7]), np.array([-0.3]), np.array([0.1]), 87.6*1e9*2*np.pi, float(2)) # [m], [m], [rad], [rad/s], harmonic number

# transfer intensity
# Is = np.zeros(R.size)
# for i in range(2,R.size):
#     Is[i] = Is[i-1] + (ams[i-1]*Ibb(Rp[i-1],zp[i-1],omega) - ams[i-1]*Is[i-1])*ds[i]
# %    Is(i) = Is(i-1) + (ams(i-1)*Ibb(Rp(i-1),zp(i-1),omega) - ams(i-1)*Is(i-1))/N2XCsq(Rp(i-1),zp(i-1),theta(i-1))*ds(i); ################## Nsq effect on intensity
# end
# % Is = Is.*N2XCsq(Rp,zp,theta); ################## Nsq effect on intensity
