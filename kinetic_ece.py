"""
General-distribution ECE emission and absorption.

For ANY electron momentum distribution (distributions.py) this module returns, at one point and one frequency,
    j(omega)      emission coefficient divided by N_r^2   [W m^-3 sr^-1 Hz^-1]   (the transported quantity)
    alpha(omega)  absorption coefficient                   [1/m]
for the cold-plasma O- or X-mode, summed over harmonics, with
  * fully relativistic resonance  gamma - N_par u_par = m Y       (Y = omega_ce/omega, N_par = N cos(theta))
  * exact Bessel functions (no small-Larmor-radius expansion)
  * the mode polarisation e and refractive index N from cold-plasma theory (Stix)
  * absorption from detailed balance with the photon momentum (identical to the anti-Hermitian part of the
    relativistic dielectric tensor):  alpha = -((2 pi c)^2/omega^3) Int d^3u eta D F,
        D F = (omega/m_e c^2) [ (m Y/u_perp) dF/du_perp + N_par dF/du_par ]   on the resonance
    so that for a Juttner distribution j = alpha * I_BB(T) holds exactly (Kirchhoff), I_BB = (omega/2 pi c)^2 T.

Single-particle kernel (per Hz, already divided by N_r^2):
    eta = 2 pi * e^2 omega^2/(8 pi^2 eps0 c) * (1/N) |e* . b_m|^2/|e_t|^2 * delta(omega(1 - N_par beta_par) - m Omega)
    b_m = ( beta_perp (J_{m-1}+J_{m+1})/2 ,  i beta_perp J_m' ,  beta_par J_m ),   argument  xi = N sin(theta) u_perp / Y
  frame: z || B, k = k (sin theta, 0, cos theta); e_t = transverse part of e (Poynting flux ~ N |e_t|^2).
  In vacuum (N = 1, transverse e) this is Bekefi's single-particle formula.
Resonance integral:  Int d^3u delta(...) G = (2 pi/omega) Int du_par gamma^2 G  on  u_perp^2 = (mY + N_par u_par)^2 - 1 - u_par^2.

Numerics: the resonance ellipse is integrated with Simpson's rule in t, u_par = u_c + u_s sinh(t), centred on the
peak of the emission integrand (found on a coarse grid), so a narrow thermal core and a long tail (kappa) are both
resolved.  Constants are SYNDIA's (eceint.py) so that Kirchhoff holds to rounding error inside SYNDIA.
"""
import math
import numpy as np
from scipy.special import jv

try:
    from compat import simpson
except ImportError:
    from scipy import integrate as _i
    simpson = getattr(_i, 'simpson', None) or getattr(_i, 'simps')

e = 1.602*1e-19
me = 9.109*1e-31
eps = 8.854*1e-12
c = 299792458
mc2 = me*c**2

MODES = ('O', 'X')


# ------------------------------------------------------------------------------------------------
# cold plasma: refractive index and polarisation
# ------------------------------------------------------------------------------------------------
def cold_N2(X, Y, theta, mode):
    """Appleton-Hartree N^2 (electrons only).  mode 'O' (+ sign) or 'X' (- sign)."""
    s2, c2 = np.sin(theta)**2, np.cos(theta)**2
    Delta = np.sqrt(0.25*Y**4*s2**2 + Y**2*(1 - X)**2*c2)
    sgn = 1.0 if mode == 'O' else -1.0
    den = 1 - X - 0.5*Y**2*s2 + sgn*Delta
    return 1 - X*(1 - X)/den


def cold_polarization(X, Y, theta, N2):
    """Unit polarisation vector (complex, |e| = 1) of the cold-plasma mode with refractive index^2 N2,
    frame z || B, k in the x-z plane, time dependence exp(-i omega t) (Stix).  Null vector of
    [[S - N^2 c^2, -iD, N^2 s c], [iD, S - N^2, 0], [N^2 s c, 0, P - N^2 s^2]]."""
    s, cth = np.sin(theta), np.cos(theta)
    if abs(1 - Y**2) < 1e-9:          # exactly at the fundamental: cold S, D diverge; take the limit from below
        Y = 1 - 1e-9
    S = 1 - X/(1 - Y**2)
    D = -X*Y/(1 - Y**2)
    P = 1 - X
    M = np.array([[S - N2*cth**2, -1j*D, N2*s*cth],
                  [1j*D, S - N2, 0],
                  [N2*s*cth, 0, P - N2*s**2]], dtype=complex)
    _, _, vh = np.linalg.svd(M)
    ev = vh[-1].conj()
    k = np.argmax(np.abs(ev))
    ev = ev*np.exp(-1j*np.angle(ev[k]))
    return ev/np.linalg.norm(ev)


def mode_props(omega, B, ne, theta, mode):
    """(N, e, e_t^2) for the cold-plasma mode; N = 0 when evanescent."""
    X = ne*e**2/(eps*me)/omega**2
    Y = e*B/me/omega
    N2 = cold_N2(X, Y, theta, mode)
    if not np.isfinite(N2) or N2 <= 0:
        return 0.0, None, 0.0
    ev = cold_polarization(X, Y, theta, N2)
    kh = np.array([np.sin(theta), 0.0, np.cos(theta)])
    et = ev - np.dot(kh, ev)*kh
    return float(np.sqrt(N2)), ev, float(np.real(np.vdot(et, et)))


def I_BB(omega, T_J):
    """Rayleigh-Jeans black body per Hz per sr (SYNDIA's convention) [W m^-2 sr^-1 Hz^-1]."""
    return (omega/(2*np.pi*c))**2*T_J


# ------------------------------------------------------------------------------------------------
# one harmonic
# ------------------------------------------------------------------------------------------------
def _resonance_range(m, Y, Npar):
    a = 1 - Npar**2
    disc = (m*Y)**2 - a
    if a <= 0 or disc <= 0:
        return None
    sq = np.sqrt(disc)
    return (m*Y*Npar - sq)/a, (m*Y*Npar + sq)/a


def _kernel(m, Y, Npar, Nb, sth, ev, upar, need_alpha, dist):
    """Integrand pieces on the resonance at the points upar: returns (gamma^2 |e*.b|^2 F, gamma^2 |e*.b|^2 DF')
    with DF' = (m Y/u_perp) dF/du_perp + N_par dF/du_par (the factor omega/mc^2 is applied by the caller)."""
    gam = m*Y + Npar*upar
    uperp = np.sqrt(np.clip(gam**2 - 1 - upar**2, 0, None))
    bperp, bpar = uperp/gam, upar/gam
    xi = Nb*sth*uperp/Y
    Jm, Jm1, Jp1 = jv(m, xi), jv(m - 1, xi), jv(m + 1, xi)
    bx = bperp*0.5*(Jm1 + Jp1)          # beta_perp m J_m(xi)/xi
    by = 1j*bperp*0.5*(Jm1 - Jp1)       # i beta_perp J_m'(xi)
    bz = bpar*Jm
    proj = np.conj(ev[0])*bx + np.conj(ev[1])*by + np.conj(ev[2])*bz
    w = gam**2*np.abs(proj)**2
    Fv = dist.F(uperp, upar)
    gj = w*Fv
    if not need_alpha:
        return gj, None
    dperp, dpar = dist.dF(uperp, upar)
    with np.errstate(divide='ignore', invalid='ignore'):
        rat = np.where(uperp > 0, dperp/uperp, 0.0)
    ga = w*(m*Y*rat + Npar*dpar)
    return gj, ga


def _harmonic(m, Y, Npar, Nb, sth, ev, dist, need_alpha=True, n_coarse=401, n_fine=601):
    rng = _resonance_range(m, Y, Npar)
    if rng is None:
        return 0.0, 0.0
    lo, hi = rng
    g0, _ = _kernel(m, Y, Npar, Nb, sth, ev, np.linspace(lo, hi, 101), False, dist)
    if not np.any(g0 > 0):                  # quick exit: nothing on this resonance (e.g. far from a thermal line)
        return 0.0, 0.0
    xs = np.linspace(lo, hi, n_coarse)
    gj, _ = _kernel(m, Y, Npar, Nb, sth, ev, xs, False, dist)
    gmax = np.max(gj)
    if not np.isfinite(gmax) or gmax <= 0:
        return 0.0, 0.0
    i = int(np.argmax(gj))
    sig = np.where(gj >= 1e-3*gmax)[0]
    width = max(xs[sig[-1]] - xs[sig[0]], (hi - lo)/n_coarse)
    uc, us = xs[i], max(width/20.0, 1e-9)
    t = np.linspace(np.arcsinh((lo - uc)/us), np.arcsinh((hi - uc)/us), n_fine)
    up = uc + us*np.sinh(t)
    up[0], up[-1] = lo, hi
    jac = us*np.cosh(t)
    gj, ga = _kernel(m, Y, Npar, Nb, sth, ev, up, need_alpha, dist)
    Ij = simpson(gj*jac, x=t)
    Ia = simpson(ga*jac, x=t) if need_alpha else 0.0
    return Ij, Ia


# ------------------------------------------------------------------------------------------------
# public
# ------------------------------------------------------------------------------------------------
def emission_absorption(dist, omega, B, ne, theta, mode, harmonics=None, n_harm=3, N_res='cold',
                        need_alpha=True, n_fine=601, return_info=False):
    """Emission (j/N^2, per Hz) and absorption coefficient at one point.
      dist     : distributions.Distribution (normalised)
      omega    : [rad/s];  B [T];  ne [m^-3];  theta = angle between k and B [rad]
      mode     : 'O' or 'X'
      harmonics: iterable of m; default = the lowest resonant harmonic and the next n_harm-1
      N_res    : 'cold' -> N_par = N cos(theta) in the resonance;  'vacuum' -> N_par = cos(theta) (as SYNDIA's
                 line shape).  The Bessel argument, flux factor 1/N and polarisation always use the cold N.
    Returns (j, alpha) or (j, alpha, info)."""
    if mode not in MODES:
        raise ValueError('mode must be O or X')
    N, ev, et2 = mode_props(omega, B, ne, theta, mode)
    info = {'N': N, 'harmonics': {}}
    if N == 0 or et2 <= 0:
        return (0.0, 0.0, info) if return_info else (0.0, 0.0)
    Y = e*B/me/omega
    Npar = (N if N_res == 'cold' else 1.0)*np.cos(theta)
    sth = np.sin(theta)
    if harmonics is None:
        if Npar**2 >= 1:
            raise ValueError('N_par >= 1: resonance curve is open; not supported')
        m0 = max(1, int(math.ceil(np.sqrt(1 - Npar**2)/Y - 1e-12)))
        harmonics = range(m0, m0 + n_harm)
    K = 2*np.pi*e**2*omega**2/(8*np.pi**2*eps*c)/N/et2      # per Hz, includes 1/N (j already / N^2)
    j = alpha = 0.0
    for m in harmonics:
        Ij, Ia = _harmonic(m, Y, Npar, N, sth, ev, dist, need_alpha, n_fine=n_fine)
        jm = ne*K*(2*np.pi/omega)*Ij
        am = -((2*np.pi*c)**2/omega**3)*ne*K*(2*np.pi/omega)*(omega/mc2)*Ia
        info['harmonics'][m] = (jm, am)
        j += jm
        alpha += am
    return (j, alpha, info) if return_info else (j, alpha)


def transport(s, j, alpha):
    """Radiative transfer along s (from the far side to the observer at s[-1]) with SYNDIA's discretisation:
    tau(s) = optical depth between s and the observer; returns (I, tau_total, j*exp(-tau))."""
    try:
        from compat import trapezoid, cumulative_trapezoid
    except ImportError:
        from scipy.integrate import trapezoid, cumulative_trapezoid
    s, j, alpha = (np.asarray(v, dtype=float) for v in (s, j, alpha))
    dtau = trapezoid(alpha, x=s) - np.append(0, cumulative_trapezoid(alpha, x=s))
    ds = np.r_[0, np.diff(s)]
    jms = j*np.exp(-dtau)
    return float(np.sum(jms*ds)), float(trapezoid(alpha, x=s)), jms
