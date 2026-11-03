"""
Electron momentum distributions for the general-distribution ECE kernel (kinetic_ece.py).

Conventions
-----------
  u = p / (m_e c)   (normalised momentum, u = gamma*beta),  gamma = sqrt(1 + u^2)
  Theta = T / (m_e c^2)
  Every distribution is normalised to ONE in d^3u:  Int F d^3u = 1   (the density n_e is applied by the kernel).

Only the GYRO-AVERAGED distribution F(u_perp, u_par) enters incoherent ECE emission and absorption
(the single-particle emission of a helical orbit does not depend on the gyrophase, which is only a time
shift).  A distribution that depends on the gyrophase -- e.g. Grad's 13-moment expansion with an
off-diagonal stress pi_xy u_x u_y in the plane perpendicular to B -- is therefore reduced to its gyro-average
here.  For a distribution given in an arbitrary (lab) frame, the gyro-average is taken numerically around the
LOCAL field direction b, so the parts of the stress tensor / heat flux that project onto b survive.

Interface (what kinetic_ece needs from a distribution object)
  F(uperp, upar)            gyro-averaged distribution, vectorised; must be even in uperp
  dF(uperp, upar)           (dF/duperp, dF/dupar); a central-difference default is provided
  scale                     characteristic momentum (e.g. sqrt(Theta)); used for numerical steps
All classes are plain module-level classes with array attributes, so they pickle (SYNDIA's parallel mode).

Classes
  Juttner(Theta)                      relativistic Maxwellian (Kirchhoff exact)
  BiMaxwell(Theta_perp, t_r)          exp(-(sqrt(1 + u_perp^2 + u_par^2/t_r) - 1)/Th_perp), normalisable for all t_r
  Kappa(Theta, kappa, form)           [1 + E/((kappa-3/2) Theta)]^-(kappa+1), E = gamma-1 ('relativistic') or u^2/2 ('classic')
  Gyrotropic(func, scale)             user function F(uperp, upar)
  Tabulated(uperp, upar, table)       table on a grid (e.g. Fokker-Planck code output), spline interpolation
  Arbitrary3D(f, b, scale)            user function f(u1, u2, u3) in any Cartesian frame; numerical gyro-average about b
  Grad13(Theta, b, pi_hat, q_hat)     Grad 13-moment expansion about a Juttner distribution (an Arbitrary3D)
"""
import warnings
import numpy as np
from scipy import integrate
from scipy.special import kve
from scipy.interpolate import RectBivariateSpline

try:
    from compat import simpson
except ImportError:  # stand-alone use
    simpson = getattr(integrate, 'simpson', None) or getattr(integrate, 'simps')


def _gamma(uperp, upar):
    return np.sqrt(1.0 + uperp**2 + upar**2)


class Distribution(object):
    """Base class.  Subclasses implement F(uperp, upar) and set self.scale."""
    scale = 0.05
    gyrotropic = True
    name = 'distribution'

    def F(self, uperp, upar):
        raise NotImplementedError

    def dF(self, uperp, upar):
        """Central differences (step 1e-4 * scale); subclasses with analytic derivatives override this."""
        h = 1e-4*max(self.scale, 1e-4)
        uperp = np.asarray(uperp, dtype=float)
        upar = np.asarray(upar, dtype=float)
        dperp = (self.F(uperp + h, upar) - self.F(uperp - h, upar))/(2*h)
        dpar = (self.F(uperp, upar + h) - self.F(uperp, upar - h))/(2*h)
        return dperp, dpar

    def norm(self, umax=None, n=401):
        """Numerical check of Int F d^3u (should be 1)."""
        umax = umax or 14*self.scale
        up = np.linspace(0, umax, n)
        ul = np.linspace(-umax, umax, 2*n - 1)
        P, L = np.meshgrid(up, ul, indexing='ij')
        val = 2*np.pi*P*self.F(P.ravel(), L.ravel()).reshape(P.shape)
        tot = simpson(simpson(val, x=ul, axis=1), x=up)
        # tail check: a non-normalisable F (or a box that is too small) shows up as weight at the edge
        edge = np.abs(val[-1, :]).max()*umax + np.abs(val[:, [0, -1]]).max()*umax
        if np.isfinite(tot) and tot != 0 and edge > 1e-3*abs(tot):
            warnings.warn('{}: {:.1e} of the normalisation integral sits at the edge of the box (umax = {:.3g}); '
                          'F may not be normalisable or the box is too small'.format(self.name, edge/abs(tot), umax),
                          RuntimeWarning)
        return tot


# ------------------------------------------------------------------------------------------------
class Juttner(Distribution):
    """Relativistic Maxwellian, F = exp(-(gamma-1)/Theta) / Z,  Z = 4 pi Theta K2(1/Theta) e^(1/Theta)."""
    name = 'juttner'

    def __init__(self, Theta):
        self.Theta = float(Theta)
        self.scale = np.sqrt(self.Theta)
        self.Z = 4*np.pi*self.Theta*kve(2, 1/self.Theta)

    def F(self, uperp, upar):
        g = _gamma(uperp, upar)
        return np.exp(-(g - 1)/self.Theta)/self.Z

    def dF(self, uperp, upar):
        g = _gamma(uperp, upar)
        Fv = np.exp(-(g - 1)/self.Theta)/self.Z
        return -uperp/(g*self.Theta)*Fv, -upar/(g*self.Theta)*Fv


# ------------------------------------------------------------------------------------------------
def _radial_norm(fE, Theta, rel=True):
    """4 pi Int u^2 f(E(u)) du for an isotropic f, E = gamma-1 (rel) or u^2/2."""
    E = (lambda u: np.sqrt(1 + u*u) - 1) if rel else (lambda u: 0.5*u*u)
    s = np.sqrt(Theta)
    g = lambda u: 4*np.pi*u*u*fE(E(u))
    a = integrate.quad(g, 0, 20*s, limit=400)[0]
    b = integrate.quad(g, 20*s, np.inf, limit=400)[0]
    return a + b


class Kappa(Distribution):
    """Kappa distribution with the SAME second moment (kinetic temperature T) as a Maxwellian in the
    non-relativistic limit:
        F ~ [1 + E / ((kappa - 3/2) Theta)]^-(kappa+1)
        form='relativistic' : E = gamma - 1   (tail ~ u^-(kappa+1); normalisable for kappa > 2)
        form='classic'      : E = u^2/2       (the usual space-physics form written in u; kappa > 3/2)
    kappa -> infinity gives the Juttner (relativistic) / Maxwellian (classic) distribution."""
    name = 'kappa'

    def __init__(self, Theta, kappa, form='relativistic'):
        if form not in ('relativistic', 'classic'):
            raise ValueError("form must be 'relativistic' or 'classic'")
        kmin = 2.0 if form == 'relativistic' else 1.5
        if kappa <= kmin:
            raise ValueError('kappa must be > {} for form={}'.format(kmin, form))
        self.Theta, self.kappa, self.form = float(Theta), float(kappa), form
        self.rel = form == 'relativistic'
        self.scale = np.sqrt(self.Theta)
        self.w = (self.kappa - 1.5)*self.Theta
        self.Z = _radial_norm(lambda E: (1 + E/self.w)**(-(self.kappa + 1)), self.Theta, self.rel)

    def _E(self, uperp, upar):
        if self.rel:
            return _gamma(uperp, upar) - 1
        return 0.5*(uperp**2 + upar**2)

    def F(self, uperp, upar):
        return (1 + self._E(uperp, upar)/self.w)**(-(self.kappa + 1))/self.Z

    def dF(self, uperp, upar):
        E = self._E(uperp, upar)
        Fv = (1 + E/self.w)**(-(self.kappa + 1))/self.Z
        dFdE = -(self.kappa + 1)/(self.w + E)*Fv
        g = _gamma(uperp, upar) if self.rel else 1.0
        return dFdE*uperp/g, dFdE*upar/g

    def T_eff(self, gamma):
        """Local 'slope temperature' -1/(dlnF/deps) [units of m_e c^2] at energy gamma (relativistic form)."""
        return (self.w + (gamma - 1))/(self.kappa + 1)


# ------------------------------------------------------------------------------------------------
class BiMaxwell(Distribution):
    """Relativistic bi-Maxwellian:
        F = exp(-(g~ - 1)/Th_perp) / Z,   g~ = sqrt(1 + u_perp^2 + u_par^2/t_r),   t_r = Th_par/Th_perp,
        Z = sqrt(t_r) * 4 pi Th_perp K2(1/Th_perp) e^(1/Th_perp)   (exact: substitute u_par = sqrt(t_r) w).
    t_r = 1 is exactly the Juttner distribution; the non-relativistic limit is exp(-u_perp^2/2Th_perp - u_par^2/2Th_par).
    The form is normalisable for every t_r > 0.
    (The 'Juttner-core' form exp(-(gamma-1)/Th_perp - (1/Th_par - 1/Th_perp) u_par^2/2) is NOT normalisable for t_r > 1:
     the +u_par^2 term beats the linear growth of gamma, so it is not used.)
    A relativistic bi-Maxwellian is not unique; with this form T_rad = T_perp at 90 deg holds only to O(Theta)
    (it is exact for t_r = 1)."""
    name = 'bimaxwell'

    def __init__(self, Theta_perp, t_r):
        if t_r <= 0:
            raise ValueError('t_r must be positive')
        self.Theta_perp, self.t_r = float(Theta_perp), float(t_r)
        self.Theta_par = self.t_r*self.Theta_perp
        self.scale = np.sqrt(max(self.Theta_perp, self.Theta_par))
        self.Z = np.sqrt(self.t_r)*4*np.pi*self.Theta_perp*kve(2, 1/self.Theta_perp)

    def _gt(self, uperp, upar):
        return np.sqrt(1.0 + uperp**2 + upar**2/self.t_r)

    def F(self, uperp, upar):
        return np.exp(-(self._gt(uperp, upar) - 1)/self.Theta_perp)/self.Z

    def dF(self, uperp, upar):
        gt = self._gt(uperp, upar)
        Fv = np.exp(-(gt - 1)/self.Theta_perp)/self.Z
        return -uperp/(gt*self.Theta_perp)*Fv, -upar/(self.t_r*gt*self.Theta_perp)*Fv


# ------------------------------------------------------------------------------------------------
class Gyrotropic(Distribution):
    """User-supplied gyrotropic distribution F(uperp, upar) (vectorised, even in uperp).
    func must be a module-level function (or picklable object) for SYNDIA's parallel mode.
    dfunc (optional) returns (dF/duperp, dF/dupar); otherwise central differences are used."""
    name = 'gyrotropic'

    def __init__(self, func, scale, dfunc=None, normalize=True, umax=None):
        self.func, self.dfunc, self.scale = func, dfunc, float(scale)
        self.Z = 1.0
        if normalize:
            self.Z = self.norm(umax=umax)

    def F(self, uperp, upar):
        return self.func(uperp, upar)/self.Z

    def dF(self, uperp, upar):
        if self.dfunc is None:
            return Distribution.dF(self, uperp, upar)
        a, b = self.dfunc(uperp, upar)
        return a/self.Z, b/self.Z


class Tabulated(Distribution):
    """F given on a grid uperp (>= 0, increasing) x upar (increasing), e.g. from a Fokker-Planck code.
    Bicubic spline, mirrored to uperp < 0 so that F is even in uperp; negative spline undershoots are
    clipped to zero; zero outside the table.  Normalised numerically unless normalize=False."""
    name = 'tabulated'

    def __init__(self, uperp, upar, table, normalize=True):
        uperp, upar, table = (np.asarray(v, dtype=float) for v in (uperp, upar, table))
        if uperp[0] != 0:
            raise ValueError('uperp grid must start at 0')
        up2 = np.r_[-uperp[:0:-1], uperp]
        tab2 = np.vstack([table[:0:-1], table])
        self.spl = RectBivariateSpline(up2, upar, tab2, kx=3, ky=3)
        self.umax_perp, self.lo, self.hi = uperp[-1], upar[0], upar[-1]
        m = table.max()
        mean2 = np.sum(table*(uperp[:, None]**2 + upar[None, :]**2)*uperp[:, None])/np.sum(table*uperp[:, None])
        self.scale = float(np.sqrt(max(mean2/3, 1e-6)))
        self.Z = 1.0
        if normalize:
            self.Z = self.norm(umax=max(self.umax_perp, abs(self.lo), self.hi), n=401)
        self._fmax = m

    def _inside(self, uperp, upar):
        return (np.abs(uperp) <= self.umax_perp) & (upar >= self.lo) & (upar <= self.hi)

    def F(self, uperp, upar):
        uperp, upar = np.broadcast_arrays(np.asarray(uperp, float), np.asarray(upar, float))
        v = self.spl.ev(uperp, upar)
        return np.where(self._inside(uperp, upar), np.clip(v, 0, None), 0.0)/self.Z

    def dF(self, uperp, upar):
        uperp, upar = np.broadcast_arrays(np.asarray(uperp, float), np.asarray(upar, float))
        ok = self._inside(uperp, upar) & (self.spl.ev(uperp, upar) > 0)
        return (np.where(ok, self.spl.ev(uperp, upar, dx=1), 0.0)/self.Z,
                np.where(ok, self.spl.ev(uperp, upar, dy=1), 0.0)/self.Z)


# ------------------------------------------------------------------------------------------------
def field_basis(b):
    """Orthonormal (e1, e2, b) with b the unit field vector (given in any Cartesian frame)."""
    b = np.asarray(b, dtype=float)
    b = b/np.linalg.norm(b)
    trial = np.array([1.0, 0, 0]) if abs(b[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = trial - np.dot(trial, b)*b
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(b, e1)
    return e1, e2, b


class Arbitrary3D(Distribution):
    """Any distribution f(u1, u2, u3) given in a Cartesian frame in which the local field direction is b.
    It need NOT be gyrotropic or separable (e.g. terms like u1*u2).  The kernel only needs the gyro-average
        F(uperp, upar) = (1/2pi) Int dphi f(uperp cos(phi) e1 + uperp sin(phi) e2 + upar b)
    which is computed with nphi equally spaced angles (spectrally accurate for smooth f).
      f        : vectorised function f(u1, u2, u3) (module-level for parallel runs); extra keyword
                 arguments for f can be passed as fargs (a dict)
      b        : field direction in the same frame (3-vector)
      scale    : characteristic momentum (for numerical derivatives / normalisation grid)
      negative : 'clip' (default) -> negative gyro-averaged values are set to 0 and counted (self.n_negative);
                 'raise'          -> ValueError;  'keep' -> no action
      normalize: numerical normalisation on a (161 x 321) grid; set False if f is already normalised."""
    name = 'arbitrary3d'
    gyrotropic = False

    def __init__(self, f, b, scale, nphi=64, fargs=None, normalize=True, umax=None, negative='clip'):
        self.f, self.fargs, self.scale, self.nphi = f, (fargs or {}), float(scale), int(nphi)
        self.e1, self.e2, self.b = field_basis(b)
        self.phi = 2*np.pi*np.arange(self.nphi)/self.nphi
        self.cphi, self.sphi = np.cos(self.phi), np.sin(self.phi)
        self.negative = negative
        self.n_negative = 0
        self.Z = 1.0
        if normalize:
            self.Z = self.norm(umax=umax or 14*self.scale, n=161)

    def _f_avg(self, uperp, upar):
        uperp, upar = np.broadcast_arrays(np.asarray(uperp, float), np.asarray(upar, float))
        shp = uperp.shape
        up, ul = uperp.reshape(-1, 1), upar.reshape(-1, 1)
        comp = [up*self.cphi*self.e1[i] + up*self.sphi*self.e2[i] + ul*self.b[i] for i in range(3)]
        vals = self.f(comp[0], comp[1], comp[2], **self.fargs)
        return np.mean(vals, axis=1).reshape(shp)

    def F(self, uperp, upar):
        v = self._f_avg(uperp, upar)/self.Z
        if self.negative != 'keep' and np.any(v < 0):
            if self.negative == 'raise':
                raise ValueError('gyro-averaged distribution is negative at {} points'.format(int(np.sum(v < 0))))
            self.n_negative += int(np.sum(v < 0))
            v = np.clip(v, 0, None)
        return v


def _grad13_f(u1, u2, u3, Theta, Z, pi_hat, q_hat):
    """Grad 13-moment expansion about a Juttner distribution, w = u / sqrt(Theta):
        f = f_J(u) [ 1 + (1/2) pi_hat_ij w_i w_j + q_hat_i w_i (w^2/5 - 1) ]
    pi_hat = pi/p (traceless stress over pressure), q_hat = q/(p v_t), v_t^2 = T/m (all dimensionless)."""
    g = np.sqrt(1 + u1*u1 + u2*u2 + u3*u3)
    fJ = np.exp(-(g - 1)/Theta)/Z
    s = 1/np.sqrt(Theta)
    w = (u1*s, u2*s, u3*s)
    w2 = w[0]**2 + w[1]**2 + w[2]**2
    corr = 1.0
    for i in range(3):
        if q_hat[i] != 0:
            corr = corr + q_hat[i]*w[i]*(w2/5 - 1)
        for j in range(3):
            if pi_hat[i][j] != 0:
                corr = corr + 0.5*pi_hat[i][j]*w[i]*w[j]
    return fJ*corr


class Grad13(Arbitrary3D):
    """Grad 13-moment distribution (see _grad13_f), given in a frame where the field direction is b.
    pi_hat : 3x3 traceless symmetric (pi/p), q_hat : 3-vector (q/(p v_t)).
    Already normalised (the corrections integrate to zero), so no numerical normalisation.
    Note: the truncated expansion becomes NEGATIVE at large w (|q_hat| w^3/5 ~ 1), which is exactly where
    relativistically down-shifted ECE samples; negative values are clipped (counted in n_negative)."""
    name = 'grad13'

    def __init__(self, Theta, b, pi_hat=None, q_hat=None, nphi=64, negative='clip'):
        pi_hat = np.zeros((3, 3)) if pi_hat is None else np.asarray(pi_hat, float)
        q_hat = np.zeros(3) if q_hat is None else np.asarray(q_hat, float)
        if abs(np.trace(pi_hat)) > 1e-12 or np.max(np.abs(pi_hat - pi_hat.T)) > 1e-12:
            raise ValueError('pi_hat must be symmetric and traceless')
        self.Theta = float(Theta)
        Z = 4*np.pi*self.Theta*kve(2, 1/self.Theta)
        fargs = dict(Theta=self.Theta, Z=Z, pi_hat=pi_hat.tolist(), q_hat=q_hat.tolist())
        Arbitrary3D.__init__(self, _grad13_f, b, np.sqrt(self.Theta), nphi=nphi, fargs=fargs,
                             normalize=False, negative=negative)

    def projected(self):
        """Field-aligned parts that survive gyro-averaging: pi_par = b.pi.b (pi_perp = -pi_par/2), q_par = q.b."""
        pi = np.asarray(self.fargs['pi_hat'])
        q = np.asarray(self.fargs['q_hat'])
        return float(self.b @ pi @ self.b), float(q @ self.b)


# ================================================================================================
# Distribution MODELS: build the local Distribution at each point of a SYNDIA beam path.
# A model's  at(R, z, Te, ne, Bvec, psin)  gets
#     R, z [m], Te [J] (SYNDIA's Te profile), ne [m^-3], Bvec = (B_R, B_phi, B_Z) [T] (lab frame), psin
# and returns a normalised Distribution.  Models must be picklable (module-level classes / functions,
# profiles stored as arrays) for SYNDIA's parallel mode.
# ================================================================================================
MC2_J = 9.109e-31*299792458**2      # m_e c^2 [J] with SYNDIA's constants


class Profile1D(object):
    """Picklable profile value(psin): constant (float) or linear interpolation of (psin, values),
    or a two-column text file [psin, value]."""

    def __init__(self, v):
        if isinstance(v, str):
            x, y = np.loadtxt(v, unpack=True)
            self.x, self.y, self.const = x, y, None
        elif np.ndim(v) == 0:
            self.const = float(v)
        else:
            self.x, self.y = (np.asarray(a, dtype=float) for a in v)
            self.const = None

    def __call__(self, psin):
        if self.const is not None:
            return self.const
        return float(np.interp(psin, self.x, self.y))


class DistributionModel(object):
    name = 'model'

    def at(self, R, z, Te, ne, Bvec, psin):
        raise NotImplementedError


class MaxwellModel(DistributionModel):
    """Juttner distribution with SYNDIA's Te (the kinetic reference for dist='maxwell')."""
    name = 'maxwell'

    def at(self, R, z, Te, ne, Bvec, psin):
        return Juttner(Te/MC2_J)


class KappaModel(DistributionModel):
    """Kappa distribution with kinetic temperature Te; kappa = float, (psin, kappa) or file."""
    name = 'kappa'

    def __init__(self, kappa, form='relativistic'):
        self.kappa, self.form = Profile1D(kappa), form

    def at(self, R, z, Te, ne, Bvec, psin):
        return Kappa(Te/MC2_J, self.kappa(psin), self.form)


class BiMaxwellModel(DistributionModel):
    """Relativistic bi-Maxwellian (see BiMaxwell), Te = T_perp; t_r = T_par/T_perp = float, (psin, t_r) or file."""
    name = 'bimaxwell'

    def __init__(self, t_r):
        self.t_r = Profile1D(t_r)

    def at(self, R, z, Te, ne, Bvec, psin):
        return BiMaxwell(Te/MC2_J, self.t_r(psin))


def lab_bhat(Bvec):
    """Unit field vector in the lab frame (R, phi, Z) from SYNDIA's (B_R, B_phi, B_Z)."""
    b = np.asarray(Bvec, dtype=float)
    return b/np.linalg.norm(b)


class Grad13Model(DistributionModel):
    """Grad 13-moment distribution about the local Juttner(Te).
      pi_hat : 3x3 traceless symmetric stress / pressure, q_hat : 3-vector heat flux / (p v_t)
      frame  : 'lab'   -> components in (R, phi, Z); the gyro-average is taken about the local b
               'field' -> components in (e1, e2, b) with b the local field (e1, e2 arbitrary)
      scale_profile : optional Profile1D-like (psin -> factor) multiplying pi_hat and q_hat."""
    name = 'grad13'

    def __init__(self, pi_hat=None, q_hat=None, frame='lab', scale_profile=1.0, nphi=64):
        self.pi_hat = np.zeros((3, 3)) if pi_hat is None else np.asarray(pi_hat, float)
        self.q_hat = np.zeros(3) if q_hat is None else np.asarray(q_hat, float)
        if frame not in ('lab', 'field'):
            raise ValueError("frame must be 'lab' or 'field'")
        self.frame, self.scale_profile, self.nphi = frame, Profile1D(scale_profile), nphi

    def at(self, R, z, Te, ne, Bvec, psin):
        b = lab_bhat(Bvec) if self.frame == 'lab' else np.array([0.0, 0.0, 1.0])
        k = self.scale_profile(psin)
        return Grad13(Te/MC2_J, b, k*self.pi_hat, k*self.q_hat, nphi=self.nphi)


class FunctionModel(DistributionModel):
    """Wraps a user function  func(R, z, Te, ne, Bvec, psin, **kwargs) -> Distribution.
    func must be defined at module level (not a lambda) for parallel runs."""
    name = 'function'

    def __init__(self, func, **kwargs):
        self.func, self.kwargs = func, kwargs

    def at(self, R, z, Te, ne, Bvec, psin):
        return self.func(R, z, Te, ne, Bvec, psin, **self.kwargs)
