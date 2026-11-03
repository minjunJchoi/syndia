"""
Version-compatibility layer for SYNDIA.

One code base for both older and newer scientific-Python stacks:
  * old   : NumPy 1.x,  SciPy < 1.14  (trapz / cumtrapz / simps exist)
  * new   : NumPy 2.x,  SciPy >= 1.14 (those names were removed; trapezoid / cumulative_trapezoid / simpson)
  * Python >= 3.12 has no `imp` module (use importlib)
  * MDSplus is imported lazily, so SYNDIA can be imported and run offline
    (e.g. with explicitly supplied ECEI operation parameters, or on a machine
    without the MDSplus client) -- the server is contacted only when actually needed.
Tested with NumPy 1.23.5/SciPy 1.9.3, NumPy 1.26.4/SciPy 1.13.1 and NumPy 2.4/SciPy 1.17 (Python 3.11).
"""
import numpy as np
from scipy import integrate as _integrate

# --- SciPy integration routines (renamed in 1.6, old names removed in 1.14) ---
trapezoid = getattr(_integrate, 'trapezoid', None) or getattr(_integrate, 'trapz')
cumulative_trapezoid = getattr(_integrate, 'cumulative_trapezoid', None) or getattr(_integrate, 'cumtrapz')
simpson = getattr(_integrate, 'simpson', None) or getattr(_integrate, 'simps')


def scalar(x):
    """Return a Python float from a scalar, a 0-d array or a size-1 array.
    NumPy >= 2 refuses to assign a (1,1) array (e.g. from RectBivariateSpline(x, y))
    into a float array element; NumPy 1.x silently took the single value."""
    a = np.asarray(x)
    if a.size != 1:
        raise TypeError('scalar(): expected a single value, got shape {}'.format(a.shape))
    return float(a.reshape(-1)[0])


def module_available(name):
    """Replacement for the removed imp.find_module(name)."""
    import importlib.util
    return importlib.util.find_spec(name) is not None


class LazyMDSConnection(object):
    """Behaves like MDSplus.Connection(server), but connects on first use."""
    def __init__(self, server):
        self._server = server
        self._conn = None

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(name)
        if self._conn is None:
            from MDSplus import Connection   # imported only when the server is really needed
            self._conn = Connection(self._server)
        return getattr(self._conn, name)
