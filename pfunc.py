"""
Plasma profile functions
Repo : https://github.com/minjunJchoi/syndia
Author : Minjun J. Choi (mjchoi@nfri.re.kr)
Collaborators : Jieun Lee, Yoonbum Nam
Acknowledgements : Tongryeol Rhee
"""
from geqdsk_dk import geqdsk_dk as geqdsk_dk
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

# geqdsk_fn = "data/g013728.003900"
# ne_fn = "data/ne_3900.dat" [norm. psi, 1e19 m^-3]
# Te_fn = "data/Te_3900.dat" [norm. psi, keV]
## set 2D coordinate
# Rmin = 1.2 # [m]
# Rmax = 2.4 # [m]
# zmin = -0.5 # [m]
# zmax = 0.5 # [m]
# ds = 0.01 # [m]
# R1d = np.arange(Rmin, Rmax, ds)
# z1d = np.arange(zmin, zmax, ds)
# R2d, z2d = np.meshgrid(R1d, z1d)
# R2d[0,:]
# z2d[:,0] # -> row vector; if you want a column vector, z2d[:,[0]]

class ProfFunc(object):
    def __init__(self, geqdsk_fn, Te_fn, ne_fn, bfactor=1.0):
        ## B = B*bfactor; considering B-field uncertainty 
        self.bfactor = bfactor

        ## read equilibrium data
        # geqdsk file
        self.geq = geqdsk_dk(filename=geqdsk_fn)
        # density data: psin_ne, ne [10^-19 m^-3]; well behaving data
        with open(ne_fn, 'r') as f:
            self.psin_ne, self.ne = np.loadtxt(f, unpack=True)
        # temperature data: psin_Te, Te [keV]; well behaving data
        with open(Te_fn, 'r') as f:
            self.psin_Te, self.Te = np.loadtxt(f, unpack=True)

        print('geqdsk file = {}'.format(geqdsk_fn))
        print('Te file = {}'.format(Te_fn))
        print('ne file = {}'.format(ne_fn))

        ## Interpolators 
        # 1) psi(R,z) interpolator from geqdsk (used in F_psin/F_ne/F_Te)
        self._psi_interp = self.geq.get_psi_normal()
        # 2) ne(psin) and Te(psin) 1D interpolators (used in F_ne/F_Te)
        self._ne_interp = interpolate.interp1d(self.psin_ne, self.ne, kind='linear', bounds_error=False, fill_value=(self.ne[0], 0.01))
        self._Te_interp = interpolate.interp1d(self.psin_Te, self.Te, kind='linear', bounds_error=False, fill_value=(self.Te[0], 0.01))
        # 3) (normalized perturbations) ne_pert((R, z, time)) and Te_pert((R, z, time)) 2D interpolators
        self._pert_time = 0
        self._Te_pert_interp = None
        self._ne_pert_interp = None

    # B = f(R, z) [T] ######################## different interpolation may result in difference
    def F_B(self, R, z):
        if type(R) is np.ndarray:
            B = np.zeros(R.shape)
            it = np.nditer(R, flags=['multi_index'], op_flags=['readonly'])
            while not it.finished:
                idx = it.multi_index
                B[idx] = self.geq.B_abs(R[idx], z[idx])*self.bfactor
                it.iternext()
        else:
            B = float(self.geq.B_abs(R, z))*self.bfactor
        return B

    # Bvec = [Br, Bz, Bt] [T]
    def F_Bvec(self, R, z):
        if type(R) is np.ndarray:
            Bvec = np.zeros((3,) + R.shape)
            it = np.nditer(R, flags=['multi_index'], op_flags=['readonly'])
            while not it.finished:
                idx = it.multi_index
                B_field = self.geq.B_field(R[idx], z[idx])
                Bvec[0,idx] = B_field[0]
                Bvec[1,idx] = B_field[1]
                Bvec[2,idx] = B_field[2]
                it.iternext()
        else:
            Bvec = self.geq.B_field(R, z)
        # Bvec = self.geq.B_fields(R, z)
        return Bvec

    # psin = f(R, z) ## different interpolation may result in difference
    def F_psin(self, R, z):
        return self._psi_interp((R, z))

    # ne [m^-3] = f(psin)
    def F_ne_psin(self, psin):
        val = self._ne_interp(psin)
        return val*1e19

    # delta ne [normalized] = f(R, z)
    def F_ne_pert(self, R, z):
        # In normalized unit total ne_pert = (ne - ne_0) / ne_0
        time = self._pert_time * np.ones_like(R)  # Perturbation time is constant for all R, z
        return self._ne_pert_interp((R, z, time)) if self._ne_pert_interp is not None else 0.0

    # ne [m^-3] = f(R, z) [m, m]
    def F_ne(self, R, z):
        # Equilibrium density
        ne = self.F_ne_psin(self.F_psin(R, z))

        # Total density
        ne = ne + ne * self.F_ne_pert(R, z)  # Add normalized perturbation
        return ne

    # Te [J] = f(psin)
    def F_Te_psin(self, psin):
        val = self._Te_interp(psin)
        return val*1000*1.602*1e-19

    # delta Te [normalized] = f(R, z)
    def F_Te_pert(self, R, z):
        # In normalized unit total Te_pert = (Te - Te_0) / Te_0
        time = self._pert_time * np.ones_like(R)  # Perturbation time is constant for all R, z
        return self._Te_pert_interp((R, z, time)) if self._Te_pert_interp is not None else 0.0

    # Te [J] = f(R, z) [m, m]
    def F_Te(self, R, z):
        # Equilibrium temperature
        Te = self.F_Te_psin(self.F_psin(R, z))

        # Total temperature
        Te = Te + Te * self.F_Te_pert(R, z)  # Add normalized perturbation
        return Te


# R = np.arange(1.8,2.0,0.01)
# z = np.arange(0.1,0.3,0.01)
# print F_B(R,z)
# Bvec = F_Bvec(R,z)
# print Bvec.shape
# B = np.sqrt(np.sum(Bvec * Bvec, axis=0))
# print B

# print np.sqrt(F_Bvec(R[0],z[0]).dot(F_Bvec(R[0],z[0])))
# print F_Bvec(R[0],z[0])
# print np.sqrt(F_Bvec(R,z).dot(F_Bvec(R,z)))
# print F_ne(R,z)
# print F_Te(R,z)
# psin = np.arange(0.1,0.2,0.01)
# print F_Te_psin(psin)
