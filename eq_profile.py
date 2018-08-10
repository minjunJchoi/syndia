# otbain 2D profiles
from geqdsk_dk import geqdsk_dk as geqdsk_dk
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

geqdsk_fn = "data/g013728.003900"
ne_fn = "data/ne_3900.dat"
Te_fn = "data/Te_3900.dat"
Rmin = 1.2 # [m]
Rmax = 2.4 # [m]
zmin = -0.5 # [m]
zmax = 0.5 # [m]
ds = 0.01 # [m]


## set 2D coordinate
R1d = np.arange(Rmin, Rmax, ds)
z1d = np.arange(zmin, zmax, ds)
R2d, z2d = np.meshgrid(R1d, z1d)
# R2d[0,:]
# z2d[:,0] # -> row vector; if you want a column vector, z2d[:,[0]]


## read files
# geqdsk file
geq = geqdsk_dk(filename=geqdsk_fn)
# density data: psin_ne, ne [10^-19 m^-3]; well behaving data
with open(ne_fn, 'r') as f:
    psin_ne, ne = np.loadtxt(f, unpack=True)
# temperature data: psin_Te, Te [keV]; well behaving data
with open(Te_fn, 'r') as f:
    psin_Te, Te = np.loadtxt(f, unpack=True)


# B = f(R, z) ######################## fix ?
def F_B(R, z):
    if type(R) is np.ndarray:
        B = np.zeros(R.shape)
        it = np.nditer(R, flags=['multi_index'], op_flags=['readonly'])
        while not it.finished:
            idx = it.multi_index
            B[idx] = geq.B_abs(R[idx], z[idx])
            it.iternext()
    else:
        B = geq.B_abs(R, z)
    return B
    # return geq.B_abs(R,z)

# psin = f(R, z) ######################## fix ?
def F_psin(R, z):
    f = geq.get_psi_normal() # interpolate.interp2d (return 2d data with two 1d coordinates)
    return f(R, z)

# ne [m^-3] = f(psin)
def F_ne_psin(psin):
    f = interpolate.interp1d(psin_ne, ne, kind='linear')
    return f(psin)*1e19

# ne [m^-3] = f(R, z) [m, m]
def F_ne(R, z):
    if type(R) is np.ndarray:
        ne = np.zeros(R.shape)
        it = np.nditer(R, flags=['multi_index'], op_flags=['readonly'])
        while not it.finished:
            idx = it.multi_index
            ne[idx] = F_ne_psin(F_psin(R[idx], z[idx]))
            it.iternext()
    else:
        ne = F_ne_psin(F_psin(R,z))
    return ne

# Te [J] = f(psin)
def F_Te_psin(psin):
    f = interpolate.interp1d(psin_Te, Te, kind='linear')
    return f(psin)*1000*1.602*1e-19

# Te [J] = f(R, z) [m, m]
def F_Te(R, z):
    if type(R) is np.ndarray:
        Te = np.zeros(R.shape)
        it = np.nditer(R, flags=['multi_index'], op_flags=['readonly'])
        while not it.finished:
            idx = it.multi_index
            Te[idx] = F_Te_psin(F_psin(R[idx], z[idx]))
            it.iternext()
    else:
        Te = F_Te_psin(F_psin(R,z))
    return Te

# R = np.arange(1.8,2.0,0.01)
# z = np.arange(0.1,0.3,0.01)
# print F_B(R,z)
# print F_ne(R,z)
# print F_Te(R,z)
psin = np.arange(0.1,0.2,0.01)
print F_Te_psin(psin)
