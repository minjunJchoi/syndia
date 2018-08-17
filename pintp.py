import numpy as np
import math
import cmath
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from scipy import interpolate

from pfunc import *

from time import strftime

# M.J. Choi (mjchoi@nfri.re.kr)
# CC BY-NC-SA

# all mks units except Te
# Rp : R coordinates on beam path [m]
# zp : z coordinate on beam path [m]
# th : angle between field and beam path [rad]
# omega : emission frequency [rad/s]
# m : harmonic number
# B : B field function in R-z coordinate [T]
# F_Te : 2d TriScatteredInterp Te function in R-z coordinates [J]
# F_ne : 2d TriScatteredInterp ne function in R-z coordinates [m^-3]
# @(x,y) : (R,z) coordinate

e = 1.602*1e-19
me = 9.109*1e-31
eps = 8.854*1e-12
c = 299792458
mc2 = me*c**2

def intp_prof(Rp, zp, th, step): # [m], [m], [rad], [rad/s], harmonic number
    Rps = Rp[::step]
    zps = zp[::step]
    B_path = np.zeros(Rps.size)
    Te_path = np.zeros(Rps.size)
    ne_path = np.zeros(Rps.size)
    for i in range(Rps.size):
        B_path[i] = F_B(Rps[i], zps[i])
        Te_path[i] = F_Te(Rps[i], zps[i])
        ne_path[i] = F_ne(Rps[i], zps[i])
    #F_Bs = interpolate.interp2d(Rps, zps, B_path, kind='linear')
    #F_nes = interpolate.interp2d(Rps, zps, ne_path, kind='linear')
    #F_Tes = interpolate.interp2d(Rps, zps, Te_path, kind='linear')
    Btck = interpolate.bisplrep(Rps, zps, B_path)
    F_Bs = lambda R,z: interpolate.bisplev(R,z,Btck)
    Tetck = interpolate.bisplrep(Rps, zps, Te_path)
    F_Tes = lambda R,z: interpolate.bisplev(R,z,Tetck)
    netck = interpolate.bisplrep(Rps, zps, ne_path/1e19)
    F_nes = lambda R,z: interpolate.bisplev(R,z,netck)*1e19


    plt.plot(Rps, ne_path, 'go-')
    plt.plot(Rps, F_nes(Rps, zps))
    plt.show()

    return F_Bs, F_Tes, F_nes
