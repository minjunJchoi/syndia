"""
Intepolated plasma profile along beam path
Repo : https://github.com/minjunJchoi/syndia 
Author : Minjun J. Choi (mjchoi@nfri.re.kr)
Collaborators : Jieun Lee, Yoonbum Nam
Acknowledgements : Tongryeol Rhee 
"""
import numpy as np
import math
import cmath
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from scipy import interpolate

from time import strftime

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

def intp_prof(Rp, zp, th, pf, step): # [m], [m], [rad]
    # define path from the inside s=0
    s = np.zeros(Rp.size)
    ds = np.zeros(Rp.size)
    for i in range(1,Rp.size):
        ds[i] = np.sqrt((Rp[i] - Rp[i-1])**2 + (zp[i] - zp[i-1])**2)
        s[i] = s[i-1] + ds[i]

    if step > 1:
        Rps = np.append(Rp[::step],Rp[-1])
        zps = np.append(zp[::step],zp[-1])
        ss = np.append(s[::step],s[-1])
    else:
        Rps = Rp
        zps = zp
        ss = s

    Bs = np.zeros(Rps.size)
    Tes = np.zeros(Rps.size)
    nes = np.zeros(Rps.size)
    for i in range(Rps.size):
        Bs[i] = pf.F_B(Rps[i], zps[i])
        Tes[i] = pf.F_Te(Rps[i], zps[i])
        nes[i] = pf.F_ne(Rps[i], zps[i])

    Btck = interpolate.splrep(ss, Bs)
    F_Bs = lambda s: interpolate.splev(s, Btck)
    Tetck = interpolate.splrep(ss, Tes/1e-16)
    F_Tes = lambda s: interpolate.splev(s ,Tetck)*1e-16
    netck = interpolate.splrep(ss, nes/1e19)
    F_nes = lambda s: interpolate.splev(s ,netck)*1e19


    #plt.plot(ss, Tes, 'go-')
    #plt.plot(s, F_Tes(s))
    #plt.show()

    return s, F_Bs, F_Tes, F_nes
