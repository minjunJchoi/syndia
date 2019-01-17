"""
ECE intensity calculation
Repo : https://github.com/minjunJchoi/syndia
Author : Minjun J. Choi (mjchoi@nfri.re.kr)
Collaborators : Jieun Lee, Yoonbum Nam
References : Hutchinson, Bornatici, Rathegeber
"""
import numpy as np
import math
import cmath
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from scipy import interpolate

#from pfunc import *

from time import strftime

#[Rm, zm, s, tau, jms, theta_max, Iece] = ece_intensity(Rp, zp, th, Rc, omega, m, F_B, F_Te, F_ne)
# M.J. Choi (mjchoi@nfri.re.kr)
# CC BY-NC-SA

# all mks units except Te
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

# it will calculate emission profile of frequency omega along s
# Rm : maximum emission position [m]
# zm : maximum emission position [m]
# Iece : ece intensity [W/m^2 Hz]
## all mks units, Temperature in [J]


e = 1.602*1e-19
me = 9.109*1e-31
eps = 8.854*1e-12
c = 299792458
mc2 = me*c**2

#[Rm, zm, s, tau, jms, theta_max, Iece] = ece_intensity(Rp, zp, th, Rc, omega, m, F_B, F_Te, F_ne)

# for verification with matlab code
#F_B = lambda s: 3.16804
#F_Te = lambda s: 1.81824e-16
#F_ne = lambda s: 2.08039e+19

def ece_intensity(s, Rp, zp, th, omega, m, F_B, F_Te, F_ne): # [m], [m], [rad], [rad/s], harmonic number
    m = float(m)

    # characteristic frequencies
    wce = lambda s: e*F_B(s)/me # [rad/s]
    wpe = lambda s: np.sqrt(F_ne(s)*e**2/(eps*me)) # [rad/s]

    # characteristic functions
    zeta = lambda s: mc2/F_Te(s)/2
    mu = lambda s,w: (w/(m*wce(s)))**2

    # refractive index
    X = lambda s: ( wpe(s)/(m*wce(s)) )**2
    if m == 1: # O-mode
    #     Nperp2C = @(x,y,w) 1 - (wpe(x,y)./w).^2; % B(3.1.22)
        N1OCsq = lambda s,th: 1 - X(s)*(1-X(s))/(1 - X(s) - np.sin(th)**2/(2*m**2) + cmath.sqrt((np.sin(th)**2/(2*m**2))**2 - (1 - X(s))*np.cos(th)**2/m**2))
        N1OCre = lambda s,th: (np.lib.scimath.sqrt(N1OCsq(s,th))).real
    elif m == 2: # X-mode
    #     Nperp2C = @(x,y,w) 1 - (wpe(x,y)./w).^2.*(w.^2 - wpe(x,y).^2)./(w.^2 - wce(x,y).^2 - wpe(x,y).^2); % B(3.1.12)
        N2XCsq = lambda s,th: 1 - X(s)*(1-X(s))/(1 - X(s) - np.sin(th)**2/(2*m**2) - cmath.sqrt((np.sin(th)**2/(2*m**2))**2 - (1 - X(s))*np.cos(th)**2/m**2)) # H(5.2.48)
        N2XCre = lambda s,th: (np.lib.scimath.sqrt(N2XCsq(s,th))).real

    # absorption coefficient
    if m == 1: # O-mode
        amO = lambda s,th: np.pi/(2*c)*wpe(s)**2*N1OCre(s,th)*(1 + 2*np.cos(th)**2)**2*np.sin(th)**4/(1 + np.cos(th)**2)**3*F_Te(s)/mc2 # H(5.2.52)

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

    # num integration results in some difference (compared to matlab code)
    shape = lambda s,w,th: 2.0*np.pi*zeta(s)**(3.5)*w/(np.sqrt(np.pi)*(m*wce(s))**2.0) \
        *integrate.quad(lambda ba: (1 - ba**2.0 - (1 - ba*np.cos(th))**2.0*mu(s,w))**2.0 \
        *(1 - ba*np.cos(th))*np.exp(-zeta(s)*(1 - (1 - ba*np.cos(th))**2.0*mu(s,w))), \
        ba1(s,w,th), ba2(s,w,th))[0]

    # perfect blackbody intensity
    Ibb = lambda s,w: (w/(2*np.pi*c))**2*F_Te(s) # H(5.2.37)

    ## start calculation
    # absorption coefficent along beam path s (from hfs to lfs)
    ams = np.zeros(s.size)
    for i in range(s.size):
        if mu(s[i], omega)*np.sin(th[i])**2 < 1: # the integral blows up with the imaginary beta
            ams[i] = (amRz(s[i], th[i])).real*(shape(s[i], omega, th[i])).real
        else:
            ams[i] = 0

    # # define path from the inside
    # s = np.zeros(Rp.size)
    ds = np.zeros(Rp.size)
    for i in range(1,Rp.size):
        ds[i] = np.sqrt((Rp[i] - Rp[i-1])**2 + (zp[i] - zp[i-1])**2)
    #     s[i] = s[i-1] + ds[i]

    # calculate differential optical depth dtau
    tau = integrate.trapz(ams,x=s) - np.append(0,integrate.cumtrapz(ams,x=s)) # from zero; a single point matters

    # emissitivty after reabsorption
    jms = np.zeros(s.size)
    for i in range(s.size):
        jms[i] = ams[i]*Ibb(s[i],omega)*np.exp(-tau[i]) # emissivity after reabsorption. B(2.2.13), B(2.2.15)

    #plt.plot(s,ams)
    #plt.plot(Rp,jms)
    
    # maximum emissivity position
    #midx = np.where(jms == jms.max())
    #Rm = np.mean(Rp[midx]) # for multiple maximum cases
    #zm = np.mean(zp[midx])
    #thm = np.mean(th[midx])
    
    #plt.plot(Rm, jms[midx], 'x')

    # mean emissivity position
    njms = jms/np.sum(jms) 
    Rm = np.sum(Rp * njms)
    zm = np.sum(zp * njms)
    thm = np.sum(th * njms)

    #plt.plot(Rm, jms[midx], 'o')

    #plt.show()

    # total intensity measured at outisde
    # ece_int = integrate.simps(jms,x=s)
    ece_int = np.sum(jms*ds)

    return ece_int, Rm, zm, thm, s, jms, ams, tau


#ece_int, Rm, zm, thm, s, jms, ams = ece_intensity(np.array([1.7]), np.array([-0.3]), np.array([0.1]), 87.6*1e9*2*np.pi, float(2)) # [m], [m], [rad], [rad/s], harmonic number

# transfer intensity
# Is = np.zeros(R.size)
# for i in range(2,R.size):
#     Is[i] = Is[i-1] + (ams[i-1]*Ibb(Rp[i-1],zp[i-1],omega) - ams[i-1]*Is[i-1])*ds[i]
# %    Is(i) = Is(i-1) + (ams(i-1)*Ibb(Rp(i-1),zp(i-1),omega) - ams(i-1)*Is(i-1))/N2XCsq(Rp(i-1),zp(i-1),theta(i-1))*ds(i); ################## Nsq effect on intensity
# end
# % Is = Is.*N2XCsq(Rp,zp,theta); ################## Nsq effect on intensity
