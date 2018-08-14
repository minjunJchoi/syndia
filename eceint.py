import numpy as np
import math
import cmath
import scipy.integrate as integrate

from pfunc import *

#[Rm, zm, s, tau, jms, theta_max, Iece] = ece_intensity(Rp, zp, th, Rc, omega, m, F_B, F_Te, F_ne)
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
#F_B = lambda R,z: 3.16804
#F_Te = lambda R,z: 1.81824e-16
#F_ne = lambda R,z: 2.08039e+19

def ece_intensity(Rp, zp, th, omega, m): # [m], [m], [rad], [rad/s], harmonic number
    m = float(m)

    # characteristic frequencies
    wce = lambda R,z: e*F_B(R,z)/me # [rad/s]
    wpe = lambda R,z: np.sqrt(F_ne(R,z)*e**2/(eps*me)) # [rad/s]

    # characteristic functions
    zeta = lambda R,z: mc2/F_Te(R,z)/2
    mu = lambda R,z,w: (w/(m*wce(R,z)))**2

    # refractive index
    X = lambda R,z: ( wpe(R,z)/(m*wce(R,z)) )**2
    if m == 1: # O-mode
    #     Nperp2C = @(x,y,w) 1 - (wpe(x,y)./w).^2; % B(3.1.22)
        N1OCsq = lambda R,z,th: 1 - X(R,z)*(1-X(R,z))/(1 - X(R,z) - np.sin(th)**2/(2*m**2) + cmath.sqrt((np.sin(th)**2/(2*m**2))**2 - (1 - X(R,z))*np.cos(th)**2/m**2))
        N1OCre = lambda R,z,th: (np.lib.scimath.sqrt(N1OCsq(R,z,th))).real
    elif m == 2: # X-mode
    #     Nperp2C = @(x,y,w) 1 - (wpe(x,y)./w).^2.*(w.^2 - wpe(x,y).^2)./(w.^2 - wce(x,y).^2 - wpe(x,y).^2); % B(3.1.12)
        N2XCsq = lambda R,z,th: 1 - X(R,z)*(1-X(R,z))/(1 - X(R,z) - np.sin(th)**2/(2*m**2) - cmath.sqrt((np.sin(th)**2/(2*m**2))**2 - (1 - X(R,z))*np.cos(th)**2/m**2)) # H(5.2.48)
        N2XCre = lambda R,z,th: (np.lib.scimath.sqrt(N2XCsq(R,z,th))).real

    # absorption coefficient
    if m == 1: # O-mode
        amO = lambda R,z,th: np.pi/(2*c)*wpe(R,z)**2*N1OCre(R,z,th)*(1 + 2*np.cos(th)**2)**2*np.sin(th)**4/(1 + np.cos(th)**2)**3*F_Te(R,z)/mc2 # H(5.2.52)

        amRz = lambda R,z,th: amO(R,z,th)
    elif m == 2: # X-mode
        am = lambda R,z: e**2*F_ne(R,z)/(4*c*me*eps) * m**(2*m-1)/math.factorial(m-1) * (F_Te(R,z)/(2*mc2))**(m-1) # H(5.2.39)

        a2sq = lambda R,z,th: (1 + (1 - X(R,z))*N2XCre(R,z,th)**2*np.cos(th)**2/(1 - X(R,z) - N2XCre(R,z,th)**2*np.sin(th)**2)**2 \
        * m**2*(1 - (m**2 -1)/m**2/X(R,z)*(1 - N2XCre(R,z,th)**2))**2)**2*np.sin(th)**2
        
        b2sq = lambda R,z,th: (1 + (1 - X(R,z))/(1 - X(R,z) - N2XCre(R,z,th)**2*np.sin(th)**2) \
        * m**2*(1 - (m**2 -1)/m**2/X(R,z)*(1 - N2XCre(R,z,th)**2))**2)**2*np.cos(th)**2

        eta2X = lambda R,z,th: N2XCre(R,z,th)**(2*m-3)*(m - 1)**2*(1 - (m+1)/m/X(R,z)*(1 - N2XCre(R,z,th)**2))**2 \
        / ((1 + np.cos(th)**2)*(a2sq(R,z,th) + b2sq(R,z,th))**(1.0/2.0)) # H(5.2.47)

        amX = lambda R,z,th: am(R,z) * eta2X(R,z,th) # H(5.2.54)

        amRz = lambda R,z,th: amX(R,z,th)

    # shape Maxwellian (relativistic + Doppler)
    ba1 = lambda R,z,w,th: (mu(R,z,w)*np.cos(th) - np.lib.scimath.sqrt(1 - mu(R,z,w)*np.sin(th)**2))/(1 + mu(R,z,w)*np.cos(th)**2)
    ba2 = lambda R,z,w,th: (mu(R,z,w)*np.cos(th) + np.lib.scimath.sqrt(1 - mu(R,z,w)*np.sin(th)**2))/(1 + mu(R,z,w)*np.cos(th)**2)

    # num integration results in some difference (compared to matlab code)
    shape = lambda R,z,w,th: 2.0*np.pi*zeta(R,z)**(3.5)*w/(np.sqrt(np.pi)*(m*wce(R,z))**2.0) \
        *integrate.quad(lambda ba: (1 - ba**2.0 - (1 - ba*np.cos(th))**2.0*mu(R,z,w))**2.0 \
        *(1 - ba*np.cos(th))*np.exp(-zeta(R,z)*(1 - (1 - ba*np.cos(th))**2.0*mu(R,z,w))), \
        ba1(R,z,w,th), ba2(R,z,w,th))[0]   

    # perfect blackbody intensity
    Ibb = lambda R,z,w: (w/(2*np.pi*c))**2*F_Te(R,z) # H(5.2.37)

    ## start calculation
    # absorption coefficent along beam path s (from hfs to lfs)
    ams = np.zeros(Rp.size)
    for i in range(Rp.size):
        if mu(Rp[i], zp[i], omega)*np.sin(th[i])**2 < 1: # the integral blows up with the imaginary beta
            ams[i] = (amRz(Rp[i], zp[i], th[i])).real*(shape(Rp[i], zp[i], omega, th[i])).real
        else:
            ams[i] = 0

    # define path from the inside
    s = np.zeros(Rp.size)
    ds = np.zeros(Rp.size)
    for i in range(1,Rp.size):
        ds[i] = np.sqrt((Rp[i] - Rp[i-1])**2 + (zp[i] - zp[i-1])**2)
        s[i] = s[i-1] + ds[i]

    # calculate differential optical depth dtau
    tau = integrate.trapz(ams,x=s) - integrate.cumtrapz(ams,x=s)

    # emissitivty after reabsorption
    jms = np.zeros(Rp.size)
    for i in range(Rp.size-1): # tau.size = Rp.size -1 
        jms[i] = ams[i]*Ibb(Rp[i],zp[i],omega)*np.exp(-tau[i]) # emissivity after reabsorption. B(2.2.13), B(2.2.15)

    # total intensity measured at outisde
    ece_int = integrate.simps(jms,x=s)

    # maximum emissivity position
    midx = np.where(jms == jms.max())
    Rm = np.mean(Rp[midx]) # for multiple maximum cases
    zm = np.mean(zp[midx])
    thm = np.mean(th[midx])

    return ece_int, Rm, zm, thm, s, jms, ams


#ece_int, Rm, zm, thm, s, jms, ams = ece_intensity(np.array([1.7]), np.array([-0.3]), np.array([0.1]), 87.6*1e9*2*np.pi, float(2)) # [m], [m], [rad], [rad/s], harmonic number

# transfer intensity
# Is = np.zeros(R.size)
# for i in range(2,R.size):
#     Is[i] = Is[i-1] + (ams[i-1]*Ibb(Rp[i-1],zp[i-1],omega) - ams[i-1]*Is[i-1])*ds[i]
# %    Is(i) = Is(i-1) + (ams(i-1)*Ibb(Rp(i-1),zp(i-1),omega) - ams(i-1)*Is(i-1))/N2XCsq(Rp(i-1),zp(i-1),theta(i-1))*ds(i); ################## Nsq effect on intensity
# end
# % Is = Is.*N2XCsq(Rp,zp,theta); ################## Nsq effect on intensity
