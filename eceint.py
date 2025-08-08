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

def ece_intensity(s, Rp, zp, th, omega, m, F_B, F_Te, F_ne, select='mean', cache_shape=True): # [m], [m], [rad], [rad/s], harmonic number
    # characteristic arrays along s (vectorized)
    wce_arr = e*F_B(s)/me # [rad/s]
    wpe_arr = np.sqrt(F_ne(s)*e**2/(eps*me)) # [rad/s]
    zeta_arr = mc2/F_Te(s)/2
    mu_arr = (omega/(m*wce_arr))**2

    # refractive index helpers (keep original formulas, evaluate per element as needed)
    X = lambda i: ( wpe_arr[i]/(m*wce_arr[i]) )**2
    if m == 1: # O-mode
        N1OCsq = lambda i,thi: 1 - X(i)*(1-X(i))/(1 - X(i) - np.sin(thi)**2/(2*m**2) + cmath.sqrt((np.sin(thi)**2/(2*m**2))**2 - (1 - X(i))*np.cos(thi)**2/m**2))
        N1OCre = lambda i,thi: (np.lib.scimath.sqrt(N1OCsq(i,thi))).real
        amRz = lambda i,thi: (np.pi/(2*c))*wpe_arr[i]**2*N1OCre(i,thi)*(1 + 2*np.cos(thi)**2)**2*np.sin(thi)**4/(1 + np.cos(thi)**2)**3*F_Te(s[i])/mc2
    elif m == 2: # X-mode
        N2XCsq = lambda i,thi: 1 - X(i)*(1-X(i))/(1 - X(i) - np.sin(thi)**2/(2*m**2) - cmath.sqrt((np.sin(thi)**2/(2*m**2))**2 - (1 - X(i))*np.cos(thi)**2/m**2)) # H(5.2.48)
        N2XCre = lambda i,thi: (np.lib.scimath.sqrt(N2XCsq(i,thi))).real
        am_core = lambda i: e**2*F_ne(s[i])/(4*c*me*eps) * m**(2*m-1)/math.factorial(m-1) * (F_Te(s[i])/(2*mc2))**(m-1) # H(5.2.39)
        a2sq = lambda i,thi: (1 + (1 - X(i))*N2XCre(i,thi)**2*np.cos(thi)**2/(1 - X(i) - N2XCre(i,thi)**2*np.sin(thi)**2)**2 \
        * m**2*(1 - (m**2 -1)/m**2/X(i)*(1 - N2XCre(i,thi)**2))**2)**2*np.sin(thi)**2
        b2sq = lambda i,thi: (1 + (1 - X(i))/(1 - X(i) - N2XCre(i,thi)**2*np.sin(thi)**2) \
        * m**2*(1 - (m**2 -1)/m**2/X(i)*(1 - N2XCre(i,thi)**2))**2)**2*np.cos(thi)**2
        eta2X = lambda i,thi: N2XCre(i,thi)**(2*m-3)*(m - 1)**2*(1 - (m+1)/m/X(i)*(1 - N2XCre(i,thi)**2))**2 \
        / ((1 + np.cos(thi)**2)*(a2sq(i,thi) + b2sq(i,thi))**(1.0/2.0)) # H(5.2.47)
        amRz = lambda i,thi: am_core(i) * eta2X(i,thi) # H(5.2.54)
    else:
        # Default to O-mode-like handling
        N1OCsq = lambda i,thi: 1 - X(i)*(1-X(i))/(1 - X(i) - np.sin(thi)**2/(2*m**2) + cmath.sqrt((np.sin(thi)**2/(2*m**2))**2 - (1 - X(i))*np.cos(thi)**2/m**2))
        N1OCre = lambda i,thi: (np.lib.scimath.sqrt(N1OCsq(i,thi))).real
        amRz = lambda i,thi: (np.pi/(2*c))*wpe_arr[i]**2*N1OCre(i,thi)*(1 + 2*np.cos(thi)**2)**2*np.sin(thi)**4/(1 + np.cos(thi)**2)**3*F_Te(s[i])/mc2

    # shape Maxwellian (relativistic + Doppler) helpers; cache to avoid repeating quad for similar (mu, th, zeta)
    def _ba1(mu, thi):
        return (mu*np.cos(thi) - np.lib.scimath.sqrt(1 - mu*np.sin(thi)**2))/(1 + mu*np.cos(thi)**2)
    def _ba2(mu, thi):
        return (mu*np.cos(thi) + np.lib.scimath.sqrt(1 - mu*np.sin(thi)**2))/(1 + mu*np.cos(thi)**2)

    shape_cache = {} if cache_shape else None
    def _shape_integral(mu, thi, zeta):
        if cache_shape:
            # quantize keys to reduce duplicates while keeping accuracy
            key = (round(float(mu), 3), round(float(thi), 3), round(float(zeta), 3))
            cached = shape_cache.get(key)
            if cached is not None:
                return cached
        ba1 = _ba1(mu, thi)
        ba2 = _ba2(mu, thi)
        # Guard: if limits are complex or invalid, return 0
        if not np.isfinite(ba1) or not np.isfinite(ba2):
            val = 0.0
        else:
            integrand = lambda ba: (1 - ba**2.0 - (1 - ba*np.cos(thi))**2.0*mu)**2.0 \
                * (1 - ba*np.cos(thi)) * np.exp(-zeta*(1 - (1 - ba*np.cos(thi))**2.0*mu))
            val = integrate.quad(integrand, ba1, ba2, epsabs=1e-6, epsrel=1e-6, limit=100)[0]
        if cache_shape:
            shape_cache[key] = val
        return val

    # perfect blackbody intensity factor per s
    Ibb_pref = (omega/(2*np.pi*c))**2

    # absorption coefficent along beam path s (from hfs to lfs)
    n = s.size
    ams = np.zeros(n)

    for i in range(n):
        thi = th[i]
        mu_i = mu_arr[i]
        # skip unphysical region (integral blows up with imaginary beta)
        if mu_i*np.sin(thi)**2 < 1:
            # shape factor (scalar)
            zeta_i = zeta_arr[i]
            m_wce = (m*wce_arr[i])
            if m_wce == 0:
                shp = 0.0
            else:
                shp_pref = 2.0*np.pi*zeta_i**3.5*omega/(np.sqrt(np.pi)*(m_wce)**2.0)
                shp = shp_pref * _shape_integral(mu_i, thi, zeta_i)
            ams[i] = (amRz(i, thi)).real * shp
        else:
            ams[i] = 0.0

    # path differentials and cumulative s (use provided s to integrate)
    ds = np.zeros(Rp.size)
    for i in range(1,Rp.size):
        ds[i] = np.sqrt((Rp[i] - Rp[i-1])**2 + (zp[i] - zp[i-1])**2)

    # differential optical depth dtau
    tau_total = integrate.trapz(ams,x=s)
    dtau = tau_total - np.append(0,integrate.cumtrapz(ams,x=s))

    # emissivity after reabsorption
    jms = ams * (Ibb_pref*F_Te(s)) * np.exp(-dtau)
   
    # max or mean emissivity position
    if select == 'max': 
        midx = np.where(jms == jms.max())
        Rm = np.mean(Rp[midx])
        zm = np.mean(zp[midx])
        thm = np.mean(th[midx])
    else: # 'mean'
        njms = jms/np.sum(jms) if np.sum(jms) != 0 else np.zeros_like(jms)
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

    # total intensity measured at outside
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
