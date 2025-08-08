"""
ECE beam path calculation
Repo : https://github.com/minjunJchoi/syndia
Author : Minjun J. Choi (mjchoi@nfri.re.kr)
Collaborators : Jieun Lee, Yoonbum Nam
Acknowledgements : Minho Woo
"""
from collections import OrderedDict
import os
import subprocess
from scipy import interpolate
import math
import h5py
import numpy as np
import matplotlib.pyplot as plt

TB_path_eq = "/home/users/mjchoi/torbeam_ifortran/eqdsk2topfile/"
TB_path_run = "path-to-your-data"

e = 1.602*1e-19
me = 9.109*1e-31
c = 299792458.0 # [m/s]

VNT = 24

def torbeam_prof(geqdsk_fn, pf):
    # set TB path
    global TB_path_run
    TB_path_run = os.path.split(geqdsk_fn)[0] + '/'

    ## save files for TORBEAM
    # read geqdsk file from the selected EFIT run time and convert it to topfile for TORBEAM
    args = "{}readeqdsk<{}".format(TB_path_eq, geqdsk_fn)
    re = subprocess.call(args, shell=True)
    args = "mv topfile {}".format(TB_path_run)
    re = subprocess.call(args, shell=True)

    # save ne.dat for TORBEAM
    x = np.sqrt(pf.psin_ne)
    y = pf.ne
    with open("{}ne.dat".format(TB_path_run), 'w') as f:
        f.write('{}\n'.format(x.size))
        for i in range(x.size):
            f.write('{} {}\n'.format(x[i], y[i]))

    # save Te.dat for TORBEAM
    x = np.sqrt(pf.psin_Te)
    y = pf.Te
    with open("{}Te.dat".format(TB_path_run), 'w') as f:
        f.write('{}\n'.format(x.size))
        for i in range(x.size):
            f.write('{} {}\n'.format(x[i], y[i]))

    print('TORBEAM profile data saved at {}'.format(TB_path_run))


def write_inbeam(nmod, xf, xpoldeg, xzb, xxb):
    with open('{}inbeam.dat'.format(TB_path_run), 'w') as f:
        f.write("&edata\n")

        if nmod == 2:
            nmod = -1

        inbeam = OrderedDict()
        inbeam['xrtol'] = 3e-07 # % required rel. error
        inbeam['xatol'] = 3e-07 # % required abs. error
        inbeam['xstep'] = 1.00 # % integration step
        inbeam['npow'] = 1 # % power absorption on(1)/off(0)
        inbeam['ncd'] = 1 # % current drive calc. on(1)/off(0)
        inbeam['ianexp'] = 2 # % analytic (1) or experimental (2) equilibrium
        inbeam['ndns'] = 2 # % analytic (1) or interpolated (2) density profiles
        inbeam['nte'] = 2 # % analyt. (1) or interp. (2) electron temp. prof.
        inbeam['ncdroutine'] = 2 # % 0->Curba, 1->Lin-Liu, 2->Lin-Liu+momentum conservation
        inbeam['nshot'] = 25485 # % shot number
        inbeam['xtbeg'] = 4.500000 # % tbegin
        inbeam['xtend'] = 4.500000 # % tende
        inbeam['nmaxh'] = 5 # % maximum harmonics
        inbeam['nrela'] = 0 # % weakly relativistic (0) full relativistic (1)
        inbeam['nabsroutine'] = 0 # %
        inbeam['noout'] = 0 # % screen out (0)
        inbeam['nmod'] = nmod # % mode selection: O-mode (1), X-mode (-1)
        inbeam['xf'] = xf # % frequency om=2*pi*xf kstar
        inbeam['xtordeg'] = 0.00000 # % geom. optics injection angle
        inbeam['xpoldeg'] = -xpoldeg # % geom. optics injection angle (- : up, + : down)
        inbeam['xxb'] = xxb # % beam launching position [cm] for tracing calc. kstar
        inbeam['xyb'] = 0.00000 # % beam launching position [cm] for tracing calc. kstar
        inbeam['xzb'] = xzb # % beam launching position [cm] for tracing calc. kstar
        inbeam['xryyb'] = 10000.000 # % initial principal curvature
        inbeam['xrzzb'] = 10000.000 # % initial principal curvature
        inbeam['xwyyb'] = 2.50000 # % initial principal width kstar
        inbeam['xwzzb'] = 2.50000 # % initial principal width kstar
        inbeam['xpw0'] = 1.00000 # % initial beam power [MW]
        inbeam['xrmaj'] = 180.00000 # % major radius
        inbeam['xrmin'] = 50.00000 # % minor radius   ########## parameters below used for the analytic calculation (not necessary) #########
        inbeam['xb0'] = 1.990000 # % central toroidal magnetic field [T']
        inbeam['xdns'] = 2.7e+13 # % electron density [cm**(-3)] - core
        inbeam['edgdns'] = 1.e+12 # % electron density [cm**(-3)] - edge
        inbeam['xe1'] = 2.000000 # % exponent in the density profile - a
        inbeam['xe2'] = 1.000000 # % exponent in the density profile - b
        inbeam['xte0'] = 2.000000 # % electron temperature [keV] - core
        inbeam['xteedg'] = 0.1000000 # % electron temperature [keV] - edge
        inbeam['xe1t'] = 2.000000 # % exponent in the temperature profile - a
        inbeam['xe2t'] = 1.000000 # % exponent in the temperature profile - b
        inbeam['xdel0'] = 0.000000 # % Shafranov shift on the axis
        inbeam['xdeled'] = 0.000000 # % Shafranov shift on the edge
        inbeam['xelo0'] = 1.000000 # % elongation on the axis
        inbeam['xeloed'] = 1.000000 # % elongation on the edge kstar
        inbeam['xq0'] = 1.000000 # % safety factor on the axis
        inbeam['xqedg'] = 2.000000 # % safety factor on the edge

        for key in inbeam:
            if key == 'xqedg':
                f.write("{} = {:g}\n".format(key, inbeam[key]))
                f.write("/\n")
            else:
                f.write("{} = {:g},\n".format(key, inbeam[key]))

        print('TORBEAM inbeam.dat written at {}'.format(TB_path_run))



def run_torbeam(hn, freq, ainit, zinit, Rinit):
    ## TORBEAM
    # initial parameters
    write_inbeam(hn, freq*1e9, ainit/np.pi*180, zinit*100.0, Rinit*100.0)
    # run
    args = "cd {} && ./run_torbeam.out".format(TB_path_run)
    re = subprocess.call(args, shell=True)
    # obtain output
    with open("{}t1_LIB.dat".format(TB_path_run), 'r') as f:
        Rp, zp = np.loadtxt(f, usecols=(0,1), unpack=True)
        Rp = Rp/100 # R [m] beam path, from R_vec_end to core; decreasing
        zp = zp/100 # z [m] beam path
        # it can be curved and get double values. Find the first resonance position

    return Rp, zp


def ray_tracing(hn, freq, ainit, zinit, Rinit, pf):
    ds = 0.005 # 5 mm grad step : similar to TORBEAM
    dt = 1.0/(freq*1e9) # time step

    omega = 2*np.pi*freq*1e9 # [rad/s]

    # Small cache for repeated (R,z) evaluations to reduce expensive pf calls
    # Quantize to 1e-4 m (~0.1 mm) grid to improve hit rate without losing accuracy
    _cache = {}
    def _q(v):
        return round(float(v), 4)
    def _get_wpe2(R, z):
        key = ('wpe2', _q(R), _q(z))
        val = _cache.get(key)
        if val is None:
            val = (5.64e4)**2*(pf.F_ne(R, z)*1.0e-6)
            _cache[key] = val
        return val
    def _get_wce2(R, z):
        key = ('wce2', _q(R), _q(z))
        val = _cache.get(key)
        if val is None:
            val = (1.76e7*pf.F_B(R, z)*1e4)**2
            _cache[key] = val
        return val

    # ray tracing for each vertical channel
    if hn == 1: # O-mode
        denom = lambda R,z: 1.0
        numer = lambda R,z: _get_wpe2(R,z)/omega**2.0
    elif hn == 2: # X-mode
        denom = lambda R,z: 1.0 + _get_wpe2(R,z)*_get_wce2(R,z)/(omega**2.0 - _get_wpe2(R,z) - _get_wce2(R,z))**2
        numer = lambda R,z: _get_wpe2(R,z)/omega**2.0*(omega**2.0 - _get_wpe2(R,z))/(omega**2.0 - _get_wpe2(R,z) - _get_wce2(R,z))

    dnumerdr = lambda R,z: (numer(R+ds,z) - numer(R-ds,z))/(2.0*ds)
    dnumerdz = lambda R,z: (numer(R,z+ds) - numer(R,z-ds))/(2.0*ds)

    dkRdt = lambda R,z: -omega/2.0*dnumerdr(R,z)/denom(R,z)
    dkzdt = lambda R,z: -omega/2.0*dnumerdz(R,z)/denom(R,z)

    l = 1
    nmax = 1000

    Rp_list = [Rinit]
    zp_list = [zinit]

    R = float(Rinit) # [m]
    z = float(zinit) # [m]

    kR = -omega/c*math.cos(ainit)
    kz = omega/c*math.sin(ainit)

    # bounds
    Rmin, Rmax = 1.265, 2.4
    zmin, zmax = -0.49, 0.49

    while l <= nmax and Rmin < R < Rmax and zmin < z < zmax:
        if l > 1:
            Rp_list.append(R)
            zp_list.append(z)

        dRdt = (c**2.0/omega)*kR/denom(R,z)
        dzdt = (c**2.0/omega)*kz/denom(R,z)

        R = R + dRdt*dt
        z = z + dzdt*dt

        kR += dkRdt(R,z)*dt
        kz += dkzdt(R,z)*dt

        l += 1

    Rp = np.array(Rp_list)
    zp = np.array(zp_list)

    return Rp, zp


def set_beam_path(Rp, zp, hn, freq, pstart, pend, pint, pf, verbose=False):
    # functions needed maybe
    wce = lambda R,z: e*pf.F_B(R,z)/me # [rad/s]

    # find the proper range
    fRz = np.zeros(Rp.size)
    for i in range(Rp.size): # increasing during loop
        if i > 0 and (Rp[i] > Rp[i-1]): # no turning back
            break
        fRz[i] = wce(Rp[i], zp[i])/(2*np.pi*1e9)*hn # EC frequency [GHz]

    ridx = np.where((fRz >= freq+pend) & (fRz <= freq+pstart))
    idx1 = ridx[0][0]
    idx2 = ridx[0][-1]

    # Rp, zp between idx1 idx2 from lfs to hfs; calculate angle between emission direction and B-field
    Rp = Rp[idx1:(idx2+1)]
    zp = zp[idx1:(idx2+1)]
    theta = np.zeros(Rp.size)
    for i in range(1,Rp.size):
        Rvec = np.array([-(Rp[i]-Rp[i-1]), -(zp[i]-zp[i-1]), 0]) # opposite direction for emission path
        Bvec = pf.F_Bvec(Rp[i], zp[i])
        theta[i] = math.acos( Bvec.dot(Rvec) / ( np.sqrt(Bvec.dot(Bvec)) * np.sqrt(Rvec.dot(Rvec)) ) ) # [rad]
    theta[0] = theta[1] + (theta[1]-theta[2])

    if verbose:
        # Rcold
        cidx = np.where(np.abs(fRz - freq) == np.abs(fRz - freq).min())
        Rcold = Rp[cidx]
        zcold = zp[cidx]
        print('Rcold, Rpath ends = ', Rcold, Rp[0], Rp[-1])

    # interpolation (for better accuracy) and change direction from hfs to lfs
    idx = np.arange(Rp.size-1,-1,-1)
    nidx = np.arange(Rp.size-1,0,-pint)
    fRp = interpolate.interp1d(idx, Rp[idx], kind='linear')
    fzp = interpolate.interp1d(idx, zp[idx], kind='linear')
    fth = interpolate.interp1d(idx, theta[idx], kind='linear')
    Rp = fRp(nidx)
    zp = fzp(nidx)
    theta = fth(nidx)

    return Rp, zp, theta


def vac_beam_path(ecei, rpos, vn, dz):
    # IN :  R posistion [m], vertical channel number
    # OUT : a ray vertical position and angle at rpos [m] [rad]
    # this will find a ray vertical position and angle at rpos [m]
    # ray starting from the array box posistion

    abcd = ecei.get_abcd(ecei.sf, ecei.sz, rpos)

    # vertical position from the reference axis (vertical center of all lens, z=0 line) at ECEI array box
    zz = (np.arange(VNT,0,-1) - 12.5)*14 + dz # [mm]
    # angle against the reference axis at ECEI array box
    aa = np.zeros(np.size(zz))

    # vertical posistion and angle at rpos
    za = np.dot(abcd, [zz, aa])
    zpos = za[0][vn-1]/1000  # zpos [m]
    apos = za[1][vn-1]  # angle [rad] positive means the (z+) up-directed (divering from array to plasma)

    return zpos, apos
