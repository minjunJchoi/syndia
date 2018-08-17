from collections import OrderedDict
import subprocess
from scipy import interpolate
import math
import h5py
import numpy as np
import matplotlib.pyplot as plt

TB_path = "/home/users/mjchoi/torbeam_ifortran/"
TB_path_run = "/home/users/mjchoi/torbeam_ifortran/run_torbeam/"
ECEI_data_path = '/eceidata/exp_2015/'

e = 1.602*1e-19
me = 9.109*1e-31

VNT = 24

def torbeam_prof(geqdsk_fn, pf):
    ## save for TORBEAM
    # read geqdsk file from the selected EFIT run time and convert it to topfile for TORBEAM
    args = "{}readeqdsk<{}".format(TB_path_eq, geqdsk_fn)
    re = subprocess.check_output(args, shell=True)
    args = "mv topfile {}".format(TB_path_run)
    re = subprocess.check_output(args, shell=True)

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

    print 'profile data for TORBEAM are saved in {}'.format(TB_path_run)


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

        print 'inbeam file was written in {}'.format(TB_path_run)



def run_torbeam(hn, freq, ainit, zinit, Rinit):
    ## TORBEAM
    # initial parameters
    write_inbeam(hn, freq*1e9, ainit/np.pi*180, zinit*100.0, Rinit*100.0)
    # run
    args = "{}run.sh".format(TB_path)
    re = subprocess.check_output(args, shell=True)
    # obtain output
    with open("{}t1_LIB.dat".format(TB_path_run), 'r') as f:
        Rp, zp = np.loadtxt(f, usecols=(0,1), unpack=True)
        Rp = Rp/100 # R [m] beam path, from R_vec_end to core; decreasing
        zp = zp/100 # z [m] beam path
        # it can be curved and get double values. Find the first resonance position

    return Rp, zp


def set_beam_path(Rp, zp, hn, freq, pstart, pend, pint, pf):
    # functions needed maybe
    wce = lambda R,z: e*pf.F_B(R,z)/me # [rad/s]

    # find the proper range
    for i in range(Rp.size): # increasing during loop
        fRz = wce(Rp[i], zp[i])/(2*np.pi*1e9)*hn # EC frequency [GHz]

        if np.abs(freq + pend - fRz) < 0.3:
            idx1 = i # no need to be very accurate
        if np.abs(freq - fRz) < 0.3:
            Rcold = Rp[i] # EC resonance position [cm] no need to be accurate
        if np.abs(freq + pstart - fRz) < 0.3:
            idx2 = i # no need to be very accurate
            break

    # Rp, zp between idx1 idx2; calculate angle between emission direction and B-field
    Rp = Rp[idx1:(idx2+1)]
    zp = zp[idx1:(idx2+1)]
    theta = np.zeros(Rp.size)
    for i in range(1,Rp.size):
        Rvec = np.array([-(Rp[i]-Rp[i-1]), -(zp[i]-zp[i-1]), 0]) # opposite direction for emission path
        Bvec = pf.F_Bvec(Rp[i], zp[i])
        theta[i] = math.acos( Bvec.dot(Rvec) / ( np.sqrt(Bvec.dot(Bvec)) * np.sqrt(Rvec.dot(Rvec)) ) ) # [rad]
    theta[0] = theta[1] + (theta[1]-theta[2])

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


def vac_beam_path(shot, dev, rpos, vn, dz):
    # IN : shot, device name, R posistion [m], vertical channel number
    # OUT : a ray vertical position and angle at rpos [m] [rad]
    # this will find a ray vertical position and angle at rpos [m]
    # ray starting from the array box posistion

    if shot < 19392:
        fname = "{:s}{:06d}/ECEI.{:06d}.{:s}FS.h5".format(ECEI_data_path, shot, shot, dev)
    else:
        fname = "{:s}{:06d}/ECEI.{:06d}.{:s}.h5".format(ECEI_data_path, shot, shot, dev)

    with h5py.File(fname, 'r') as f:
        # get attributes
        dset = f['ECEI']
        sf = dset.attrs['LensFocus']
        sz = dset.attrs['LensZoom']

        rpos = rpos*1000  # [m] -> [mm] # for rpos = [1600:50:2300]

        # ABCD matrix for LFS, HFS, GFS
        if dev == 'L' or dev == 'GR':
            sp = 3350 - rpos

            abcd = np.array([[1,250+sp],[0,1]]).dot(
                   np.array([[1,0],[(1.52-1)/(-730),1.52]])).dot(
                   np.array([[1,135],[0,1]])).dot(
                   np.array([[1,0],[(1-1.52)/(2700*1.52),1/1.52]])).dot(
                   np.array([[1,1265-sz],[0,1]])).dot(
                   np.array([[1,0],[(1.52-1)/1100,1.52]])).dot(
                   np.array([[1,40],[0,1]])).dot(
                   np.array([[1,0],[(1-1.52)/(-1100*1.52),1/1.52]])).dot(
                   np.array([[1,sz],[0,1]])).dot(
                   np.array([[1,0],[0,1.52]])).dot(
                   np.array([[1,65],[0,1]])).dot(
                   np.array([[1,0],[(1-1.52)/(800*1.52),1/1.52]])).dot(
                   np.array([[1,710-sf+140],[0,1]])).dot(
                   np.array([[1,0],[(1.52-1)/(-1270),1.52]])).dot(
                   np.array([[1,90],[0,1]])).dot(
                   np.array([[1,0],[(1-1.52)/(1270*1.52),1/1.52]])).dot(
                   np.array([[1,539+35+sf],[0,1]]))
        elif dev == 'H' or dev == 'HT':
            sp = 3350 - rpos

            abcd = np.array([[1,250+sp],[0,1]]).dot(
                   np.array([[1,0],[(1.52-1)/(-730),1.52]])).dot(
                   np.array([[1,135],[0,1]])).dot(
                   np.array([[1,0],[(1-1.52)/(2700*1.52),1/1.52]])).dot(
                   np.array([[1,1265-sz],[0,1]])).dot(
                   np.array([[1,0],[(1.52-1)/1100,1.52]])).dot(
                   np.array([[1,40],[0,1]])).dot(
                   np.array([[1,0],[(1-1.52)/(-1100*1.52),1/1.52]])).dot(
                   np.array([[1,sz],[0,1]])).dot(
                   np.array([[1,0],[0,1.52]])).dot(
                   np.array([[1,65],[0,1]])).dot(
                   np.array([[1,0],[(1-1.52)/(800*1.52),1/1.52]]))
            if shot > 12297:  # since 2015 campaign
                abcd = abcd.dot(
                   np.array([[1,520-sf+590-9.2],[0,1]])).dot(
                   np.array([[1,0],[(1.52-1)/(-1100),1.52]])).dot(
                   np.array([[1,88.4],[0,1]])).dot(
                   np.array([[1,0],[(1-1.52)/(1100*1.52),1/1.52]])).dot(
                   np.array([[1,446+35+sf-9.2],[0,1]]))
            else:
                abcd = abcd.dot(
                   np.array([[1,520-sf+590],[0,1]])).dot(
                   np.array([[1,0],[(1.52-1)/(-1400),1.52]])).dot(
                   np.array([[1,70],[0,1]])).dot(
                   np.array([[1,0],[(1-1.52)/(1400*1.52),1/1.52]])).dot(
                   np.array([[1,446+35+sf],[0,1]]))
        elif dev == 'G' or dev == 'GT':
            sp = 3150 - rpos

            abcd = np.array([[1,1350-sz+sp],[0,1]]).dot(
                   np.array([[1,0],[0,1.545]])).dot(
                   np.array([[1,100],[0,1]])).dot(
                   np.array([[1,0],[(1-1.545)/(900*1.545),1/1.545]])).dot(
                   np.array([[1,1430-sf+660+sz+470],[0,1]])).dot(
                   np.array([[1,0],[0,1.545]])).dot(
                   np.array([[1,70],[0,1]])).dot(
                   np.array([[1,0],[(1-1.545)/(800*1.545),1/1.545]])).dot(
                   np.array([[1,sf-470],[0,1]])).dot(
                   np.array([[1,0],[0,1.545]])).dot(
                   np.array([[1,80],[0,1]])).dot(
                   np.array([[1,0],[(1-1.545)/(800*1.545),1/1.545]])).dot(
                   np.array([[1,390],[0,1]]))

        # vertical position from the reference axis (vertical center of all lens, z=0 line) at ECEI array box
        zz = (np.arange(VNT,0,-1) - 12.5)*14 + dz # [mm]
        # angle against the reference axis at ECEI array box
        aa = np.zeros(np.size(zz))

        # vertical posistion and angle at rpos
        za = np.dot(abcd, [zz, aa])
        zpos = za[0][vn-1]/1000  # zpos [m]
        apos = za[1][vn-1]  # angle [rad] positive means the (z+) up-directed (divering from array to plasma)

        return zpos, apos
