# ece intensity forward modeling
import numpy as np
import h5py
import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt

from bpath import vac_beam_path, tb_beam_path
from eceint import ece_intensity
from pfunc import *
# def radiation_temperature(self, shot, dev, clist):

e = 1.602*1e-19
me = 9.109*1e-31
eps = 8.854*1e-12
c = 299792458
mc2 = me*c**2

shot = 13728

Nz = 1 # number of vertical rays for a single channel
zstart = 0 # [mm] first ray at the mini lens -14
zend = 0 # [mm] last ray at the mini lens 14
Nf = 1 # number of frequency rays for a single channel
fstart = 0 # [GHz] frequency bandwidth -0.3
fend = 0 # [GHz] frequency bandwidth +0.3

pstart = 7.8 # [GHz] cal start point (hfs) = cold resonance + pstart
pend = -2 # [GHz] cal end point (lfs) = cold resonance + pend
pint = 0.1 # ECE integration path inter step. 0.1 = 10%

Rinit = 2.39 # where vacuum region ends and plasma region starts

## ECEI specific for each shot
data_path = '/eceidata/exp_2015/'
dev = 'G'
clist = ['ECEI_G0101']

Lcz = 9 # e^2 fallding practical vertical width for minilens [mm]
Bcf = 0.3 # [GHz]

if shot < 19392:
    fname = "{:s}{:06d}/ECEI.{:06d}.{:s}FS.h5".format(data_path, shot, shot, dev)
    cnidx1 = 6
else:
    fname = "{:s}{:06d}/ECEI.{:06d}.{:s}.h5".format(data_path, shot, shot, dev)
    cnidx1 = 7

with h5py.File(fname, 'r') as f:
    # get attributes
    dset = f['ECEI']
    mode = dset.attrs['Mode'].strip()
    if mode is 'O':
        hn = 1  # harmonic number
    elif mode is 'X':
        hn = 2
    lo = dset.attrs['LoFreq']

###################### run tbprof; ready for TORBEAM

## loop for channels
cnum = len(clist)
int_meas = np.zeros(cnum) # ECE intensity measured;
# you use this to make a 'synthetic' ECE image from simulation data to compare with the measured ECE image
rad_temp = np.zeros(cnum) # radiation temperature [keV] of each channel;
# it is same with real abs_temp for black body;
# you determine correctness of temperature/density profile (from Thomson, etc) (variable inputs of syndia) by comparing this output with what 'well calibrated' ECE gives you
abs_temp = np.zeros(cnum) # abs_temp from F_Te (input of syndia)
Rch = np.zeros(cnum)  # R [m] of each channel
zch = np.zeros(cnum)  # z [m] of each channel
ach = np.zeros(cnum)  # angle [rad] of each channel
dz = np.linspace(zstart, zend, Nz) # dz [mm] of sub z rays at minilens
for cn in range(0, cnum):
    # ECEI channel
    vn = int(clist[cn][(cnidx1):(cnidx1+2)])
    fn = int(clist[cn][(cnidx1+2):(cnidx1+4)])

    # define sub rays
    fsub = np.linspace((fn-1)*0.9 + 2.6 + lo + fstart, (fn-1)*0.9 + 2.6 + lo + fend, Nf) # frequency [GHz] of sub rays
    zsub = np.zeros(dz.size)
    asub = np.zeros(dz.size)
    S = 0
    ## loop over sub rays of a single channel
    for i in range(dz.size):
        # vacuum approximation until Rinit for ECEI
        zsub[i], asub[i] = vac_beam_path(shot, dev, Rinit, vn, dz[i]) # vertical position [m] and rangle [rad] of sub rays at Rinit

        for j in range(fsub.size):
            # find beam path
            Rp, zp, theta = tb_beam_path(hn, fsub[j], asub[i], zsub[i], Rinit, pstart, pend, pint) # [GHz], [rad], [m], [m]

            # calculate ECE intensity along path
            ece_int, Rm, zm, thm, s, jms, ams, tau = ece_intensity(Rp, zp, theta, 2*np.pi*fsub[j]*1e9, hn) # [m], [m], [rad], [rad/s], harmonic number

            print 'ece_int Iece = {:g}'.format(ece_int)
            print 'Rm = {:g}'.format(Rm)
            print 'zm = {:g}'.format(zm)

            # channel response in optics and IF
            dS = np.exp(-2*(dz[i]/Lcz)**4) * np.exp(-2*( (fsub[j]-np.mean(fsub))/Bcf )**4)
            S = S + dS

            int_meas[cn] = int_meas[cn] + ece_int * dS
            Rch[cn] = Rch[cn] + Rm * dS
            zch[cn] = zch[cn] + zm * dS

    # average over response
    int_meas[cn] = int_meas[cn] / S
    Rch[cn] = Rch[cn] / S
    zch[cn] = zch[cn] / S

    # radiation temperature
    rad_temp[cn] = int_meas[cn] / (np.mean(fsub)*2.0*np.pi*1e9/(2.0*np.pi*c))**2.0 / (1000.0*e) # [keV]
    abs_temp[cn] = F_Te(Rch[cn], zch[cn]) / (1000.0*e) # [keV]

    print 'S = {:g}'.format(S)
    print 'Rch = {:g}'.format(Rch[cn])
    print 'zch = {:g}'.format(zch[cn])
    print 'imeas = {:g}'.format(int_meas[cn])
    print 'rad_temp = {:g}'.format(rad_temp[cn])
    print 'abs_temp = {:g}'.format(abs_temp[cn])

#plt.plot(s,ams)
#plt.show()

#ece_int = integrate.simps(jms,x=s)
#print ece_int

#print dz
#print fsub
#print zsub
#print asub

#print Rp
#print zp

#plt.plot(Rp,zp)
#plt.show()
