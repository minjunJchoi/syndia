# ece intensity forward modeling
import h5py
import math
import matplotlib.pyplot as plt

from vbp import vac_beam_path
from bpath import beam_path
# def radiation_temperature(self, shot, dev, clist):

data_path = '/eceidata/exp_2015/'
shot = 13728
dev = 'G'
clist = 'ECEI_G0101'

e = 1.602*1e-19
me = 9.109*1e-31
eps = 8.854*1e-12
c = 299792458
mc2 = me*c**2

Nz = 1 # number of vertical rays for a single channel
zstart = 0 # [mm] first ray at the mini lens -14
zend = 0 # [mm] last ray at the mini lens 14
Nf = 1 # number of frequency rays for a single channel
fstart = 0 # [GHz] frequency bandwidth -0.3
fend = 0 # [GHz] frequency bandwidth +0.3

pstart = 7.5 # [GHz] cal start point (hfs) = cold resonance + pstart
pend = -2 # [GHz] cal end point (lfs) = cold resonance + pend
pint = 0.1 # ECE integration path inter step. 0.1 = 10%

Rinit = 2.39 # where vacuum region ends and plasma region starts


## shot information
if shot < 19392:
    fname = "{:s}{:06d}/ECEI.{:06d}.{:s}FS.h5".format(data_path, shot, shot, dev)
    cnidx1 = 6
else:
    fname = "{:s}{:06d}/ECEI.{:06d}.{:s}.h5".format(data_path, shot, shot, dev)
    cnidx1 = 7

with h5py.File(fname, 'r') as f:
    # get attributes
    dset = f['ECEI']
    mode = dset.attrs['Mode']
    if 'O' in mode:
        hn = 1  # harmonic number
    elif 'X' in mode:
        hn = 2
    lo = dset.attrs['LoFreq']

###################### run tbprof; ready for TORBEAM

## loop for channels
cnum = len(clist)
intmeas = np.zeros(cnum) # ECE intensity measured
radtemp = np.zeros(cnum) # radiation temperature [keV] of each channel
Rch = np.zeros(cnum)  # R [m] of each channel
zch = np.zeros(cnum)  # z [m] of each channel
ach = np.zeros(cnum)  # angle [rad] of each channel
dz = np.linspace(zstart, zend, Nz) # dz [mm] of sub z rays at minilens
for c in range(0, cnum):
    vn = int(clist[c][(cnidx1):(cnidx1+2)])
    fn = int(clist[c][(cnidx1+2):(cnidx1+4)])

    # define sub rays
    fsub = np.linspace((fn-1)*0.9 + 2.6 + lo + fstart, (fn-1)*0.9 + 2.6 + lo + fend, Nf) # frequency [GHz] of sub rays
    zsub = np.zeros(dz.size)
    asub = np.zeros(dz.size)
    ## loop over sub rays of a single channel
    for i in range(dz.size):
        # vacuum approximation until Rinit
        zsub[i], asub[i] = vac_beam_path(shot, dev, Rinit, vn, dz[i]) # vertical position [m] and rangle [rad] of sub rays at Rinit

        for j in range(fsub.size):
            # find beam path
            Rp, zp, theta = beam_path(hn, fsub[j], asub[i], zsub[i], Rinit, pstart, pend, pint) # [GHz], [rad], [m], [m]

            # calculate ECE intensity along path


            #
            # intmeas[c]
            # radtemp[c]
            # Rch[c]
            # zch[c]

print dz
print fsub
print zsub
print asub

plt.plot(Rp, zp)
plt.show()
