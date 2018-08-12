# ece intensity forward modeling

from vbp import vac_beam_path
from inbeam import write_inbeam
from pfunc import *
import subprocess
import math
# def radiation_temperature(self, shot, dev, clist):

TB_path = "/home/users/mjchoi/torbeam_ifortran/"

shot = 13728
dev = 'G'
clist = 'ECEI_G1208'

e = 1.602*1e-19
me = 9.109*1e-31
eps = 8.854*1e-12
c = 299792458
mc2 = me*c**2

Nz = 5 # number of vertical rays for a single channel
zstart = -10 # [mm] first ray at the mini lens -14
zend = 10 # [mm] last ray at the mini lens 14
Nf = 5 # number of frequency rays for a single channel
fstart = -0.2 # [GHz] frequency bandwidth -0.3
fend = 0.2 # [GHz] frequency bandwidth +0.3

pstart = 7.5 # [GHz] cal start point (hfs) = cold resonance + pstart
pend = -2 # [GHz] cal end point (lfs) = cold resonance + pend
pint = 0.1 # ECE integration path inter step. 0.1 = 10%

R_vac_end = 2.39 # where vacuum region ends and plasma region starts

## functions needed
wce = lambda R,z: e*F_B(R,z)/me # [rad/s]

## shot information
if shot < 19392:
    fname = "{:s}{:06d}/ECEI.{:06d}.{:s}FS.h5".format(DIR, shot, shot, dev)
    cnidx1 = 6
else:
    fname = "{:s}{:06d}/ECEI.{:06d}.{:s}.h5".format(DIR, shot, shot, dev)
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
    fs = np.linspace((fn-1)*0.9 + 2.6 + lo + fstart, (fn-1)*0.9 + 2.6 + lo + fend, Nf) # frequency [GHz] of sub rays
    zs = np.zeros(dz.size)
    as = np.zeros(dz.size)
    ## loop over sub rays of a single channel
    for i in range(dz.size):
        # vacuum approximation until R_vac_end
        zs[i], as[i] = vac_beam_path(shot, dev, R_vac_end, vn, dz[i]) # vertical position [m] and rangle [rad] of sub rays at R_vac_end

        for j in range(fs.size):
            ## TORBEAM
            # initial parameters
            write_inbeam(hn, fs[j]*1e9, as[i]/np.pi*180, zs[i]*100, R_vac_end*100)
            # run
            args = "{}run.sh".format(TB_path)
            re = subprocess.check_output(args, shell=True)
            # obtain output
            with open("{}t1_LIB.dat".format(TB_path_run), 'r') as f:
                Rp, zp = np.loadtxt(f, usecols=(0,1), unpack=True)
                Rp = Rp/100 # R [m] beam path, from R_vec_end to core; decreasing
                zp = zp/100 # z [m] beam path
                # it can be curved and get double values. Find the first resonance position

            # find the proper range
            for i in range(Rp.size): # increasing during loop
                fRz = wce(Rp[i], zp[i])/(2*np.pi*1e9)*hn # EC frequency [GHz]
                if np.abs(fs[j] + pend - fRz) < 0.3:
                    idx1 = j # no need to be very accurate
                if np.abs(fs[j] - fRz) < 0.3:
                    Rcold = Rp[i] # EC resonance position [cm] no need to be accurate
                if np.abs(fs[j] + pstart - fRz) < 0.3:
                    idx2 = j # no need to be very accurate
                    break

            # Rp, zp between idx1 idx2; calculate angle between emission direction and B-field
            Rp = Rp[idx1:(idx2+1)]
            zp = zp[idx1:(idx2+1)]
            theta = np.zeros(Rp.size)
            for i in range(1,Rp.size):
                Rvec = [-(Rp[i]-Rp[i-1]), -(zp[i]-zp[i-1]), 0] # opposite direction for emission path
                Bvec = F_Bvec(Rp[i], zp[i])
                theta[i] = math.acos( Bvec.dot(Rvec) / ( np.sqrt(Bvec.dot(Bvec)) * np.sqrt(Rvec.dot(Rvec)) ) ) # [rad]

            # interpolation (for better accuracy) and change direction from hfs to lfs
            idx = np.range(idx2, idx1-1, -1)
            nidx = np.range(idx2, idx1, -pint)
            fRp = interpolate.interp1d(idx, Rp[idx], kind='linear')
            fzp = interpolate.interp1d(idx, zp[idx], kind='linear')
            fth = interpolate.interp1d(idx, theta[idx], kind='linear')
            Rp = fRp(nidx)
            zp = fzp(nidx)
            theta = fth(nidx)

            # calculate ECE intensity along path


            #
            intmeas[c]
            radtemp[c]
            Rch[c]
            zch[c]
