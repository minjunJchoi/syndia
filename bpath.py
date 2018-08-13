import subprocess
from scipy import interpolate

from inbeam import write_inbeam
from pfunc import *

TB_path = "/home/users/mjchoi/torbeam_ifortran/"
TB_path_run = "/home/users/mjchoi/torbeam_ifortran/run_torbeam/"

## functions needed
wce = lambda R,z: e*F_B(R,z)/me # [rad/s]

def beam_path(hn, freq, ainit, zinit, Rinit, pstart, pend, pint):
    ## TORBEAM
    # initial parameters
    write_inbeam(hn, freq*1e9, ainit/np.pi*180, zinit*100, Rinit*100)
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
        if np.abs(freq + pend - fRz) < 0.3:
            idx1 = j # no need to be very accurate
        if np.abs(freq - fRz) < 0.3:
            Rcold = Rp[i] # EC resonance position [cm] no need to be accurate
        if np.abs(freq + pstart - fRz) < 0.3:
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


beam_path(hn, freq, ainit, zinit, Rinit, pstart, pend, pint)
