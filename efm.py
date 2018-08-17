# ece intensity forward modeling
import subprocess
import h5py
import numpy as np
import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt

from pfunc import *
from bpath import vac_beam_path, tb_beam_path
from eceint import ece_intensity
from pintp import intp_prof

from time import strftime

# def radiation_temperature(self, shot, dev, clist):

e = 1.602*1e-19
me = 9.109*1e-31
eps = 8.854*1e-12
c = 299792458
mc2 = me*c**2

ECEI_data_path = "/eceidata/exp_2015/"
TB_path_eq = "/home/users/mjchoi/torbeam_ifortran/eqdsk2topfile/"
TB_path_run = "/home/users/mjchoi/torbeam_ifortran/run_torbeam/"

# geqdsk_fn = "data/g013728.003900"
# ne_fn = "data/ne_3900.dat"
# Te_fn = "data/Te_3900.dat"
class EceFwdMod(object):
    def __init__(self):
        self.Dlist = []

    def get_profile(self, geqdsk_fn, Te_fn, ne_fn):
        self.geqdsk_fn = geqdsk_fn
        self.pf = ProfFunc(geqdsk_fn, Te_fn, ne_fn)

    def torbeam_prof(self):
        ## save for TORBEAM
        # read geqdsk file from the selected EFIT run time and convert it to topfile for TORBEAM
        args = "{}readeqdsk<{}".format(TB_path_eq, self.geqdsk_fn)
        re = subprocess.check_output(args, shell=True)
        args = "mv topfile {}".format(TB_path_run)
        re = subprocess.check_output(args, shell=True)

        # save ne.dat for TORBEAM
        x = np.sqrt(self.pf.psin_ne)
        y = self.pf.ne
        with open("{}ne.dat".format(TB_path_run), 'w') as f:
            f.write('{}\n'.format(x.size))
            for i in range(x.size):
                f.write('{} {}\n'.format(x[i], y[i]))

        # save Te.dat for TORBEAM
        x = np.sqrt(self.pf.psin_Te)
        y = self.pf.Te
        with open("{}Te.dat".format(TB_path_run), 'w') as f:
            f.write('{}\n'.format(x.size))
            for i in range(x.size):
                f.write('{} {}\n'.format(x[i], y[i]))

        print 'profile data for TORBEAM are saved in {}'.format(TB_path_run)

    def set_channel(self, shot, clist):
        self.shot = shot
        self.clist = expand_clist(clist)

        ## If ECEI
        self.Lcz = 9 # e^2 fallding practical vertical width for minilens [mm]
        self.Bcf = 0.3 # e^2 fallding for IF response [GHz]

        if shot < 19392:
            cnidx1 = 6
            dev = self.clist[0][5]
            fname = "{:s}{:06d}/ECEI.{:06d}.{:s}FS.h5".format(ECEI_data_path, shot, shot, dev)
        else:
            cnidx1 = 7
            dev = self.clist[0][5:7]
            fname = "{:s}{:06d}/ECEI.{:06d}.{:s}.h5".format(ECEI_data_path, shot, shot, dev)

        with h5py.File(fname, 'r') as f:
            # get attributes
            dset = f['ECEI']
            mode = dset.attrs['Mode'].strip()
            if mode is 'O':
                self.hn = 1  # harmonic number
            elif mode is 'X':
                self.hn = 2
            self.lo = dset.attrs['LoFreq']

    def rad_temp(self, fstart=-0.35, fend=0.35, Nf=10, zstart=-14, zend=14, Nz=10, ToR=0):
        ## bpath interp TORBEAM or Ray tracing
        ## pintp
        ## eceint
        # Nf = 1 # number of frequency rays for a single channel
        # fstart = 0 # [GHz] frequency bandwidth -0.3
        # fend = 0 # [GHz] frequency bandwidth +0.3
        # Nz = 1 # number of vertical rays for a single channel
        # zstart = 0 # [mm] first ray at the mini lens -14
        # zend = 0 # [mm] last ray at the mini lens 14
        # ToR = 0 # 0 : TORBEAM, 1 : Ray tracing

        pstart = 7.8 # [GHz] cal start point (hfs) = cold resonance + pstart
        pend = -2 # [GHz] cal end point (lfs) = cold resonance + pend
        pint = 0.1 # ECE integration path inter step. 0.1 = 10%
        Rinit = 2.39 # where vacuum region ends and plasma region starts

        # ready for TORBEAM
        if ToR == 0:
            self.torbeam_prof()

        return

        ## loop for channels
        cnum = len(self.clist)
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

                    print 'start time = {}'.format(strftime("%y%m%d-%H%M%S"))
                    # calculate ECE intensity along path with profile function along path
                    s, F_Bs, F_Tes, F_nes = intp_prof(Rp, zp, theta, 0)
                    ece_int, Rm, zm, thm, s, jms, ams, tau = ece_intensity(s, Rp, zp, theta, 2*np.pi*fsub[j]*1e9, hn, F_Bs, F_Tes, F_nes) # [m,m,m,rad,rad/s,hn,funcs]
                    print 'end time ={}'.format(strftime("%y%m%d-%H%M%S"))

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

def expand_clist(clist):
    # IN : List of channel names (e.g. 'ECEI_G1201-1208' or 'ECEI_GT1201-1208').
    # OUT : Expanded list (e.g. 'ECEI_G1201', ..., 'ECEI_G1208')

    # KSTAR ECEI
    exp_clist = []
    for c in range(len(clist)):
        if 'ECEI' in clist[c] and len(clist[c]) == 15: # before 2018
            vi = int(clist[c][6:8])
            fi = int(clist[c][8:10])
            vf = int(clist[c][11:13])
            ff = int(clist[c][13:15])

            for v in range(vi, vf+1):
                for f in range(fi, ff+1):
                    exp_clist.append(clist[c][0:6] + '%02d' % v + '%02d' % f)
        elif 'ECEI' in clist[c] and len(clist[c]) == 16: # since 2018
            vi = int(clist[c][7:9])
            fi = int(clist[c][9:11])
            vf = int(clist[c][12:14])
            ff = int(clist[c][14:16])

            for v in range(vi, vf+1):
                for f in range(fi, ff+1):
                    exp_clist.append(clist[c][0:6] + '%02d' % v + '%02d' % f)
        else:
            exp_clist.append(clist[c])
    clist = exp_clist

    return clist
