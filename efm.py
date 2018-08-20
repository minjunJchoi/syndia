# ece intensity forward modeling
import subprocess
import h5py
import numpy as np
import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt

from kei import *
from pfunc import *
from bpath import torbeam_prof, write_inbeam, run_torbeam, set_beam_path, vac_beam_path
from pintp import intp_prof
from eceint import ece_intensity

from time import strftime

# def radiation_temperature(self, shot, dev, clist):

e = 1.602*1e-19
me = 9.109*1e-31
eps = 8.854*1e-12
c = 299792458
mc2 = me*c**2

ECEI_data_path = "/eceidata/exp_2015/"

# geqdsk_fn = "data/g013728.003900"
# ne_fn = "data/ne_3900.dat"
# Te_fn = "data/Te_3900.dat"
class EceFwdMod(object):
    def __init__(self):
        self.Dlist = []

    def get_profile(self, geqdsk_fn, Te_fn, ne_fn):
        self.geqdsk_fn = geqdsk_fn
        self.pf = ProfFunc(geqdsk_fn, Te_fn, ne_fn)

    def set_channel(self, shot, clist):
        self.shot = shot
        self.clist = expand_clist(clist)

        ## If ECEI
        self.Lcz = 9 # e^2 fallding practical vertical width for minilens [mm]
        self.Bcf = 0.3 # e^2 fallding for IF response [GHz]
        self.ecei = KstarEceiInfo(shot, clist)

    def run(self, fstart=-0.35, fend=0.35, Nf=10, zstart=-14, zend=14, Nz=10, ToR=0):
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
            torbeam_prof(self.geqdsk_fn, self.pf)

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

        # sub ray dz
        dz = np.linspace(zstart, zend, Nz) # dz [mm] of sub z rays at minilens

        for cn in range(0, cnum):
            ## If ECEI
            # channel numbers
            vn = int(self.clist[cn][(self.ecei.cnidx1):(self.ecei.cnidx1+2)])
            fn = int(self.clist[cn][(self.ecei.cnidx1+2):(self.ecei.cnidx1+4)])
            # define sub rays
            fsub = np.linspace((fn-1)*0.9 + 2.6 + self.ecei.lo + fstart, (fn-1)*0.9 + 2.6 + self.ecei.lo + fend, Nf) # frequency [GHz] of sub rays
            zsub = np.zeros(dz.size)
            asub = np.zeros(dz.size)
            S = 0
            ## loop over sub rays of a single channel
            for i in range(dz.size):
                # for ECEI; vacuum approximation until Rinit
                zsub[i], asub[i] = vac_beam_path(self.ecei, Rinit, vn, dz[i]) # vertical position [m] and rangle [rad] of sub rays at Rinit

                # for ECE on midplane
                # zsub[i] = dz[i]
                # asub[i] = 0

                for j in range(fsub.size):
                    # find beam path
                    if ToR == 0:
                        Rp, zp = run_torbeam(self.ecei.hn, fsub[j], asub[i], zsub[i], Rinit)

                    # take proper range of beam path
                    Rp, zp, theta = set_beam_path(Rp, zp, self.ecei.hn, fsub[j], pstart, pend, pint, self.pf)

                    # profile function along path
                    s, F_Bs, F_Tes, F_nes = intp_prof(Rp, zp, theta, self.pf, 0)

                    # calculate ECE intensity along path
                    ece_int, Rm, zm, thm, s, jms, ams, tau = ece_intensity(s, Rp, zp, theta, 2*np.pi*fsub[j]*1e9, self.ecei.hn, F_Bs, F_Tes, F_nes) # [m,m,m,rad,rad/s,hn,funcs]

                    print 'ece_int Iece = {:g}'.format(ece_int)
                    print 'Rm = {:g}'.format(Rm)
                    print 'zm = {:g}'.format(zm)

                    # channel response in optics and IF
                    dS = np.exp(-2*(dz[i]/self.Lcz)**4) * np.exp(-2*( (fsub[j]-np.mean(fsub))/self.Bcf )**4)
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
            abs_temp[cn] = self.pf.F_Te(Rch[cn], zch[cn]) / (1000.0*e) # [keV]

            print 'S = {:g}'.format(S)
            print 'Rch = {:g}'.format(Rch[cn])
            print 'zch = {:g}'.format(zch[cn])
            print 'imeas = {:g}'.format(int_meas[cn])
            print 'rad_temp = {:g}'.format(rad_temp[cn])
            print 'abs_temp = {:g}'.format(abs_temp[cn])

        return Rch, zch, int_meas, rad_temp, abs_temp

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
        elif 'ECEI' in clist[c] and len(clist[c]) == 16: # ECEI_GT0101 since 2018
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
