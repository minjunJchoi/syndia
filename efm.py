"""
ECE intensity forward modeling; modified for both ECE and ECEI (2022.01.12)
Repo : https://github.com/minjunJchoi/syndia
Author : Minjun J. Choi (mjchoi@nfri.re.kr)
Collaborators : Jieun Lee, Yoonbum Nam
"""
import subprocess
import h5py
import numpy as np
import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import time

from kstar_ece_info import *
from kstar_ecei_info import *
from kstar_eceirt_info import *
from pfunc import *
from bpath import torbeam_prof, write_inbeam, run_torbeam, ray_tracing, set_beam_path, vac_beam_path
from pintp import intp_prof
from eceint import ece_intensity


e = 1.602*1e-19
me = 9.109*1e-31
eps = 8.854*1e-12
c = 299792458
mc2 = me*c**2

class EceFwdMod(object):
    def __init__(self):
        self.Dlist = []

    def set_profile(self, geqdsk_fn, Te_fn, ne_fn, bfactor=1.0):
        self.geqdsk_fn = geqdsk_fn
        self.pf = ProfFunc(geqdsk_fn, Te_fn, ne_fn, bfactor=bfactor)

    def set_channel(self, shot, clist):
        self.shot = shot

        #### diagnostics specific parameters 
        if 'ECEI' in clist[0] and shot < 35000:
            self.Lcz = 9 # e^2 fallding practical vertical width for minilens [mm]
            self.Bcf = 0.3 # e^2 fallding for IF response [GHz]
            self.diag = KstarEceiInfo(shot, clist)
        elif 'ECEI' in clist[0] and shot > 35000:
            self.Lcz = 9 # e^2 fallding practical vertical width for minilens [mm]
            self.Bcf = 0.3 # e^2 fallding for IF response [GHz]
            self.diag = KstarEceiRemoteInfo(shot, clist)            
        else:
            self.Lcz = 9 # e^2 fallding practical vertical width for minilens [mm]
            self.Bcf = 0.3 # e^2 fallding for IF response [GHz]
            self.diag = KstarEceInfo(shot, clist)

    def run(self, fstart=-0.35, fend=0.35, Nf=10, zstart=-14, zend=14, Nz=10, pstart=7.8, pend=-2, pint=0.1, Rinit=2.39, torbeam=0, select='mean'):
        ## bpath interp TORBEAM or Ray tracing
        ## pintp
        ## eceint
        # Nf = 1 # number of frequency rays for a single channel
        # fstart = 0 # [GHz] frequency bandwidth -0.3
        # fend = 0 # [GHz] frequency bandwidth +0.3
        # Nz = 1 # number of vertical rays for a single channel
        # zstart = 0 # [mm] first ray at the mini lens -14
        # zend = 0 # [mm] last ray at the mini lens 14

        # pstart = 7.8 # [GHz] cal start point (hfs) = cold resonance + pstart
        # pend = -2 # [GHz] cal end point (lfs) = cold resonance + pend
        # pint = 0.1 # ECE integration path inter step. 0.1 = 10%
        # Rinit = 2.39 # where vacuum region ends and plasma region starts

        # torbeam = 1 # 1 : TORBEAM, 0 : Ray tracing

        # ready for TORBEAM
        if torbeam == 1:
            torbeam_prof(self.geqdsk_fn, self.pf)

        ## loop for channels
        cnum = len(self.diag.clist)

        int_meas = np.zeros(cnum) # ECE intensity measured;
        # you use this to make a 'synthetic' ECE image from simulation data to compare with the measured ECE image
        rad_temp = np.zeros(cnum) # radiation temperature [keV] of each channel;
        # it is same with real abs_temp for black body;
        # you determine correctness of temperature/density profile (from Thomson, etc) (variable inputs of syndia) by comparing this output with what 'well calibrated' ECE gives you
        abs_temp = np.zeros(cnum) # abs_temp from F_Te (input of syndia)
        tau = np.zeros(cnum) # optical depth

        Rch = np.zeros(cnum)  # R [m] of each channel
        zch = np.zeros(cnum)  # z [m] of each channel
        ach = np.zeros(cnum)  # angle [rad] of each channel

        # sub ray dz
        dz = np.linspace(zstart, zend, Nz) # dz [mm] of sub z rays at minilens

        for cn in range(cnum):
            # define sub rays
            fsub = np.linspace(self.diag.ece_freq[cn] + fstart, self.diag.ece_freq[cn] + fend, Nf) # frequency [GHz] of sub rays
            zsub = np.zeros(dz.size)
            asub = np.zeros(dz.size)
            S = 0
            ## loop over sub rays of a single channel
            for i in range(dz.size):
                if 'ECEI' in self.diag.clist[0]: ## If ECEI; vacuum approximation until Rinit
                    zsub[i], asub[i] = vac_beam_path(self.diag, Rinit, self.diag.vn_list[cn], dz[i]) # vertical position [m] and rangle [rad] of sub rays at Rinit
                else: # for ECE on midplane
                    zsub[i] = dz[i]
                    asub[i] = 0

                for j in range(fsub.size):
                    # find beam path
                    if torbeam == 1:
                        Rp, zp = run_torbeam(self.diag.hn, fsub[j], asub[i], zsub[i], Rinit)
                        #plt.plot(Rp, zp, 'go-')
                        #Rprt, zprt = ray_tracing(self.diag.hn, fsub[j], asub[i], zsub[i], Rinit, self.pf)
                        #plt.plot(Rprt,zprt, 'bx-')
                        #plt.show()
                    else:
                        Rp, zp = ray_tracing(self.diag.hn, fsub[j], asub[i], zsub[i], Rinit, self.pf)

                    # find proper range of beam path and make interpolated beam path
                    Rp, zp, theta = set_beam_path(Rp, zp, self.diag.hn, fsub[j], pstart, pend, pint, self.pf)

                    # profile function along path
                    s, F_Bs, F_Tes, F_nes = intp_prof(Rp, zp, theta, self.pf, 0)

                    # calculate ECE intensity along path
                    ece_int, Rm, zm, thm, s, jms, ams = ece_intensity(s, Rp, zp, theta, 2*np.pi*fsub[j]*1e9, self.diag.hn, F_Bs, F_Tes, F_nes, select=select) # [m,m,m,rad,rad/s,hn,funcs]

                    print('{:d}-ray tau = {:g}'.format(j, integrate.trapz(ams,x=s)))
                    print('{:d}-ray ece_int Iece = {:g}'.format(j, ece_int))
                    print('{:d}-ray Rm = {:g}'.format(j, Rm))
                    print('{:d}-ray zm = {:g}'.format(j, zm))
                    print('{:d}-ray thm = {:g}'.format(j, thm))

                    # channel response in optics and IF
                    dS = np.exp(-2*(dz[i]/self.Lcz)**4) * np.exp(-2*( (fsub[j]-np.mean(fsub))/self.Bcf )**4)
                    S = S + dS

                    int_meas[cn] = int_meas[cn] + ece_int * dS
                    tau[cn] = tau[cn] + integrate.trapz(ams,x=s) * dS                    
                    Rch[cn] = Rch[cn] + Rm * dS
                    zch[cn] = zch[cn] + zm * dS

            # average over response
            int_meas[cn] = int_meas[cn] / S
            tau[cn] = tau[cn] / S            
            Rch[cn] = Rch[cn] / S
            zch[cn] = zch[cn] / S

            # radiation temperature
            rad_temp[cn] = int_meas[cn] / (np.mean(fsub)*2.0*np.pi*1e9/(2.0*np.pi*c))**2.0 / (1000.0*e) # [keV]
            abs_temp[cn] = self.pf.F_Te(Rch[cn], zch[cn]) / (1000.0*e) # [keV]

            # print('S = {:g}'.format(S))
            print('Rch = {:g}'.format(Rch[cn]))
            print('zch = {:g}'.format(zch[cn]))
            print('imeas = {:g}'.format(int_meas[cn]))
            print('rad_temp = {:g}'.format(rad_temp[cn]))
            print('abs_temp = {:g}'.format(abs_temp[cn]))
            print('tau = {:g}'.format(tau[cn]))
            #print('abs_dens = {:g}'.format(self.pf.F_ne(Rch[cn], zch[cn])/1e19))

        return Rch, zch, int_meas, tau, rad_temp, abs_temp

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
