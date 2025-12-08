"""
ECE intensity forward modeling; modified for both ECE and ECEI (2022.01.12)
Repo : https://github.com/minjunJchoi/syndia
Author : Minjun J. Choi (mjchoi@kfe.re.kr)
Collaborators : Jieun Lee, Yoonbum Nam, Github Copilot (Claude Sonnet 4)
"""
import subprocess
import h5py
import numpy as np
import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

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

def process_single_channel(args):
    """Process a single channel for parallel computation"""
    cn, diag, pf, fstart, fend, Nf, zstart, zend, Nz, pstart, pend, pint, Rinit, torbeam, select, Lcz, Bcf = args
    
    # sub ray dz
    dz = np.linspace(zstart, zend, Nz) # dz [mm] of sub z rays at minilens
    
    # define sub rays
    fsub = np.linspace(diag.ece_freq[cn] + fstart, diag.ece_freq[cn] + fend, Nf) # frequency [GHz] of sub rays
    zsub = np.zeros(dz.size)
    asub = np.zeros(dz.size)
    S = 0
    I = 0
    
    int_meas = 0
    tau = 0
    Rch = 0
    zch = 0
    
    ## loop over sub rays of a single channel
    for i in range(dz.size):
        if 'ECEI' in diag.clist[0]: ## If ECEI; vacuum approximation until Rinit
            zsub[i], asub[i] = vac_beam_path(diag, Rinit, diag.vn_list[cn], dz[i]) # vertical position [m] and rangle [rad] of sub rays at Rinit
        else: # for ECE on midplane
            zsub[i] = dz[i]
            asub[i] = 0

        for j in range(fsub.size):
            # find beam path
            if torbeam == 1:
                Rp, zp = run_torbeam(diag.hn, fsub[j], asub[i], zsub[i], Rinit)
            else:
                Rp, zp = ray_tracing(diag.hn, fsub[j], asub[i], zsub[i], Rinit, pf)

            # find proper range of beam path and make interpolated beam path
            Rp, zp, theta = set_beam_path(Rp, zp, diag.hn, fsub[j], pstart, pend, pint, pf)

            # profile function along path
            s, F_Bs, F_Tes, F_nes = intp_prof(Rp, zp, theta, pf, 0)

            # calculate ECE intensity along path
            ece_int, Rm, zm, thm, s, jms, ams = ece_intensity(s, Rp, zp, theta, 2*np.pi*fsub[j]*1e9, diag.hn, F_Bs, F_Tes, F_nes, select=select)

            # channel response in optics and IF
            dS = np.exp(-2*(dz[i]/Lcz)**4) * np.exp(-2*( (fsub[j]-np.mean(fsub))/Bcf )**4)
            S = S + dS
            I = I + ece_int * dS

            int_meas = int_meas + ece_int * dS
            tau = tau + integrate.trapz(ams,x=s) * dS                    
            Rch = Rch + Rm * ece_int * dS
            zch = zch + zm * ece_int * dS

    # average over response
    int_meas = int_meas / S
    tau = tau / S            
    Rch = Rch / I
    zch = zch / I

    # radiation temperature
    rad_temp = int_meas / (np.mean(fsub)*2.0*np.pi*1e9/(2.0*np.pi*c))**2.0 / (1000.0*e) # [keV]
    abs_temp = pf.F_Te(Rch, zch) / (1000.0*e) # [keV]

    return cn, Rch, zch, int_meas, tau, rad_temp, abs_temp

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

    def run(self, fstart=-0.35, fend=0.35, Nf=10, zstart=-14, zend=14, Nz=10, pstart=7.8, pend=-2, pint=0.1, Rinit=2.39, torbeam=0, select='mean', parallel=True, n_processes=None):
        """
        Run ECE forward modeling
        
        Parameters:
        -----------
        Nf = 1 # number of frequency rays for a single channel
        fstart = 0 # [GHz] frequency bandwidth -0.3
        fend = 0 # [GHz] frequency bandwidth +0.3
        Nz = 1 # number of vertical rays for a single channel
        zstart = 0 # [mm] first ray at the mini lens -14
        zend = 0 # [mm] last ray at the mini lens 14
        pstart = 7.8 # [GHz] cal start point (hfs) = cold resonance + pstart
        pend = -2 # [GHz] cal end point (lfs) = cold resonance + pend
        pint = 0.1 # ECE integration path inter step. 0.1 = 10%
        Rinit = 2.39 # where vacuum region ends and plasma region starts
        torbeam = 1 # 1 : TORBEAM, 0 : Ray tracing

        parallel = bool, default=True
            If True, use parallel processing; Parallelize over channels using multiprocessing
        n_processes = int, default=None
            Number of processes to use. If None, use number of CPU cores
        """

        # ready for TORBEAM
        if torbeam == 1:
            torbeam_prof(self.geqdsk_fn, self.pf)

        ## loop for channels
        cnum = len(self.diag.clist)

        int_meas = np.zeros(cnum) # ECE intensity measured;
        rad_temp = np.zeros(cnum) # radiation temperature [keV] of each channel;
        abs_temp = np.zeros(cnum) # abs_temp from F_Te (input of syndia)
        tau = np.zeros(cnum) # optical depth

        Rch = np.zeros(cnum)  # R [m] of each channel
        zch = np.zeros(cnum)  # z [m] of each channel
        ach = np.zeros(cnum)  # angle [rad] of each channel

        if parallel:
            # Parallel processing
            if n_processes is None:
                n_processes = mp.cpu_count()
            
            print(f"Using channel-level parallel processing with {n_processes} processes for {cnum} channels")
            
            # Prepare arguments for each channel
            args_list = []
            for cn in range(cnum):
                args = (cn, self.diag, self.pf, fstart, fend, Nf, zstart, zend, Nz, 
                        pstart, pend, pint, Rinit, torbeam, select, self.Lcz, self.Bcf)
                args_list.append(args)
            
            # Process channels in parallel
            with ProcessPoolExecutor(max_workers=n_processes) as executor:
                futures = [executor.submit(process_single_channel, args) for args in args_list]
                
                for future in as_completed(futures):
                    cn, Rch_cn, zch_cn, int_meas_cn, tau_cn, rad_temp_cn, abs_temp_cn = future.result()
                    
                    Rch[cn] = Rch_cn
                    zch[cn] = zch_cn
                    int_meas[cn] = int_meas_cn
                    tau[cn] = tau_cn
                    rad_temp[cn] = rad_temp_cn
                    abs_temp[cn] = abs_temp_cn
                    
                    print(f'Channel {cn}: Rch = {Rch[cn]:.6g}, zch = {zch[cn]:.6g}')
                    print(f'Channel {cn}: imeas = {int_meas[cn]:.6g}, rad_temp = {rad_temp[cn]:.6g}, abs_temp = {abs_temp[cn]:.6g}, tau = {tau[cn]:.6g}')
        else:
            # Sequential processing (original code)
            self._run_sequential(cnum, fstart, fend, Nf, zstart, zend, Nz, pstart, pend, pint, Rinit, torbeam, select,
                               int_meas, rad_temp, abs_temp, tau, Rch, zch)

        return Rch, zch, int_meas, tau, rad_temp, abs_temp

    def _run_sequential(self, cnum, fstart, fend, Nf, zstart, zend, Nz, pstart, pend, pint, Rinit, torbeam, select,
                       int_meas, rad_temp, abs_temp, tau, Rch, zch):
        """Sequential version using process_single_channel function"""
        print(f"Using sequential processing for {cnum} channels")
        
        for cn in range(cnum):
            # Prepare arguments for the channel
            args = (cn, self.diag, self.pf, fstart, fend, Nf, zstart, zend, Nz, 
                    pstart, pend, pint, Rinit, torbeam, select, self.Lcz, self.Bcf)
            
            # Process single channel
            cn_result, Rch_cn, zch_cn, int_meas_cn, tau_cn, rad_temp_cn, abs_temp_cn = process_single_channel(args)
            
            # Store results
            Rch[cn] = Rch_cn
            zch[cn] = zch_cn
            int_meas[cn] = int_meas_cn
            tau[cn] = tau_cn
            rad_temp[cn] = rad_temp_cn
            abs_temp[cn] = abs_temp_cn
            
            print(f'Channel {cn}: Rch = {Rch[cn]:.6g}, zch = {zch[cn]:.6g}')
            print(f'Channel {cn}: imeas = {int_meas[cn]:.6g}, rad_temp = {rad_temp[cn]:.6g}, abs_temp = {abs_temp[cn]:.6g}, tau = {tau[cn]:.6g}')

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
