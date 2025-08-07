"""
Enhanced parallel processing functions for ECE forward modeling
Supporting both channel-level and sub-ray level parallelization
"""
import numpy as np
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import scipy.integrate as integrate
from bpath import vac_beam_path, run_torbeam, ray_tracing, set_beam_path
from pintp import intp_prof
from eceint import ece_intensity

def process_sub_ray(args):
    """Process a single sub-ray (i, j combination) for ultra-fine parallelization"""
    i, j, dz_i, fsub_j, diag, cn, Rinit, torbeam, pf, pstart, pend, pint, select, Lcz, Bcf, fsub_mean = args
    
    # Calculate zsub and asub for this i
    if 'ECEI' in diag.clist[0]:
        zsub_i, asub_i = vac_beam_path(diag, Rinit, diag.vn_list[cn], dz_i)
    else:
        zsub_i = dz_i
        asub_i = 0
    
    # Find beam path
    if torbeam == 1:
        Rp, zp = run_torbeam(diag.hn, fsub_j, asub_i, zsub_i, Rinit)
    else:
        Rp, zp = ray_tracing(diag.hn, fsub_j, asub_i, zsub_i, Rinit, pf)

    # Find proper range of beam path and make interpolated beam path
    Rp, zp, theta = set_beam_path(Rp, zp, diag.hn, fsub_j, pstart, pend, pint, pf)

    # Profile function along path
    s, F_Bs, F_Tes, F_nes = intp_prof(Rp, zp, theta, pf, 0)

    # Calculate ECE intensity along path
    ece_int, Rm, zm, thm, s, jms, ams = ece_intensity(
        s, Rp, zp, theta, 2*np.pi*fsub_j*1e9, diag.hn, 
        F_Bs, F_Tes, F_nes, select=select
    )
    
    # Channel response in optics and IF
    dS = np.exp(-2*(dz_i/Lcz)**4) * np.exp(-2*((fsub_j-fsub_mean)/Bcf)**4)
    
    tau_contrib = integrate.trapz(ams, x=s)
    
    return {
        'i': i, 'j': j,
        'ece_int': ece_int, 'Rm': Rm, 'zm': zm, 'thm': thm,
        'tau_contrib': tau_contrib, 'dS': dS
    }

def process_single_channel_fine_parallel(args):
    """
    Process a single channel with sub-ray level parallelization
    This can be useful when you have few channels but many sub-rays
    """
    cn, diag, pf, fstart, fend, Nf, zstart, zend, Nz, pstart, pend, pint, Rinit, torbeam, select, Lcz, Bcf, n_threads = args
    
    # sub ray dz
    dz = np.linspace(zstart, zend, Nz)
    
    # define sub rays
    fsub = np.linspace(diag.ece_freq[cn] + fstart, diag.ece_freq[cn] + fend, Nf)
    fsub_mean = np.mean(fsub)
    
    # Prepare all sub-ray combinations
    sub_ray_args = []
    for i in range(len(dz)):
        for j in range(len(fsub)):
            args_ij = (i, j, dz[i], fsub[j], diag, cn, Rinit, torbeam, pf, 
                      pstart, pend, pint, select, Lcz, Bcf, fsub_mean)
            sub_ray_args.append(args_ij)
    
    # Process sub-rays in parallel using threads (good for I/O bound tasks)
    results = []
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = [executor.submit(process_sub_ray, args) for args in sub_ray_args]
        for future in as_completed(futures):
            results.append(future.result())
    
    # Aggregate results
    S = sum(r['dS'] for r in results)
    int_meas_cn = sum(r['ece_int'] * r['dS'] for r in results) / S
    tau_cn = sum(r['tau_contrib'] * r['dS'] for r in results) / S
    Rch_cn = sum(r['Rm'] * r['dS'] for r in results) / S
    zch_cn = sum(r['zm'] * r['dS'] for r in results) / S
    
    # Radiation temperature
    rad_temp_cn = int_meas_cn / (fsub_mean*2.0*np.pi*1e9/(2.0*np.pi*299792458))**2.0 / (1000.0*1.602*1e-19)
    abs_temp_cn = pf.F_Te(Rch_cn, zch_cn) / (1000.0*1.602*1e-19)
    
    return cn, Rch_cn, zch_cn, int_meas_cn, tau_cn, rad_temp_cn, abs_temp_cn
