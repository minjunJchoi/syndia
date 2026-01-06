#!/usr/bin/env python3
"""
Example script 
"""

import pickle
import time
import numpy as np
from efm import EceFwdMod
import matplotlib.pyplot as plt

# Shot info
shot = 13728
mean_time = 3900
bfactor = 1.01 # to consider B field uncertainty; B=B*bfactor
geqdsk_fn = 'example_data/g013728.003900'
Te_fn = 'example_data/Te_3900.dat' # normalized psi(:), Te(:) [keV]
ne_fn = 'example_data/ne_3900.dat' # normalized psi(:), ne(:) [1e19 m^-3]

# Channel info
clist = ['ECEI_G0101-2408']

# Parallel run with multiprocessing or not
parallel = True
n_processes = 16  # Number of processes for parallel execution

# Input Parameters for single central ray per each channel (position)
input_params = {
    'fstart': 0, 'fend': 0, 'Nf': 1, 'zstart': 0, 'zend': 0, 'Nz': 1,
    'pstart': 7.8, 'pend': -2, 'pint': 0.1, 'Rinit': 2.39, 'torbeam': 0, 'select': 'max'
}
# # Input parameters for multiple rays per each channel (synthetic data)
# input_params = {
#     'fstart': -0.35, 'fend': 0.35, 'Nf': 10, 'zstart': -14, 'zend': 14, 'Nz': 10, 
#     'pstart': 7.8, 'pend': -2, 'pint': 0.1, 'Rinit': 2.39, 'torbeam': 0, 'select': 'max'
# }


def run_equilibrium(fname=None):
    """No perturbations; ne0, Te0, I_E0"""   
   
    # Initialize the forward model
    efm = EceFwdMod()
    
    # Set up profile without perturbations
    efm.set_profile(geqdsk_fn, Te_fn, ne_fn, bfactor=bfactor)  

    # Set up channels 
    efm.set_channel(shot, clist)

    print("="*60)
    print("Run ECE Forward Modeling: Equilibrium without perturbations")
    print("="*60)
    
    # Parallel run
    start_time = time.time()
    results_par = efm.run(parallel=parallel, n_processes=n_processes, **input_params)
    par_time = time.time() - start_time
    print(f"Execution time: {par_time:.2f} seconds")

    # Save results
    if fname is not None:
        Rch, zch, _, tau, rad_temp, abs_temp = results_par
        with open(fname, 'wb') as fout:
            pickle.dump([clist, Rch, zch, tau, rad_temp, abs_temp], fout)
        print(f"Results saved to {fname}")


def verify_results():
    print("="*60)
    print("Verifying results consistency...")
    print("="*60)

    select_param = input_params['select']
    fname = f'example_data/ecei_pos_{shot}_{clist[0]}_{mean_time}ms_b{bfactor}_{select_param}_EQ_PRL0.pkl'
    with open(fname, 'rb') as fin:
        [_, Rch, zch, tau, rad_temp, abs_temp] = pickle.load(fin)
        print(f"Loaded results from {fname}")

    fname = f'example_data/ecei_pos_{shot}_{clist[0]}_{mean_time}ms_b{bfactor}_{select_param}_EQ_PRL.pkl'
    with open(fname, 'rb') as fin:
        [_, Rch_CMP, zch_CMP, tau_CMP, rad_temp_CMP, abs_temp_CMP] = pickle.load(fin)
        print(f"Loaded results from {fname}")

    fig, (a1, a2, a3, a4, a5) = plt.subplots(5, 1, figsize=(6, 6))
    a1.plot(Rch, '-x')
    a1.plot(Rch_CMP, '--')
    a1.set_title(f'Rch: max diff {np.max(np.abs(Rch - Rch_CMP)):.2e}')
    a2.plot(zch, '-x')
    a2.plot(zch_CMP, '--')
    a2.set_title(f'zch: max diff {np.max(np.abs(zch - zch_CMP)):.2e}')
    a3.plot(tau, '-x')
    a3.plot(tau_CMP, '--')
    a3.set_title(f'tau: max diff {np.max(np.abs(tau - tau_CMP)):.2e}')
    a4.plot(rad_temp, '-x')
    a4.plot(rad_temp_CMP, '--')
    a4.set_title(f'rad_temp: max diff {np.max(np.abs(rad_temp - rad_temp_CMP)):.2e}')
    a5.plot(abs_temp, '-x')
    a5.plot(abs_temp_CMP, '--')
    a5.set_title(f'abs_temp: max diff {np.max(np.abs(abs_temp - abs_temp_CMP)):.2e}')
    plt.tight_layout()
    plt.show()


def check_channels(fname=None):
    with open(fname, 'rb') as fin:
        [_, Rch, zch, tau, rad_temp, abs_temp] = pickle.load(fin)    

    # Create a 2x2 subplot grid and flatten axes for easy unpacking
    fig, axes = plt.subplots(2, 2, figsize=(6, 10), sharex=True, sharey=True)
    ax1, ax2, ax3, ax4 = axes.ravel()

    sc1 = ax1.scatter(Rch, zch, 100, rad_temp, marker='s')
    fig.colorbar(sc1, ax=ax1)
    ax1.set_title('Radiation Temperature')
    ax1.set_xlabel('R (m)')
    ax1.set_ylabel('Z (m)')
    ax1.set_aspect('equal', adjustable='box')

    sc2 = ax2.scatter(Rch, zch, 100, abs_temp, marker='s')
    fig.colorbar(sc2, ax=ax2)
    ax2.set_xlabel('R (m)')
    ax2.set_ylabel('Z (m)')
    ax2.set_aspect('equal', adjustable='box')
    ax2.set_title('Absolute Temperature')
    
    pdata = 1 - np.abs(abs_temp - rad_temp)/abs_temp
    sc3 = ax3.scatter(Rch, zch, 100, pdata, marker='s', vmin=0, vmax=1)
    fig.colorbar(sc3, ax=ax3)
    ax3.set_title('Validity')
    ax3.set_xlabel('R [m]')
    ax3.set_ylabel('z [m]')
    ax3.set_aspect('equal', adjustable='box')

    pdata = tau
    sc4 = ax4.scatter(Rch, zch, 100, pdata, marker='s', vmin=0, vmax=3)
    fig.colorbar(sc4, ax=ax4)
    ax4.set_title('Optical depth')
    ax4.set_xlabel('R [m]')
    ax4.set_ylabel('z [m]')
    ax4.set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # run equilibrium case
    fname = f'example_data/ecei_pos_{shot}_{clist[0]}_{mean_time}ms_b{bfactor}_{input_params["select"]}.pkl'
    run_equilibrium(fname)

    # # compare different runs
    # verify_results()

    # check channels 
    fname = f'example_data/ecei_pos_{shot}_{clist[0]}_{mean_time}ms_b{bfactor}_{input_params["select"]}.pkl'
    check_channels(fname)


