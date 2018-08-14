# otbain 2D profiles
from geqdsk_dk import geqdsk_dk as geqdsk_dk
import numpy as np
import subprocess

geqdsk_fn = "data/g013728.003900"
ne_fn = "data/ne_3900.dat"
Te_fn = "data/Te_3900.dat"
TB_path_eq = "/home/users/mjchoi/torbeam_ifortran/eqdsk2topfile/"
TB_path_run = "/home/users/mjchoi/torbeam_ifortran/run_torbeam/"
## read files
# geqdsk file
geq = geqdsk_dk(filename=geqdsk_fn)
# density data: psin_ne, ne [10^-19 m^-3]; well behaving data
with open(ne_fn, 'r') as f:
    psin_ne, ne = np.loadtxt(f, unpack=True)
# temperature data: psin_Te, Te [keV]; well behaving data
with open(Te_fn, 'r') as f:
    psin_Te, Te = np.loadtxt(f, unpack=True)

## save for TORBEAM
# read geqdsk file from the selected EFIT run time and convert it to topfile for TORBEAM
args = "{}readeqdsk<{}".format(TB_path_eq, gfn)
re = subprocess.check_output(args, shell=True)
args = "mv topfile {}".format(TB_path_run)
re = subprocess.check_output(args, shell=True)

# save ne.dat for TORBEAM
x = np.sqrt(psin_ne)
y = ne
with open("{}ne.dat".format(TB_path_run), 'w') as f:
    f.write('{}\n'.format(x.size))
    for i in range(x.size):
        f.write('{} {}\n'.format(x[i], y[i]))

# save Te.dat for TORBEAM
x = np.sqrt(psin_Te)
y = Te
with open("{}Te.dat".format(TB_path_run), 'w') as f:
    f.write('{}\n'.format(x.size))
    for i in range(x.size):
        f.write('{} {}\n'.format(x[i], y[i]))
