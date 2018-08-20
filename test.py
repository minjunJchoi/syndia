import pickle
import scipy.io as sio

from efm import *

geqdsk_fn = 'data/g013728.003900'
Te_fn = 'data/Te_3900.dat'
ne_fn = 'data/ne_3900.dat'

A = EceFwdMod()

A.get_profile(geqdsk_fn, Te_fn, ne_fn)

A.set_channel(13728,['ECEI_G0101-2408'])

Rch, zch, _, _, abs_temp = A.run(fstart=0,fend=0,Nf=1,zstart=0,zend=0,Nz=1,torbeam=1)

# Rch, zch, int_meas, rad_temp, abs_temp = A.run(fstart=-0.35, fend=0.35, Nf=10, zstart=-14, zend=14, Nz=10, torbeam=1)

Rch = Rch.reshape(24,8)
zch = zch.reshape(24,8)
abs_temp = abs_temp.reshape(24,8)

# int_meas = int_meas.reshape(24,8)
# rad_temp = rad_temp.reshape(24,8)


sio.savemat('data/results.mat', {'Rch':Rch, 'zch':zch, 'abs_temp':abs_temp})

# sio.savemat('data/results.mat', {'Rch':Rch, 'zch':zch, 'abs_temp':abs_temp, 'int_meas':int_meas, 'rad_temp':rad_temp})
