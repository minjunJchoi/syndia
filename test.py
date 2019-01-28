import pickle
import scipy.io as sio
import matplotlib.pyplot as plt
import time

from efm import *

geqdsk_fn = '/home/users/mjchoi/syndia/data/g021693.005000'
Te_fn = '/home/users/mjchoi/syndia/data/21693_5000_Te.dat' # normalized psi(:), Te(:) [keV]
ne_fn = '/home/users/mjchoi/syndia/data/21693_5000_ne.dat' # normalized psi(:), ne(:) [1e19 m^-3]

A = EceFwdMod()

A.set_profile(geqdsk_fn, Te_fn, ne_fn)


## compare TORBEAM and ray tracing
# ray tracing parameters are tuned against TORBEAM
#A.set_channel(19323,['ECEI_G0101'])

#st = time.time()
#RchTB, zchTB, _, _, _ = A.run(fstart=0,fend=0,Nf=1,zstart=0,zend=0,Nz=1,torbeam=1)
#print 'TB time = {}'.format(time.time() - st)

#st = time.time()
#RchRT, zchRT, _, _, _ = A.run(fstart=0,fend=0,Nf=1,zstart=0,zend=0,Nz=1,torbeam=0)
#print 'RT time = {}'.format(time.time() - st)



## channel posistion and Te for calibration
A.set_channel(21693,['ECEI_GR1201'])
Rch, zch, _, _, abs_temp = A.run(fstart=0,fend=0,Nf=1,zstart=0,zend=0,Nz=1,pstart=7.8,pend=-2,pint=0.5,Rinit=2.39,torbeam=0)

Rch = Rch.reshape(24,8)
zch = zch.reshape(24,8)
abs_temp = abs_temp.reshape(24,8)
sio.savemat('data/results_21693_5000_GR.mat', {'Rch':Rch, 'zch':zch, 'abs_temp':abs_temp})



## synthetic ECE image generation or radiation temperature check
# A.set_channel(13728,['ECEI_G0101-2408'])
# Rch, zch, int_meas, rad_temp, abs_temp = A.run(fstart=-0.35, fend=0.35, Nf=10, zstart=-14, zend=14, Nz=10, torbeam=1)

# Rch = Rch.reshape(24,8)
# zch = zch.reshape(24,8)
# int_meas = int_meas.reshape(24,8)
# rad_temp = rad_temp.reshape(24,8)
# abs_temp = abs_temp.reshape(24,8)
# sio.savemat('data/results.mat', {'Rch':Rch, 'zch':zch, 'int_meas':int_meas, 'rad_temp':rad_temp, 'abs_temp':abs_temp})



# sio.savemat('data/results.mat', {'Rch':Rch, 'zch':zch, 'abs_temp':abs_temp, 'int_meas':int_meas, 'rad_temp':rad_temp})
