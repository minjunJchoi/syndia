import os, sys
sys.path.insert(0, os.pardir)
import pickle
import scipy.io as sio
import matplotlib.pyplot as plt
import time

from efm import *

geqdsk_fn = '/home/users/mjchoi/syndia/run/data/g026896.005500'
Te_fn = '/home/users/mjchoi/syndia/run/data/26896_5500_Te.dat' # normalized psi(:), Te(:) [keV]
ne_fn = '/home/users/mjchoi/syndia/run/data/26896_5500_ne.dat' # normalized psi(:), ne(:) [1e19 m^-3]

A = EceFwdMod()

A.set_profile(geqdsk_fn, Te_fn, ne_fn)


## check profiles
Raxis = np.arange(1.8,2.25,0.01)
zaxis = np.zeros(len(Raxis))
Te = A.pf.F_Te(Raxis, zaxis)/(1.602*1e-19)/1000 # [J] -> [keV]
ne = A.pf.F_ne(Raxis, zaxis)/1e19 # [m-3] -> 1e19 [m-3]
plt.plot(Raxis, Te, 'k', label='Te [keV]')
plt.plot(Raxis, ne, 'r', label='ne [1e19 m-3]')
plt.xlabel('R [m]')
plt.legend()
plt.show()


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
shot = 26896
clist = ['ECEI_GT0101-2408']
fname = '/home/users/mjchoi/syndia/run/data/ecei_pos_{:d}_{:s}.pkl'.format(shot, clist[0])
A.set_channel(shot, clist)

Rch, zch, _, rad_temp, abs_temp = A.run(fstart=0,fend=0,Nf=1,zstart=0,zend=0,Nz=1,torbeam=1)

with open(fname, 'wb') as fout:
    pickle.dump([clist, Rch, zch, rad_temp, abs_temp], fout)




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
