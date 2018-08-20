from efm import *

A = EceFwdMod()

A.get_profile('data/g013728.003900', 'data/Te_3900.dat', 'data/ne_3900.dat')

A.set_channel(13728,['ECEI_G0101'])

A.rad_temp(fstart=0,fend=0,Nf=1,zstart=0,zend=0,Nz=1,ToR=0)
