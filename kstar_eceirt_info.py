"""
KSTAR ECEI RT information
Repo : https://github.com/minjunJchoi/syndia
Author : Minjun J. Choi (mjchoi@nfri.re.kr)
Acknowledgements : Jaehyun Lee
"""
#!/usr/bin/env python2.7

import os
import h5py
import numpy as np
from MDSplus import Connection

ECEI_TREE = 'ECEI'

# Access MDSplus server to get ECEI operation parameters # Do not inherit Connection for multiprocessing compatibility
MDSPLUS_SERVER = os.environ.get('MDSPLUS_SERVER', 'mdsr.kstar.kfe.re.kr:8005')
mdsconn = Connection(MDSPLUS_SERVER)

class KstarEceiRemoteInfo(object):
    def __init__(self, shot, clist):

        self.shot = shot

        self.clist = expand_clist(clist)
        cnum = len(self.clist)

        self.cnidx1 = 7
        self.dev = self.clist[0][5:7]

        # read operation parameters from MDSplus
        mdsconn.openTree(ECEI_TREE, self.shot)

        hn_node = '\\{0}::TOP.ECEI_{1}:{2}_MODE'.format(ECEI_TREE, self.dev, self.dev) 
        self.hn = mdsconn.get(hn_node).data()
        self.hn = int(self.hn)

        itf_node = '\\{0}::TOP:{1}'.format(ECEI_TREE, 'ECEI_I_TF') 
        self.itf = mdsconn.get(itf_node).data() # [kA]
        self.bt = self.itf*0.0995556 # [kA] -> [T]

        lo_node = '\\{0}::TOP.ECEI_{1}:{2}_LOFREQ'.format(ECEI_TREE, self.dev, self.dev) 
        self.lo = mdsconn.get(lo_node).data() # [GHz]

        sf_node = '\\{0}::TOP.ECEI_{1}:{2}_LENSFOCUS'.format(ECEI_TREE, self.dev, self.dev) 
        self.sf = mdsconn.get(sf_node).data() # [mm]

        sz_node = '\\{0}::TOP.ECEI_{1}:{2}_LENSZOOM'.format(ECEI_TREE, self.dev, self.dev) 
        self.sz = mdsconn.get(sz_node).data() # [mm]

        mdsconn.closeTree(ECEI_TREE, self.shot)
    
        # get vn, fn numbers and ece frequency [GHz]
        self.vn_list = np.zeros(cnum, dtype='int16')
        self.fn_list = np.zeros(cnum, dtype='int16')
        self.ece_freq = np.zeros(cnum)
        for cn, chname in enumerate(self.clist):
            self.vn_list[cn] = int(chname[(self.cnidx1):(self.cnidx1+2)])
            self.fn_list[cn] = int(chname[(self.cnidx1+2):(self.cnidx1+4)])
            self.ece_freq[cn] = (self.fn_list[cn] - 1)*0.9 + 2.6 + self.lo

    def get_abcd(self, sf, sz, Rinit):
        # ABCD matrix
        if self.dev == 'GT':
            sp = 2300 - Rinit*1000
            abcd = np.array([[1,sp+(2025-sz)],[0,1]]).dot(
                   np.array([[1,0],[(1.52-1)/(-1000*1),1.52/1]])).dot(
                   np.array([[1,160],[0,1]])).dot(
                   np.array([[1,0],[(1-1.52)/(1000*1.52),1/1.52]])).dot(
                   np.array([[1,2280-(2025+160-sz)],[0,1]])).dot(
                   np.array([[1,0],[(1.52-1)/(1000*1),1.52/1]])).dot(
                   np.array([[1,20],[0,1]])).dot(
                   np.array([[1,0],[0,1/1.52]])).dot(
                   np.array([[1,(4343-sf)-(2280+20)],[0,1]])).dot(
                   np.array([[1,0],[(1.52-1)/(-1200*1),1.52/1]])).dot(
                   np.array([[1,140],[0,1]])).dot(
                   np.array([[1,0],[(1-1.52)/(1200*1.52),1/1.52]])).dot(
                   np.array([[1,4520-(4343+140-sf)],[0,1]])).dot(
                   np.array([[1,0],[0,1.52/1]])).dot(
                   np.array([[1,30],[0,1]])).dot(
                   np.array([[1,0],[0,1/1.52]])).dot(
                   np.array([[1,4940-(4520+30)],[0,1]]))
        elif self.dev == 'GR':
            sp = 2300 - Rinit*1000
            abcd = np.array([[1,sp+(2025-sz)],[0,1]]).dot(
                   np.array([[1,0],[(1.52-1)/(-1000),1.52/1]])).dot(
                   np.array([[1,160],[0,1]])).dot(
                   np.array([[1,0],[(1-1.52)/(1000*1.52),1/1.52]])).dot(
                   np.array([[1,2280-(2025+160-sz)],[0,1]])).dot(
                   np.array([[1,0],[(1.52-1)/(1000*1),1.52/1]])).dot(
                   np.array([[1,20],[0,1]])).dot(
                   np.array([[1,0],[0,1/1.52]])).dot(
                   np.array([[1,4343-(2280+20)-sf],[0,1]])).dot(
                   np.array([[1,0],[(1.52-1)/(-1200),1.52/1]])).dot(
                   np.array([[1,140],[0,1]])).dot(
                   np.array([[1,0],[(1-1.52)/(1200*1.52),1/1.52]])).dot(
                   np.array([[1,4520-(4343+140-sf)],[0,1]])).dot(
                   np.array([[1,0],[0,1.52]])).dot(
                   np.array([[1,30],[0,1]])).dot(
                   np.array([[1,0],[0,1/1.52]])).dot(
                   np.array([[1,4940-(4520+30)],[0,1]]))
        elif self.dev == 'HT':
            sp = 2300 - Rinit*1000
            abcd = np.array([[1,sp+2553.01],[0,1]]).dot(
                   np.array([[1,0],[(1.526-1)/(-695*1),1.526/1]])).dot(
                   np.array([[1,150],[0,1]])).dot(
                   np.array([[1,0],[0,1/1.526]])).dot(
                   np.array([[1,4500.41-(2553.01+150)-sz],[0,1]])).dot(
                   np.array([[1,0],[0,1.526/1]])).dot(
                   np.array([[1,40],[0,1]])).dot(
                   np.array([[1,0],[(1-1.526)/(-515*1.526),1/1.526]])).dot(
                   np.array([[1,6122.41-(4500.41+40-sz)-sf],[0,1]])).dot(
                   np.array([[1,0],[0,1.52/1]])).dot(
                   np.array([[1,150],[0,1]])).dot(
                   np.array([[1,0],[(1-1.526)/(630*1.526),1/1.526]])).dot(
                   np.array([[1,6478.41-(6122.41+150-sf)],[0,1]])).dot(
                   np.array([[1,0],[0,1.526/1]])).dot(
                   np.array([[1,40],[0,1]])).dot(
                   np.array([[1,0],[0,1/1.526]])).dot(
                   np.array([[1,7161.01-(6478.41+40)],[0,1]]))

        return abcd


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
        elif 'ECEI' in clist[c] and len(clist[c]) == 16: # since 2018
            vi = int(clist[c][7:9])
            fi = int(clist[c][9:11])
            vf = int(clist[c][12:14])
            ff = int(clist[c][14:16])

            for v in range(vi, vf+1):
                for f in range(fi, ff+1):
                    exp_clist.append(clist[c][0:7] + '%02d' % v + '%02d' % f)
        else:
            exp_clist.append(clist[c])
    clist = exp_clist

    return clist
