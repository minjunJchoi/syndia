# Author : Minjun J. Choi (mjchoi@nfri.re.kr)
#
# Description : This code reads the KSTAR MDSplus server data
#
# Acknowledgement : Special thanks to Dr. Y.M. Jeon
#
import os

from MDSplus import Connection
# from MDSplus import DisconnectFromMds
# from MDSplus._mdsshr import MdsException

import numpy as np
import matplotlib.pyplot as plt

# Access MDSplus server to get ECEI operation parameters # Do not inherit Connection for multiprocessing compatibility
MDSPLUS_SERVER = os.environ.get('MDSPLUS_SERVER', 'mdsr.kstar.kfe.re.kr:8005')
mdsconn = Connection(MDSPLUS_SERVER)

class KstarEceInfo(object):
    def __init__(self, shot, clist):
        # super(KstarEceInfo,self).__init__(MDSPLUS_SERVER)  # call __init__ in Connection

        self.shot = shot
        self.clist = clist
        
        self.hn = 2

        self.channel_freq()

    def channel_freq(self):
        if self.shot < 12273:  # YEAR 2014
            freqECE = [110,111,112,113,114,115,116,117,118,119,120,121,123,124,125,126,127,128,129,130,131,132,133,134,138,139,140,141,142,143,144,145,146,147,148,149,151,152,153,154,155,156,157,158,159,160,161,162,164,165,167,168,169,170,171,172,173,174,175,176,177,178,182,183,184,185,186,187,188,189,190,191,192,193,195,196]
        elif (12273 <= self.shot) and (self.shot <= 14386):  # YEAR 2015
            freqECE = [110,111,112,113,114,115,116,117,118,119,120,121,123,124,125,126,127,128,129,130,131,132,133,134,138,139,140,141,142,143,144,145,146,147,148,149,151,152,153,154,155,156,157,158,159,160,161,162]
        elif (14386 < self.shot) and (self.shot <= 17356):  # YEAR 2016
            freqECE = [110,111,112,113,114,115,116,117,118,119,120,121,123,124,125,126,127,128,129,130,131,132,133,134,138,139,140,141,142,143,144,145,146,147,148,149,151,152,153,154,155,156,157,158,159,160,161,162,78,79,81,82,83,84,85,86,87,88,89,90,91,92,96,97,98,99,100,101,102,103,104,105,106,107,109,110]
        elif (17356 < self.shot) and (self.shot <= 19399):  # YEAR 2017
            freqECE = [110,111,112,113,114,115,116,117,118,119,120,121,123,124,125,126,127,128,129,130,131,132,133,134,138,139,140,141,142,143,144,145,146,147,148,149,151,152,153,154,155,156,157,158,159,160,161,162,78,79,81,82,83,84,85,86,87,88,89,90,91,92,96,97,98,99,100,101,102,103,104,105,106,107,109,110]
        elif (19400 <= self.shot):  # YEAR 2018
            freqECE = [110,111,112,113,114,115,116,117,118,119,120,121,123,124,125,126,127,128,129,130,131,132,133,134,138,139,140,141,142,143,144,145,146,147,148,149,151,152,153,154,155,156,157,158,159,160,161,162,78,79,81,82,83,84,85,86,87,88,89,90,91,92,96,97,98,99,100,101,102,103,104,105,106,107,109,110]

        ece_freq_dict = {}
        for i, f in enumerate(freqECE):
            chname = 'ECE{:02d}'.format(i+1)
            ece_freq_dict[chname] = f # GHz

        self.ece_freq = np.zeros(len(self.clist))
        for i, chname in enumerate(self.clist):
            self.ece_freq[i] = ece_freq_dict[chname]

    def channel_position(self):  # Needs updates ####################
        # get channel position either from MDSplus server or kstardata
        cnum = len(self.clist)
        self.rpos = np.arange(cnum, dtype=np.float64)  # R [m]
        self.zpos = np.zeros(cnum)  # z [m]
        self.apos = np.arange(cnum, dtype=np.float64)  # angle [rad]

        if ('ECE' == self.clist[0][0:3]):
            try:
                self.read_mds_rpos()
            except:
                print('FAIL to read {:s} rpos from MDS server'.format(self.clist[0]))
                ece_rpos = ece_pos.get_ece_pos(self.shot)
                for c in range(cnum):
                    self.rpos[c] = ece_rpos[self.clist[c]]
                print('rpos is obtained from kstardata (2nd harmonics cold resonance)')
  
    def read_mds_rpos(self):
        # MDSplus server tree
        tree = 'KSTAR'

        # open tree
        mdsconn.openTree(tree, self.shot)
        print('OPEN MDS tree {:s} to read rpos'.format(tree))
            
        # read rnode from MDSplus 
        cnum = len(self.clist)
        for c in range(cnum):
            # set rnode 
            if 'ECE' == self.clist[0][0:3]: # ECE
                rnode = '\{:s}:RPOS2ND'.format(self.clist[c])

            # read rnode
            rpos = mdsconn.get(rnode).data()
            if hasattr(rpos, "__len__"):
                self.rpos[c] = mdsconn.get(rnode).data()[0]
            else:
                self.rpos[c] = rpos

            # print('rpos of {:s} is read from MDS tree {:s}'.format(self.clist[c], tree))

        # close tree
        mdsconn.closeTree(tree, self.shot)


class NoPosMdsError(Exception):
    def __init__(self, msg='No position in MDSplus server'):
        self.msg = msg

    def __str__(self):
        return self.msg


if __name__ == "__main__":
    pass

    # g = KstarMds(shot=17245,clist=['neAVGM'])
    # g.get_data(trange=[0,10])
    # plt.plot(g.time, g.data[0,:], color='k')
    # plt.show()
    # g.close

# DisconnectFromMds(g.socket)
