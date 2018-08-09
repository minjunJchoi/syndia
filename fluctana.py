# Author : Minjun J. Choi (mjchoi@nfri.re.kr)
#
# Description : This code calculates cross power, coherence, cross phase, etc with fusion plasma diagnostics data
#
# Acknowledgement : Dr. S. Zoletnik and Prof. Y.-c. Ghim
#
# Last updated
#  2018.03.23 : version 0.10; even nfft -> odd nfft (for symmetry)


import matplotlib.pyplot as plt
from scipy import signal

from kstarecei import *
from kstarmds import *


CM = plt.cm.get_cmap('RdYlBu_r')


class FluctAna(object):
    def __init__(self):
        self.Dlist = []

    def add_data(self, D, shot, trange, clist, norm=1, atrange=[1.0, 1.01], res=0):
        clist = expand_clist(clist)

        D.get_data(shot, trange, clist, norm=norm, atrange=atrange, res=res)
        self.Dlist.append(D)

    def del_data(self, dnum):
        del self.Dlist[dnum]

    def list_data(self):
        for i in range(len(self.Dlist)):
            print '---- DATA SET # %d for [%g, %g] s ----' % (i, self.Dlist[i].trange[0], self.Dlist[i].trange[1])
            for j, c in enumerate(self.Dlist[i].clist):
                print '{%03d:%s}' % (j, c),
                if np.mod(j+1, 4) == 0 or j == len(self.Dlist[i].clist)-1:
                    print ''
            print ''
            # print '     # %d size : %s' % (i, self.Dlist[i].data.shape)

    def add_channel(self, dnum, clist):  # re-do fftbins after add channels
        old_clist = self.Dlist[dnum].clist

        # add channels (no duplicates)
        clist = expand_clist(clist)
        clist = [c for c in clist if c not in self.Dlist[dnum].clist]

        # add data
        time, data = self.Dlist[dnum].get_data(self.Dlist[dnum].shot, self.Dlist[dnum].trange, clist, norm=norm, atrange=atrange, res=res)
        self.Dlist[dnum].data = np.concatenate((self.Dlist[dnum].data, data), axis=0)

        # update clist
        self.Dlist[dnum].clist = old_clist + clist

        self.list_data()

    def del_channel(self, dnum, clist):
        clist = expand_clist(clist)

        for i in range(len(clist)):
            # find index to be deleted
            del_idx = [j for j, s in enumerate(self.Dlist[dnum].clist) if clist[i] in s]

            # delete data
            self.Dlist[dnum].data = np.delete(self.Dlist[dnum].data, del_idx, 0)

            # delete fftdata if it has
            if hasattr(self.Dlist[dnum], 'fftdata'):
                self.Dlist[dnum].fftdata = np.delete(self.Dlist[dnum].fftdata, del_idx, 0)

            # update clist
            self.Dlist[dnum].clist = [self.Dlist[dnum].clist[k] for k in range(len(self.Dlist[dnum].clist)) if k not in del_idx]

        self.list_data()

    def change_shot(self, dnum, shot):
        if type(dnum) is list:
            for i in dnum:
                self.Dlist[i].time, self.Dlist[i].data = self.Dlist[i].get_data(shot, self.Dlist[i].trange, self.Dlist[i].clist, norm=norm, atrange=atrange, res=res)
        else:
            self.Dlist[dnum].time, self.Dlist[dnum].data = self.Dlist[dnum].get_data(shot, self.Dlist[dnum].trange, self.Dlist[dnum].clist, norm=norm, atrange=atrange, res=res)


    def change_trange(self, dnum, trange):
        if type(dnum) is list:
            for i in dnum:
                self.Dlist[i].time, self.Dlist[i].data = self.Dlist[i].get_data(self.Dlist[i].shot, trange, self.Dlist[i].clist, norm=norm, atrange=atrange, res=res)
        else:
            self.Dlist[dnum].time, self.Dlist[dnum].data = self.Dlist[dnum].get_data(self.Dlist[dnum].shot, trange, self.Dlist[dnum].clist, norm=norm, atrange=atrange, res=res)

    def change_clist(self, dnum, clist):
        clist = expand_clist(clist)

        if type(dnum) is list:
            for i in dnum:
                self.Dlist[i].time, self.Dlist[i].data = self.Dlist[i].get_data(self.Dlist[i].shot, self.Dlist[i].trange, clist, norm=norm, atrange=atrange, res=res)
        else:
            self.Dlist[dnum].time, self.Dlist[dnum].data = self.Dlist[dnum].get_data(self.Dlist[dnum].shot, self.Dlist[dnum].trange, clist, norm=norm, atrange=atrange, res=res)

    def fftbins(self, nfft, window, overlap, detrend):
        # IN : self, data set number, nfft, window name, detrend or not
        # OUT : bins x N FFT of time series data; frequency axis

        self.list_data()

        for dnum in range(len(self.Dlist)):
            # get bins and window function
            tnum = len(self.Dlist[dnum].data[0,:])
            bins, win = fft_window(tnum, nfft, window, overlap)

            # make an x-axis # (default full length) 0~fN -fN~-f1
            dt = self.Dlist[dnum].time[1] - self.Dlist[dnum].time[0]  # time step
            ax = np.fft.fftfreq(nfft, d=dt)
            if np.mod(nfft, 2) == 0:  # even nfft
                ax = np.hstack([ax[0:(nfft/2)], -(ax[nfft/2]), ax[(nfft/2):nfft]])
            self.Dlist[dnum].ax = ax

            # make fftdata # (default full length) 0~fN -fN~-f1
            cnum = len(self.Dlist[dnum].data)
            if np.mod(nfft, 2) == 0:  # even nfft
                self.Dlist[dnum].fftdata = np.zeros((cnum, bins, nfft+1), dtype=np.complex_)
            else:  # odd nfft
                self.Dlist[dnum].fftdata = np.zeros((cnum, bins, nfft), dtype=np.complex_)

            for c in range(cnum):
                x = self.Dlist[dnum].data[c,:]

                for b in range(bins):
                    idx1 = int(b*np.fix(nfft*(1 - overlap)))
                    idx2 = idx1 + nfft

                    sx = x[idx1:idx2]

                    if detrend is 1:
                        sx = signal.detrend(sx, type='linear')
                        sx = signal.detrend(sx, type='constant')
                    else:
                        sx = signal.detrend(sx, type='constant')  # subtract mean

                    sx = sx * win  # apply window function

                    # get fft
                    if np.mod(nfft, 2) == 0:  # even nfft
                        fftdata = np.fft.fft(sx, n=nfft)/nfft  # divide by the length
                        fftdata = np.hstack([fftdata[0:(nfft/2)], np.conj(fftdata[nfft/2]), fftdata[(nfft/2):nfft]])
                        self.Dlist[dnum].fftdata[c,b,:] = fftdata
                    else: # odd nfft
                        self.Dlist[dnum].fftdata[c,b,:] = np.fft.fft(sx, n=nfft)/nfft  # divide by the length

            # update attributes
            if np.mod(nfft, 2) == 0:
                self.Dlist[dnum].nfft = nfft + 1
            else:
                self.Dlist[dnum].nfft = nfft
            self.Dlist[dnum].window = window
            self.Dlist[dnum].overlap = overlap
            self.Dlist[dnum].detrend = detrend
            self.Dlist[dnum].bins = bins
            self.Dlist[dnum].win = win

            print 'dnum %d fftbins %d with %s size %d overlap %g detrend %d' % (dnum, bins, window, nfft, overlap, detrend)

    def cross_power(self, done, dtwo, dc=0):
        # IN : data number one (ref), data number two (cmp), etc
        # OUT : x-axis (ax), y-axis (val)

        self.Dlist[dtwo].vkind = 'cross_power'

        rnum = len(self.Dlist[done].data)  # number of ref channels
        cnum = len(self.Dlist[dtwo].data)  # number of cmp channels
        bins = self.Dlist[dtwo].bins  # number of bins
        nfft = self.Dlist[dtwo].nfft  # length of fft of bin data
        win_factor = np.mean(self.Dlist[dtwo].win**2)  # window factors

        # reference channel names
        self.Dlist[dtwo].rname = []

        # half (0~fN) axis
        self.Dlist[dtwo].ax = self.Dlist[dtwo].ax[0:(nfft/2+1)]
        # value dimension
        val = np.zeros((cnum, bins, nfft), dtype=np.complex_)  # (full length for calculation)
        self.Dlist[dtwo].val = np.zeros((cnum, (nfft/2+1)))  # (half length for return)

        # calculation loop for multi channels
        for c in range(cnum):
            # reference channel number
            if rnum is 1:
                self.Dlist[dtwo].rname.append(self.Dlist[done].clist[0])
            else:
                #self.Dlist[dtwo].rname.append(self.Dlist[done].clist[c])
                if c + dc < 0:
                    self.Dlist[dtwo].rname.append(self.Dlist[done].clist[c])
                else:
                    self.Dlist[dtwo].rname.append(self.Dlist[done].clist[c+dc])

            # calculate cross power for each channel and each bins
            for b in range(bins):
                if rnum is 1:  # single reference channel
                    X = self.Dlist[done].fftdata[0,b,:]
                else:  # number of ref channels = number of cmp channels
                    #X = self.Dlist[done].fftdata[c,b,:]
                    if c + dc < 0:
                        X = self.Dlist[done].fftdata[c,b,:]
                    else:
                        X = self.Dlist[done].fftdata[c+dc,b,:]

                Y = self.Dlist[dtwo].fftdata[c,b,:]

                val[c,b,:] = 2*X*np.matrix.conjugate(Y) / win_factor  # product 2 for half return

            # average over bins
            Pxy = np.mean(val[c,:,:], 0)
            # result saved in val
            self.Dlist[dtwo].val[c,:] = np.abs(Pxy[0:(nfft/2+1)]).real

            # std saved in std

    def coherence(self, done, dtwo, dc=0):
        # IN : data number one (ref), data number two (cmp), etc
        # OUT : x-axis (ax), y-axis (val)

        self.Dlist[dtwo].vkind = 'coherence'

        rnum = len(self.Dlist[done].data)  # number of ref channels
        cnum = len(self.Dlist[dtwo].data)  # number of cmp channels
        bins = self.Dlist[dtwo].bins  # number of bins
        nfft = self.Dlist[dtwo].nfft  # length of fft of bin data

        # reference channel names
        self.Dlist[dtwo].rname = []

        # half (0~fN) axis
        self.Dlist[dtwo].ax = self.Dlist[dtwo].ax[0:(nfft/2+1)]
        # value dimension
        val = np.zeros((cnum, bins, nfft), dtype=np.complex_)  # (full length for calculation)
        self.Dlist[dtwo].val = np.zeros((cnum, (nfft/2+1)))  # (half length for return)

        # calculation loop for multi channels
        for c in range(cnum):
            # reference channel names
            if rnum is 1:
                self.Dlist[dtwo].rname.append(self.Dlist[done].clist[0])
            else:
                #self.Dlist[dtwo].rname.append(self.Dlist[done].clist[c])
                if c + dc < 0:
                    self.Dlist[dtwo].rname.append(self.Dlist[done].clist[c])
                else:
                    self.Dlist[dtwo].rname.append(self.Dlist[done].clist[c+dc])

            # calculate cross power for each channel and each bins
            for b in range(bins):
                if rnum is 1:  # single reference channel
                    X = self.Dlist[done].fftdata[0,b,:]
                else:  # number of ref channels = number of cmp channels
                    #X = self.Dlist[done].fftdata[c,b,:]
                    if c + dc < 0:
                        X = self.Dlist[done].fftdata[c,b,:]
                    else:
                        X = self.Dlist[done].fftdata[c+dc,b,:]

                Y = self.Dlist[dtwo].fftdata[c,b,:]

                Pxx = X * np.matrix.conjugate(X)
                Pyy = Y * np.matrix.conjugate(Y)
                val[c,b,:] = X*np.matrix.conjugate(Y) / np.sqrt(Pxx*Pyy)
                # saturated data gives zero Pxx!!

            # average over bins
            Gxy = np.mean(val[c,:,:], 0)
            # results saved in val
            self.Dlist[dtwo].val[c,:] = np.abs(Gxy[0:(nfft/2+1)]).real

    def cross_phase(self, done, dtwo, dc=0):
        # IN : data number one (ref), data number two (cmp)
        # OUT : x-axis (ax), y-axis (val)

        self.Dlist[dtwo].vkind = 'cross_phase'

        rnum = len(self.Dlist[done].data)  # number of ref channels
        cnum = len(self.Dlist[dtwo].data)  # number of cmp channels
        bins = self.Dlist[dtwo].bins  # number of bins
        nfft = self.Dlist[dtwo].nfft  # length of fft of bin data

        # reference channel names
        self.Dlist[dtwo].rname = []

        # distance
        self.Dlist[dtwo].dist = np.zeros(cnum)

        # half (0~fN) axis
        self.Dlist[dtwo].ax = self.Dlist[dtwo].ax[0:(nfft/2+1)]
        # value dimension
        val = np.zeros((cnum, bins, nfft), dtype=np.complex_)  # (full length for calculation)
        self.Dlist[dtwo].val = np.zeros((cnum, (nfft/2+1)))  # (half length for return)

        # calculation loop for multi channels
        for c in range(cnum):
            # reference channel number and distance between ref and cmp channels
            if rnum is 1:
                self.Dlist[dtwo].rname.append(self.Dlist[done].clist[0])
                self.Dlist[dtwo].dist[c] = np.sqrt((self.Dlist[dtwo].rpos[c] - self.Dlist[done].rpos[0])**2 + \
                (self.Dlist[dtwo].zpos[c] - self.Dlist[done].zpos[0])**2)
            else:
                if c + dc < 0:
                    self.Dlist[dtwo].rname.append(self.Dlist[done].clist[c])
                    self.Dlist[dtwo].dist[c] = np.sqrt((self.Dlist[dtwo].rpos[c] - self.Dlist[done].rpos[c])**2 + \
                    (self.Dlist[dtwo].zpos[c] - self.Dlist[done].zpos[c])**2)
                else:
                    self.Dlist[dtwo].rname.append(self.Dlist[done].clist[c+dc])
                    self.Dlist[dtwo].dist[c] = np.sqrt((self.Dlist[dtwo].rpos[c] - self.Dlist[done].rpos[c+dc])**2 + \
                    (self.Dlist[dtwo].zpos[c] - self.Dlist[done].zpos[c+dc])**2)

            # calculate cross power for each channel and each bins
            for b in range(bins):
                if rnum is 1:  # single reference channel
                    X = self.Dlist[done].fftdata[0,b,:]
                else:  # number of ref channels = number of cmp channels
                    #X = self.Dlist[done].fftdata[c,b,:]
                    if c + dc < 0:
                        X = self.Dlist[done].fftdata[c,b,:]
                    else:
                        X = self.Dlist[done].fftdata[c+dc,b,:]
                Y = self.Dlist[dtwo].fftdata[c,b,:]

                val[c,b,:] = X*np.matrix.conjugate(Y)

            # average over bins
            Pxy = np.mean(val[c,:,:], 0)
            # result saved in val
            self.Dlist[dtwo].val[c,:] = np.arctan2(Pxy[0:(nfft/2+1)].imag, Pxy[0:(nfft/2+1)].real).real

            # std saved in std

    def mplot(self, dnum, cnum, type='time', **kwargs):
        if 'ylimits' in kwargs: ylimits = kwargs['ylimits']
        if 'xlimits' in kwargs: xlimits = kwargs['xlimits']

        for i in cnum:
            if type is 'time':
                pbase = self.Dlist[dnum].time
                pdata = self.Dlist[dnum].data[i,:]
            elif type is 'val':
                pbase = self.Dlist[dnum].ax/1000
                pdata = self.Dlist[dnum].val[i,:].real
                rname = self.Dlist[dnum].rname[i]
                if self.Dlist[dnum].vkind is 'coherence':
                    plt.axhline(y=1/np.sqrt(self.Dlist[dnum].bins), color='r')
            pname = self.Dlist[dnum].clist[i]
            pshot = self.Dlist[dnum].shot

            plt.plot(pbase, pdata)  # plot

            if 'ylimits' in kwargs:  # ylimits
                plt.ylim([ylimits[0], ylimits[1]])
            if 'xlimits' in kwargs:  # xlimits
                plt.xlim([xlimits[0], xlimits[1]])
            else:
                plt.xlim([pbase[0], pbase[-1]])

            if type is 'time':
                plt.title('SHOT = %d, %s' % (pshot, pname), fontsize=10)
            elif type is 'val':
                plt.title('SHOT = %d, %s-%s' % (pshot, rname, pname), fontsize=10)

            if type is 'time':
                plt.xlabel('Time [s]')
                plt.ylabel('ECEI signal')
            elif type is 'val' and self.Dlist[dnum].vkind is 'cross_power':
                plt.xlabel('Frequency [kHz]')
                plt.ylabel('Cross power spectral density')
                plt.yscale('log')
            elif type is 'val' and self.Dlist[dnum].vkind is 'coherence':
                plt.xlabel('Frequency [kHz]')
                plt.ylabel('Coherence')
            elif type is 'val' and self.Dlist[dnum].vkind is 'cross_phase':
                plt.xlabel('Frequency [kHz]')
                plt.ylabel('Cross phase [rad]')

            plt.show()

    def oplot(self, dnum, cnum, type='time', **kwargs):
        if 'ylimits' in kwargs: ylimits = kwargs['ylimits']
        if 'xlimits' in kwargs: xlimits = kwargs['xlimits']

        for i in cnum:
            if type is 'time':
                pbase = self.Dlist[dnum].time
                pdata = self.Dlist[dnum].data[i,:]
            elif type is 'val':
                pbase = self.Dlist[dnum].ax/1000
                pdata = self.Dlist[dnum].val[i,:].real
                rname = self.Dlist[dnum].rname[i]
                if i == 0 and self.Dlist[dnum].vkind is 'coherence':
                    plt.axhline(y=1/np.sqrt(self.Dlist[dnum].bins), color='r')
            pname = self.Dlist[dnum].clist[i]

            plt.plot(pbase, pdata)

            if type is 'time':
                print 'dnum %d : channel %s is plotted' % (dnum, pname)
            elif type is 'val':
                print 'dnum %d : calculation %s-%s is plotted' % (dnum, rname, pname)

            if 'ylimits' in kwargs: # ylimits
                plt.ylim([ylimits[0], ylimits[1]])
            if 'xlimits' in kwargs: # xlimits
                plt.xlim([xlimits[0], xlimits[1]])
            else:
                plt.xlim([pbase[0], pbase[-1]])

            if type is 'time':
                plt.xlabel('Time [s]')
                plt.ylabel('ECEI signal')
            elif type is 'val' and self.Dlist[dnum].vkind is 'cross_power':
                plt.xlabel('Frequency [kHz]')
                plt.ylabel('Cross power spectral density')
                plt.yscale('log')
            elif type is 'val' and self.Dlist[dnum].vkind is 'coherence':
                plt.xlabel('Frequency [kHz]')
                plt.ylabel('Coherence')
            elif type is 'val' and self.Dlist[dnum].vkind is 'cross_phase':
                plt.xlabel('Frequency [kHz]')
                plt.ylabel('Cross phase [rad]')

        plt.show()

    def spec(self, dnum, cnum, nfft=516, climits=[-160, -60], **kwargs):
        if 'flimits' in kwargs: flimits = kwargs['flimits']*1000
        if 'xlimits' in kwargs: xlimits = kwargs['xlimits']

        fs = self.Dlist[dnum].fs
        nov = nfft*0.9

        for i in cnum:
            pbase = self.Dlist[dnum].time
            pdata = self.Dlist[dnum].data[i,:]
            pname = self.Dlist[dnum].clist[i]
            pshot = self.Dlist[dnum].shot

            pxx, freq, time, cax = plt.specgram(pdata, NFFT=nfft, Fs=fs, noverlap=nov, xextent=[pbase[0], pbase[-1]],
                                                vmin=climits[0], vmax=climits[1], cmap=CM)  # spectrum

            plt.colorbar(cax)

            if 'flimits' in kwargs:  # flimits
                plt.ylim([flimits[0], flimits[1]])
            if 'xlimits' in kwargs:  # xlimits
                plt.ylim([xlimits[0], xlimits[1]])
            else:
                plt.xlim([pbase[0], pbase[-1]])

            plt.title(pname, fontsize=10)  # labeling
            plt.xlabel('Time [s]')
            plt.ylabel('Frequency [Hz]')

            plt.show()

    def iplot(self, dnum, snum=0, vlimits=[-0.1, 0.1], **kwargs):
        # keyboard interactive iplot
        # (intp='none', climits=[-0.1,0.1], **kwargs)

        # data filtering

        c = raw_input('automatic, or manual [a,m]: ')
        tidx1 = 0  # starting index
        if c == 'a':
            # make axes
            fig = plt.figure(facecolor='w', figsize=(5,10))
            ax1 = fig.add_axes([0.1, 0.75, 0.7, 0.2])  # [left bottom width height]
            ax2 = fig.add_axes([0.1, 0.1, 0.7, 0.60])
            ax3 = fig.add_axes([0.83, 0.1, 0.03, 0.6])
            axs = [ax1, ax2, ax3]

            tstep = int(input('time step [idx]: '))  # jumping index # tstep = 10
            for tidx in range(tidx1, len(self.Dlist[dnum].time), tstep):
                # prepare data
                pdata = self.Dlist[dnum].data[:,tidx]

                # position
                rpos = self.Dlist[dnum].rpos[:]
                zpos = self.Dlist[dnum].zpos[:]

                # plot
                axs[0].cla()
                axs[1].cla()
                axs[2].cla()
                plt.ion()

                axs[0].plot(self.Dlist[dnum].time, self.Dlist[dnum].data[snum,:])  # ax1.hold(True)
                axs[0].axvline(x=self.Dlist[dnum].time[tidx], color='g')
                sc = axs[1].scatter(rpos, zpos, 500, pdata, marker='s', vmin=vlimits[0], vmax=vlimits[1], cmap=CM)
                axs[1].set_aspect('equal')
                plt.colorbar(sc, cax=axs[2])

                axs[1].set_xlabel('R [m]')
                axs[1].set_ylabel('z [m]')
                axs[1].set_title('ECE image')

                plt.show()
                plt.pause(0.1)

            plt.ioff()
            plt.close()

        elif c == 'm':
            tidx = tidx1
            while True:
                # make axes
                fig = plt.figure(facecolor='w', figsize=(5,10))
                ax1 = fig.add_axes([0.1, 0.75, 0.7, 0.2])  # [left bottom width height]
                ax2 = fig.add_axes([0.1, 0.1, 0.7, 0.60])
                ax3 = fig.add_axes([0.83, 0.1, 0.03, 0.6])
                axs = [ax1, ax2, ax3]

                # prepare data
                pdata = self.Dlist[dnum].data[:,tidx]

                # position
                rpos = self.Dlist[dnum].rpos[:]
                zpos = self.Dlist[dnum].zpos[:]

                # plot
                axs[0].cla()
                axs[1].cla()
                axs[2].cla()
                plt.ion()

                axs[0].plot(self.Dlist[dnum].time, self.Dlist[dnum].data[snum,:])  # ax1.hold(True)
                axs[0].axvline(x=self.Dlist[dnum].time[tidx], color='g')
                sc = axs[1].scatter(rpos, zpos, 500, pdata, marker='s', vmin=vlimits[0], vmax=vlimits[1], cmap=CM)
                axs[1].set_aspect('equal')
                plt.colorbar(sc, cax=axs[2])

                axs[1].set_xlabel('R [m]')
                axs[1].set_ylabel('z [m]')
                axs[1].set_title('ECE image')

                plt.show()

                k = raw_input('set time step [idx][+,-,0]: ')
                try:
                    tstep = int(k)
                    if tstep == 0:
                        plt.ioff()
                        plt.close()
                        break
                except:
                    pass

                if tidx + tstep < len(self.Dlist[dnum].time) - 1 and 0 < tidx + tstep:
                    tidx = tidx + tstep

                plt.ioff()
                plt.close()

        self.Dlist[dnum].pdata = pdata

    def cplot(self, dnum, snum=0, frange=[0, 100], vlimits=[0, 1], **kwargs):
        if 'ylimits' in kwargs: ylimits = kwargs['ylimits']
        if 'xlimits' in kwargs: xlimits = kwargs['xlimits']
        # calculate mean coherence image
        # or cross power rms image
        # or group velocity image

        # axis
        pbase = self.Dlist[dnum].ax/1000  # [kHz]

        # fidx
        idx = np.where((pbase >= frange[0])*(pbase <= frange[1]))
        idx1 = int(idx[0][0])
        idx2 = int(idx[0][-1]+1)

        # data
        if self.Dlist[dnum].vkind is 'cross_power':  # rms
            pdata = np.sqrt(np.sum(self.Dlist[dnum].val[:,idx1:idx2], 1))
        elif self.Dlist[dnum].vkind is 'coherence':  # mean coherence
            pdata = np.mean(self.Dlist[dnum].val[:,idx1:idx2], 1)
        elif self.Dlist[dnum].vkind is 'cross_phase':  # group velocity
            base = self.Dlist[dnum].ax[idx1:idx2]  # [Hz]
            pdata = np.zeros(len(self.Dlist[dnum].val))
            for c in range(len(self.Dlist[dnum].val)):
                data = self.Dlist[dnum].val[c,idx1:idx2]
                pfit = np.polyfit(base, data, 1)
                fitdata = np.polyval(pfit, base)
                pdata[c] = 2*np.pi*self.Dlist[dnum].dist[c]/pfit[0]/1000.0  # [km/s]

                # chisq = np.sum((data - fitdata)**2)
                if c == snum:
                    sbase = base/1000  # [kHz]
                    sdata = fitdata

        # position
        rpos = self.Dlist[dnum].rpos[:]
        zpos = self.Dlist[dnum].zpos[:]

        # prepare figure
        fig = plt.figure(facecolor='w', figsize=(5,10))
        ax1 = fig.add_axes([0.1, 0.75, 0.7, 0.2])  # [left bottom width height]
        ax2 = fig.add_axes([0.1, 0.1, 0.7, 0.60])
        ax3 = fig.add_axes([0.83, 0.1, 0.03, 0.6])
        axs = [ax1, ax2, ax3]

        # sample plot
        axs[0].plot(pbase, self.Dlist[dnum].val[snum,:])  # ax1.hold(True)
        if self.Dlist[dnum].vkind is 'cross_phase':
            axs[0].plot(sbase, sdata)
        axs[0].axvline(x=pbase[idx1], color='g')
        axs[0].axvline(x=pbase[idx2], color='g')

        if self.Dlist[dnum].vkind is 'cross_power':
            axs[0].set_yscale('log')
        if 'ylimits' in kwargs: # ylimits
            axs[0].set_ylim([ylimits[0], ylimits[1]])
        if 'xlimits' in kwargs: # xlimits
            axs[0].set_xlim([xlimits[0], xlimits[1]])
        else:
            axs[0].set_xlim([pbase[0], pbase[-1]])

        # pdata plot
        sc = axs[1].scatter(rpos, zpos, 500, pdata, marker='s', vmin=vlimits[0], vmax=vlimits[1], cmap=CM)
        axs[1].set_aspect('equal')

        # color bar
        plt.colorbar(sc, cax=axs[2])

        axs[1].set_xlabel('R [m]')
        axs[1].set_ylabel('z [m]')
        if self.Dlist[dnum].vkind is 'cross_power':
            axs[1].set_title('Cross power rms')
        elif self.Dlist[dnum].vkind is 'coherence':
            axs[1].set_title('Coherence mean')
        elif self.Dlist[dnum].vkind is 'cross_phase':
            axs[1].set_title('Group velocity [km/s]')

        self.Dlist[dnum].pdata = pdata


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
                    exp_clist.append(clist[c][0:6] + '%02d' % v + '%02d' % f)
        else:
            exp_clist.append(clist[c])
    clist = exp_clist

    return clist


def fft_window(tnum, nfft, window, overlap):
    # IN : full length of time series, nfft, window name, overlap ratio
    # OUT : bins, 1 x nfft window function

    # use overlapping
    bins = int(np.fix((tnum/nfft - overlap)/(1 - overlap)))

    # window function
    if window is 'rectwin':  # overlap = 0.5
        win = np.ones(nfft)
    elif window is 'hann':  # overlap = 0.5
        win = np.hanning(nfft)
    elif window is 'hamm':  # overlap = 0.5
        win = np.hamming(nfft)
    elif window is 'kaiser':  # overlap = 0.62
        win = np.kaiser(nfft, beta=3)
    elif window is 'HFT248D':  # overlap = 0.84
        z = 2*np.pi/nfft*np.arange(0,nfft)
        win = 1 - 1.985844164102*np.cos(z) + 1.791176438506*np.cos(2*z) - 1.282075284005*np.cos(3*z) + \
            0.667777530266*np.cos(4*z) - 0.240160796576*np.cos(5*z) + 0.056656381764*np.cos(6*z) - \
            0.008134974479*np.cos(7*z) + 0.000624544650*np.cos(8*z) - 0.000019808998*np.cos(9*z) + \
            0.000000132974*np.cos(10*z)

    return bins, win
