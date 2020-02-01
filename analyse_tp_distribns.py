'''
Python script to analyse the transition path probability distributions from the enhanced_kmc code

Daniel J. Sharpe
Jan 2020
'''

import numpy as np
from math import floor

class Analyse_tp_distribns(object):

    def __init__(self,stat,nbins,binw,bin_min):
        self.stat=stat
        self.nbins=nbins
        self.binw=binw
        self.bin_min=bin_min
        self.bin_max=self.bin_min+(self.nbins*self.binw)

    def get_hist_arr(self):
        hist_arr = np.zeros(self.nbins,dtype=int)
        with open("tp_distribns.dat","r") as tp_distribns_f:
            for line in tp_distribns_f.readlines():
                val=float(line.split()[stat+1])
                if val>self.bin_max or val<self.bin_min:
                    print "found bad value: ", val
                    raise RuntimeError
                hist_arr[int(floor((val-self.bin_min)/self.binw))] += 1
        return hist_arr

if __name__=="__main__":
    ### CHOOSE PARAMS ###

    # statistic to analyse
    # 0=dynamical activity, 1=time, 2=path prob, 3=entropy flow
    stat=1

    #binning params
    nbins=300
    binw=0.05
    bin_min=0.

    # run
    hist=Analyse_tp_distribns(stat,nbins,binw,bin_min)
    hist_arr = hist.get_hist_arr()
    print hist_arr
