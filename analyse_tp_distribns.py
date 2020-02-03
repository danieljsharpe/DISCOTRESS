'''
Python script to analyse the transition path probability distributions from the enhanced_kmc code

Daniel J. Sharpe
Jan 2020
'''

import numpy as np
import matplotlib.pyplot as plt
from math import floor

class Analyse_tp_distribns(object):

    def __init__(self,stat,nbins,binw,bin_min,binall):
        self.stat=stat
        self.nbins=nbins
        self.binw=binw
        self.bin_min=bin_min
        self.bin_max=self.bin_min+(self.nbins*self.binw)
        self.binall=binall
        self.ntpaths=0

    def get_hist_arr(self):
        hist_arr = np.zeros(self.nbins,dtype=int)
        with open("tp_distribns.dat","r") as tp_distribns_f:
            for line in tp_distribns_f.readlines():
                val=float(line.split()[stat+1])
                if not (val>self.bin_max or val<self.bin_min):
                    hist_arr[int(floor((val-self.bin_min)/self.binw))] += 1
                elif self.binall:
                    print "found bad value: ", val
                    raise RuntimeError
                self.ntpaths+=1
        return hist_arr

    def plot_hist(self,hist_arr):
        hist_arr=hist_arr.astype(np.float64)*1./float(self.ntpaths) # normalise
        bins=[self.bin_min+(i*self.binw) for i in range(self.nbins)]
        plt.bar(bins,hist_arr,self.binw,color='blue')
        plt.show()

if __name__=="__main__":
    ### CHOOSE PARAMS ###

    # statistic to analyse
    # 0=dynamical activity, 1=time, 2=path prob, 3=entropy flow
    stat=1

    #binning params
    nbins=300
    binw=0.01
    bin_min=0.
    binall=True # enforce that all values 

    # run
    calc_hist_obj=Analyse_tp_distribns(stat,nbins,binw,bin_min,binall)
    hist_arr = calc_hist_obj.get_hist_arr()
    print hist_arr
    print "total number of observed A<-B transition paths: ", calc_hist_obj.ntpaths
    print "total number of binned A<-B transition paths: ", np.sum(hist_arr)
    # plot
    calc_hist_obj.plot_hist(hist_arr)
