'''
Python script to analyse the transition path probability distributions from the enhanced_kmc code

Daniel J. Sharpe
Jan 2020
'''

import numpy as np
# import matplotlib.pyplot as plt
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
        self.vals=None

    def get_hist_arr(self):
        hist_arr = np.zeros(self.nbins,dtype=int)
        vals=[]
        with open("tp_distribns.dat","r") as tp_distribns_f:
            for line in tp_distribns_f.readlines():
                val=float(line.split()[stat+1])
                vals.append(val)
                if not (val>self.bin_max or val<self.bin_min):
                    hist_arr[int(floor((val-self.bin_min)/self.binw))] += 1
                elif self.binall:
                    print "found bad value: ", val
                    raise RuntimeError
                self.ntpaths+=1
        self.vals=np.array(vals,dtype=float)
        return hist_arr

    def plot_hist(self,hist_arr):
        hist_arr=hist_arr.astype(np.float64)*1./float(self.ntpaths) # normalise
        bins=[self.bin_min+(i*self.binw) for i in range(self.nbins)]
        plt.bar(bins,hist_arr,self.binw,color='blue')
        plt.show()

    def calc_mfpt(self):
        return np.sum(self.vals)/float(self.ntpaths)

if __name__=="__main__":
    ### CHOOSE PARAMS ###

    # statistic to analyse
    # 0=dynamical activity, 1=time, 2=path prob, 3=entropy flow
    stat=1

    #binning params
    nbins=100
    binw=50.
    bin_min=0.
    binall=False # enforce that all values 

    # run
    calc_hist_obj=Analyse_tp_distribns(stat,nbins,binw,bin_min,binall)
    hist_arr = calc_hist_obj.get_hist_arr()
    print hist_arr
    print "total number of observed A<-B transition paths: ", calc_hist_obj.ntpaths
    print "total number of binned A<-B transition paths: ", np.sum(hist_arr)
    print "MFPT: ", calc_hist_obj.calc_mfpt()
    # plot
#    calc_hist_obj.plot_hist(hist_arr)
