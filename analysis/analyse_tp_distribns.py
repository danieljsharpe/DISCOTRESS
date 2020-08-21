'''
Python script to analyse the transition path probability distributions from the enhanced_kmc code

Daniel J. Sharpe
Jan 2020
'''

import numpy as np
import matplotlib.pyplot as plt
from math import floor
from math import sqrt

class Analyse_tp_distribns(object):

    def __init__(self,stat,nbins,binw,bin_min,binall,logvals):
        self.stat=stat
        self.nbins=nbins
        self.binw=binw
        self.bin_min=bin_min
        self.bin_max=self.bin_min+(self.nbins*self.binw)
        self.binall=binall
        self.logvals=logvals
        self.ntpaths=0
        self.vals=None

    def get_hist_arr(self):
        hist_arr = np.zeros(self.nbins,dtype=int)
        vals=[]
        with open("tp_distribns.dat","r") as tp_distribns_f:
            for line in tp_distribns_f.readlines():
                val=float(line.split()[stat+1])
                if self.logvals: val=np.log10(val)
                vals.append(val)
                if not (val>self.bin_max or val<self.bin_min):
                    hist_arr[int(floor((val-self.bin_min)/self.binw))] += 1
                elif self.binall:
                    print "found bad value: ", val
                    raise RuntimeError
                self.ntpaths+=1
        self.vals=np.array(vals,dtype=float)
        return hist_arr

    def plot_hist(self,hist_arr,nxticks,nyticks,ymax,figfmt="pdf"):
        hist_arr=hist_arr.astype(np.float64)*1./float(self.ntpaths) # normalise
        bins=[self.bin_min+(i*self.binw) for i in range(self.nbins)]
        plt.figure(figsize=(10.,7.)) # size in inches
        plt.bar(bins,hist_arr,self.binw,color='blue')
        if self.logvals:
            plt.xlabel("$\log_{10}(t_\mathrm{FPT}\ /\ \mathrm{ns})$",fontsize=42)
            plt.ylabel("$p ( \log_{10} t_\mathrm{FPT} )$",fontsize=42)
        else:
            plt.xlabel("$t_\mathrm{FPT}$",fontsize=42)
            plt.ylabel("p(t_\mahtrm{FPT})$",fontsize=42)
        ax = plt.gca()
        ax.set_xlim([self.bin_min,self.bin_max])
        ax.set_ylim([0.,ymax])
        ax.tick_params(direction="out",labelsize=24)
        xtick_intvl=float(self.bin_max-self.bin_min)/float(nxticks)
        ytick_intvl=float(ymax)/float(nyticks)
        xtick_vals=[self.bin_min+(float(i)*xtick_intvl) for i in range(nxticks+1)]
        if xtick_intvl.is_integer(): xtick_vals = [int(xtick_val) for xtick_val in xtick_vals]
        ytick_vals=[0.+(float(i)*ytick_intvl) for i in range(nyticks+1)]
        ax.set_xticks(xtick_vals)
        ax.set_yticks(ytick_vals)
        xticklabels=["$"+str(xtick_val)+"$" for xtick_val in xtick_vals]
        yticklabels=["$"+str(ytick_val)+"$" for ytick_val in ytick_vals]
        ax.set_xticklabels(xticklabels)
        ax.set_yticklabels(yticklabels)
        plt.tight_layout()
        plt.savefig("fpt_distribn."+figfmt,format=figfmt,bbox_inches="tight")
        plt.show()

    def calc_mfpt(self):
        if not self.logvals: return np.sum(self.vals)/float(self.ntpaths)
        else: return np.sum([10**val for val in self.vals])/float(self.ntpaths)

    def calc_stderr(self,mfpt):
        stderr=0.
        for val in self.vals:
            if not self.logvals: stderr+=(val-mfpt)**2
            else: stderr+=(10**val-mfpt)**2
        return sqrt((1./float(self.ntpaths-1))*stderr)/sqrt(float(self.ntpaths))

if __name__=="__main__":
    ### CHOOSE PARAMS ###

    # statistic to analyse
    # 0=dynamical activity, 1=time, 2=path prob, 3=entropy flow
    stat=1

    #binning params
    nbins=75
    binw=0.2
    bin_min=0.
    binall=False # enforce that all values must be encompassed in the bin range
    logvals=True # take log_10 of values
    # plot params
    nxticks=15 # no. of ticks on x axis
    nyticks=12 # no. of ticks on y axis
    ymax=0.12 # max value for y (prob) axis

    # run
    calc_hist_obj=Analyse_tp_distribns(stat,nbins,binw,bin_min,binall,logvals)
    hist_arr = calc_hist_obj.get_hist_arr()
    print hist_arr
    print "total number of observed A<-B transition paths: ", calc_hist_obj.ntpaths
    print "total number of binned A<-B transition paths: ", np.sum(hist_arr)
    mfpt = calc_hist_obj.calc_mfpt()
    print "MFPT: ", mfpt
    print "standard error: ", calc_hist_obj.calc_stderr(mfpt)
    # plot
#    calc_hist_obj.plot_hist(hist_arr,nxticks,nyticks,ymax)
