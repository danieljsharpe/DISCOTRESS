'''
Python script to process data from DISCOTRESS: plot selected trajectories and analyse the time-dependent probability distributions

Daniel J. Sharpe
Mar 2020
'''

import numpy as np
import matplotlib.pyplot as plt

class Analyse_trajs(object):

    def __init__(self,first_wid,nwalkers,nsupstates,tmax,tintvl,logtime):
        self.first_wid=first_wid
        self.nwalkers=nwalkers
        self.nsupstates=nsupstates
        self.tmax=tmax
        self.tintvl=tintvl
        self.logtime=logtime
        assert (tmax/tintvl).is_integer()
        self.tbinvals=[0.+(float(i)*self.tintvl) for i in range(int(tmax/tintvl))]
        self.probdistribn=np.zeros((self.nsupstates,len(self.tbinvals)),dtype=int)
        # colours etc for plotting
        self.pltcolors=["royalblue","salmon","springgreen","gold","fuchsia","aqua","blueviolet", \
                        "seagreen","orange","lightgray"]
        self.pltlabels=["$\mathrm{I}$","$\mathrm{II}$","$\mathrm{III}$","$\mathrm{IV}$","$\mathrm{V}$", \
                        "$\mathrm{VI}$","$\mathrm{VII}$","$\mathrm{VIII}$","$\mathrm{XI}$","$\mathrm{X}$"]

    # read a single column of a file
    @staticmethod
    def read_single_col(fname,colidx=0,fmtfunc=int):
        vals=[]
        with open(fname,"r") as scfile:
            for line in scfile.readlines():
                vals.append(fmtfunc(line.split()[colidx]))
        return vals

    # calculate time-dependent probability distribution p(t) for the superstates
    def calc_timedepdistribn(self):
        superstates=Analyse_trajs.read_single_col("superstates.dat") # superstates to which nodes belong
        assert max(superstates)==self.nsupstates-1
        for wid in range(self.first_wid,self.first_wid+self.nwalkers):
            curr_tbin_idx=0
            last_supstate=None
            walker_f=open("walker.0."+str(wid)+".dat","r")
            print "reading info for walker: ", wid
            for line in walker_f.readlines():
                node_id=int(line.split()[0])
                supstate=superstates[node_id-1]
                t=float(line.split()[3]) # time (from walker file)
                if self.logtime:
                    if t<1.: t=1. # if logtime, round values <1. to 1., for easier interpretation
                    t=np.log10(t)
#                print "time: ", t, " bin_idx: ", curr_tbin_idx, " bin_val: ", self.tbinvals[curr_tbin_idx]
                if self.tbinvals[curr_tbin_idx]>t: # not yet reached next bin increment - only take the first time the bin is seen
#                    print "  skipping to next recorded node"
                    last_supstate=supstate
                    continue
                while self.tbinvals[curr_tbin_idx]<t:
#                    print "    incrementing count for previous bin: ", self.tbinvals[curr_tbin_idx], " superstate: ", last_supstate
                    self.probdistribn[last_supstate,curr_tbin_idx]+=1
                    curr_tbin_idx+=1
#                print "    incrementing count for this bin: ", self.tbinvals[curr_tbin_idx], " superstate: ", supstate
                self.probdistribn[supstate,curr_tbin_idx]+=1
                curr_tbin_idx+=1
                last_supstate=supstate
            walker_f.close()
            # the final set is absorbing, pad the remainder of the probdistribn array
            while curr_tbin_idx<len(self.tbinvals):
                self.probdistribn[last_supstate,curr_tbin_idx]+=1
                curr_tbin_idx+=1
        # convert counts to probabilities
        self.probdistribn=self.probdistribn.astype(float)
        for i in range(len(self.tbinvals)-1):
            self.probdistribn[:,i] *= 1./np.sum(self.probdistribn[:,i])
#            print np.sum(self.probdistribn[:,i])
            assert abs(np.sum(self.probdistribn[:,i])-1.)<1.E-08
#        for i in range(self.nsupstates):
#            print self.probdistribn[i,:]

    # plot trajectories against a given order parameter
    def plot_trajs(self,nxticks,nyticks,op_fname,op_minval,op_maxval,figfmt="pdf"):
        op_vals=Analyse_trajs.read_single_col(op_fname,colidx=0,fmtfunc=float) # read order param values from file
        plt.figure(figsize=(10,6)) # figure size in inches
        for i, wid in enumerate(range(self.first_wid,self.first_wid+self.nwalkers)):
            walkertdata=[] # time data for walker
            walkeropdata=[] # order param data for walker
            walker_f=open("walker.0."+str(wid)+".dat","r")
            for line in walker_f.readlines():
                node_id=int(line.split()[0])
                t=float(line.split()[3])
                if self.logtime:
                    if t<1.: t=1. # if logtime, round values <1. to 1., for easier interpretation
                    t=np.log10(t)
                op_val=op_vals[node_id-1]
                walkertdata.append(t)
                walkeropdata.append(op_val)
            walker_f.close()
            plt.plot(walkertdata,walkeropdata,linewidth=5,color=self.pltcolors[i])
        if self.logtime: plt.xlabel("$\log_{10} t$",fontsize=24)
        else: plt.xlabel("$t$",fontsize=24)
        plt.ylabel("$\mathrm{Energy}\ \mathrm{kcal\,mol}^{-1}$",fontsize=24)
        ax=plt.gca()
        ax.set_xlim([0.,self.tmax])
        ax.set_ylim([op_minval,op_maxval])
        ax.tick_params(direction="out",labelsize=18)
        xtick_intvl=self.tmax/float(nxticks)
        xtick_vals=[0.+(float(i)*xtick_intvl) for i in range(nxticks+1)]
        if xtick_intvl.is_integer(): xtick_vals = [int(xtick_val) for xtick_val in xtick_vals]
        ytick_intvl=(op_maxval-op_minval)/float(nyticks)
        ytick_vals=[op_minval+(float(i)*ytick_intvl) for i in range(nyticks+1)]
        if ytick_intvl.is_integer(): ytick_vals = [int(ytick_val) for ytick_val in ytick_vals]
        ax.set_xticks(xtick_vals)
        ax.set_yticks(ytick_vals)
        xticklabels=["$"+str(xtick_val)+"$" for xtick_val in xtick_vals]
        yticklabels=["$"+str(ytick_val)+"$" for ytick_val in ytick_vals]
        ax.set_xticklabels(xticklabels)
        ax.set_yticklabels(yticklabels)
        plt.savefig("trajs_rep."+figfmt,format=figfmt)
        plt.show()

    # plot time-dependent occupation probability distribution for superstates
    def plot_probdistribn(self,nxticks,figfmt="pdf"):
        plt.figure(figsize=(10,6)) # figure size in inches
        for i in range(self.nsupstates):
            plt.plot(self.tbinvals,self.probdistribn[i,:],linewidth=5,color=self.pltcolors[i],label=self.pltlabels[i])
        if self.logtime:
            plt.xlabel("$\log_{10} t$",fontsize=24)
            plt.ylabel("$\mathrm{Occupation\ probability}\ p(\log_{10} t)$",fontsize=24) # note "\" before space char
        else:
            plt.xlabel("$t$",fontsize=24)
            plt.ylabel("$\mathrm{Occupation\ probability}\ p(t)$",fontsize=24)
        plt.legend(loc="upper right")
        ax=plt.gca()
        ax.set_xlim([0.,self.tmax])
        ax.set_ylim([0.,1.])
        ax.tick_params(direction="out",labelsize=18)
        xtick_intvl=self.tmax/float(nxticks)
        xtick_vals=[0.+(float(i)*xtick_intvl) for i in range(nxticks+1)]
        if xtick_intvl.is_integer(): xtick_vals = [int(xtick_val) for xtick_val in xtick_vals]
        ytick_vals=[0.+(float(i)*0.1) for i in range(11)]
        ax.set_xticks(xtick_vals)
        ax.set_yticks(ytick_vals)
        xticklabels=["$"+str(xtick_val)+"$" for xtick_val in xtick_vals]
        yticklabels=["$"+str(ytick_val)+"$" for ytick_val in ytick_vals]
        ax.set_xticklabels(xticklabels)
        ax.set_yticklabels(yticklabels)
        plt.savefig("tprobdistribn."+figfmt,format=figfmt)
        plt.show()

if __name__=="__main__":
    ### CHOOSE PARAMS ###
    first_wid=0 # ID of first walker used in analysis
    nwalkers=10 # total number of walkers used in analysis
    nsupstates=10 # number of superstates (in file "superstates.dat", indexed from 0) to which nodes belong, used to calc p(t) for the superstates
    tmax=15. # maximum time (log_{10} if logtime). The minimum time is assumed to be zero (even on a log scale!)
    tintvl=0.1 # time interval for binning occupation probabilities
    logtime=True # use log_{10} of time. Note that values t<1. are rounded up to 1
    # plot params
    nxticks=15
    op_fname="energies.dat" # order parameter (single-col) file name, for plotting individual trajectories
    nyticks=10 # only used if plotting indiv trajectories
    op_minval, op_maxval = -455., -405.

    # run
    analyse_trajs_obj=Analyse_trajs(first_wid,nwalkers,nsupstates,tmax,tintvl,logtime)
#    analyse_trajs_obj.calc_timedepdistribn()
#    analyse_trajs_obj.plot_probdistribn(nxticks)
    analyse_trajs_obj.plot_trajs(nxticks,nyticks,op_fname,op_minval,op_maxval)
