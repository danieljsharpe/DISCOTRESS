'''
Python script to process data from DISCOTRESS: plot selected trajectories and analyse the time-dependent probability distributions

Daniel J. Sharpe
Mar 2020
'''

import numpy as np
import matplotlib.pyplot as plt
import os

class Analyse_trajs(object):

    def __init__(self,first_wid,nwalkers,nsupstates,tmax,tintvl,logtime,dir_rootname,nruns):
        self.first_wid=first_wid
        self.nwalkers=nwalkers
        self.nsupstates=nsupstates
        self.tmax=tmax
        self.tintvl=tintvl
        self.logtime=logtime
        assert (tmax/tintvl).is_integer()
        self.tbinvals=[0.+(float(i)*self.tintvl) for i in range(int(tmax/tintvl))]
        self.probdistribn=np.zeros((self.nsupstates,len(self.tbinvals)),dtype=int)
        self.dir_rootname=dir_rootname
        self.nruns=nruns
        self.cwd=os.getcwd() # current working directory
        # colours etc for plotting
        self.pltcolors=["royalblue","salmon","springgreen","gold","fuchsia","aqua","blueviolet", \
                        "seagreen","orange","lightgray"]
        self.pltlabels=["$\mathrm{U}$","$\mathrm{I1}$","$\mathrm{I2}$","$\mathrm{F}$","$\mathrm{V}$", \
                        "$\mathrm{VI}$","$\mathrm{VII}$","$\mathrm{VIII}$","$\mathrm{XI}$","$\mathrm{X}$"]

    # read a single column of a file
    @staticmethod
    def read_single_col(fname,dirname,colidx=0,fmtfunc=int):
        vals=[]
        with open(dirname+"/"+fname,"r") as scfile:
            for line in scfile.readlines():
                vals.append(fmtfunc(line.split()[colidx]))
        return vals

    # calculate transition path statistics averaged over several independent runs
    def calc_tp_stats(self,n_tpbins):
        qi_vals= np.full(n_tpbins,float("nan"),dtype=float) # transition path committor probability values
        ri_vals = np.full(n_tpbins,0.,dtype=float) # transition path probability density values
        for nrun in range(self.nruns):
            if self.nruns==1 and self.dir_rootname is None:
                tp_stats_dir=self.cwd
            else:
                tp_stats_dir=self.cwd+"/"+self.dir_rootname+"_"+str(nrun+1)
            qi_vals_iter=Analyse_trajs.read_single_col("tp_stats.dat",tp_stats_dir,colidx=4,fmtfunc=float)
            ri_vals_iter=Analyse_trajs.read_single_col("tp_stats.dat",tp_stats_dir,colidx=3,fmtfunc=float)
            # accumulate values
            qi_vals = np.nansum(np.stack((qi_vals,qi_vals_iter)),axis=0)
            ri_vals += ri_vals_iter
        # calculate averaged valyes
        qi_vals *= 1./float(self.nruns)
        ri_vals *= 1./float(self.nruns)
        for i in range(n_tpbins):
            if ri_vals[i]==0.: qi_vals[i]=float("nan") # no data on this state for calculation of committor
        # write averaged values to file
        with open("tp_committors.dat","w") as qi_vals_f:
            for qi in qi_vals: qi_vals_f.write("%.8f\n" % qi)
        with open("tp_densities.dat","w") as ri_vals_f:
            for ri in ri_vals: ri_vals_f.write("%.8f\n" % ri)

    # calculate time-dependent probability distribution p(t) for the superstates
    def calc_timedepdistribn(self):
        superstates=Analyse_trajs.read_single_col("superstates.dat",self.cwd) # superstates to which nodes belong
        assert max(superstates)==self.nsupstates-1
        for nrun in range(self.nruns):
            for wid in range(self.first_wid,self.first_wid+self.nwalkers):
                if self.nruns==1 and self.dir_rootname is None:
                    walker_fname="walker.0."+str(wid)+".dat"
                else:
                    walker_fname="/"+self.dir_rootname+"_"+str(nrun+1)+"/walker.0."+str(wid)+".dat"
                curr_tbin_idx=0
                last_supstate=None
                walker_f=open(self.cwd+walker_fname,"r")
                print "reading info for walker: ", wid
                for line in walker_f.readlines():
                    node_id=int(line.split()[0])
                    supstate=superstates[node_id-1]
                    t=float(line.split()[3]) # time (from walker file)
                    if self.logtime:
                        if t<1.: t=1. # if logtime, round values <1. to 1., for easier interpretation
                        t=np.log10(t)
#                    print "time: ", t, " bin_idx: ", curr_tbin_idx, " bin_val: ", self.tbinvals[curr_tbin_idx]
                    if self.tbinvals[curr_tbin_idx]>t: # not yet reached next bin increment - only take the first time the bin is seen
                        last_supstate=supstate
                        continue # skip to next recorded node
                    while self.tbinvals[curr_tbin_idx]<t:
                        self.probdistribn[last_supstate,curr_tbin_idx]+=1 # increment count for previous bin
                        curr_tbin_idx+=1
                    # increment count for this bin
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
            assert abs(np.sum(self.probdistribn[:,i])-1.)<1.E-08
#        for i in range(self.nsupstates):
#            print self.probdistribn[i,:]

    # plot trajectories against a given order parameter
    def plot_trajs(self,nxticks,nyticks,op_fname,op_minval,op_maxval,figfmt="pdf"):
        op_vals=Analyse_trajs.read_single_col(op_fname,self.cwd,colidx=0,fmtfunc=float) # read order param values from file
        plt.figure(figsize=(12,8.5)) # figure size in inches
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
        if self.logtime: plt.xlabel("$\log_{10}(t\ /\ \mathrm{s})$",fontsize=42)
        else: plt.xlabel("$t$",fontsize=42)
        plt.ylabel("$\mathrm{Potential\ energy}\ /\ \mathrm{kcal\,mol}^{-1}$",fontsize=42)
        ax=plt.gca()
        ax.set_xlim([0.,self.tmax])
        ax.set_ylim([op_minval,op_maxval])
        ax.tick_params(direction="out",labelsize=24)
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
        plt.tight_layout()
        plt.savefig("trajs_rep."+figfmt,format=figfmt,bbox_inches="tight")
        plt.show()

    # plot time-dependent occupation probability distribution for superstates
    def plot_probdistribn(self,nxticks,figfmt="pdf",supstate_b=0,supstate_a=None):
        if supstate_a is None: supstate_a=self.nsupstates-1
        plt.figure(figsize=(12,8.5)) # figure size in inches
        for i in range(self.nsupstates):
            plt.plot(self.tbinvals,self.probdistribn[i,:],linewidth=5,color=self.pltcolors[i],label=self.pltlabels[i])
            if i==supstate_b or i==supstate_a: # shade region under curves corresponding to initial and final states
                plt.fill_between(self.tbinvals,self.probdistribn[i,:],facecolor=self.pltcolors[i],alpha=0.5)
        if self.logtime:
            plt.xlabel("$\log_{10}(t\ /\ \mathrm{s})$",fontsize=42)
            plt.ylabel("$\mathrm{Occupation\ probability}\ p(\log_{10} t)$",fontsize=42) # note "\" before space char
        else:
            plt.xlabel("$t$",fontsize=42)
            plt.ylabel("$\mathrm{Occupation\ probability}\ p(t)$",fontsize=42)
        plt.legend(loc="center right",fontsize=24)
        ax=plt.gca()
        ax.set_xlim([0.,self.tmax])
        ax.set_ylim([0.,1.])
        ax.tick_params(direction="out",labelsize=24)
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
        plt.tight_layout()
        plt.savefig("tprobdistribn."+figfmt,format=figfmt,bbox_inches="tight")
        plt.show()

if __name__=="__main__":
    ### CHOOSE PARAMS ###
    first_wid=0 # ID of first walker used in analysis
    nwalkers=3 # total number of walkers used in analysis (assumed to be same for each of the nruns runs)
    nsupstates=4 # number of superstates (in file "superstates.dat", indexed from 0) to which nodes belong, used to calc p(t) for the superstates
    tmax=15. # maximum time (log_{10} if logtime). The minimum time is assumed to be zero (even on a log scale!)
    tintvl=0.1 # time interval for binning occupation probabilities
    logtime=True # use log_{10} of time. Note that values t<1. are rounded up to 1
    dir_rootname=None # if not None, then assume traj info is stored in nruns dirs of the current dir, dir_rootname_x, x=1,...,nruns
    nruns=1 # number of independent kMC runs (each assumed to be of length nwalkers)
    n_tpbins=68777 # number of bins in calculation of transition path statistics
    # plot params
    nxticks=15
    op_fname="energies.dat" # order parameter (single-col) file name, for plotting individual trajectories
    nyticks=10 # only used if plotting indiv trajectories
    op_minval, op_maxval = -455., -405.

    ### RUN ###
    analyse_trajs_obj=Analyse_trajs(first_wid,nwalkers,nsupstates,tmax,tintvl,logtime,dir_rootname,nruns)
#    analyse_trajs_obj.calc_tp_stats(n_tpbins)
    # calc and plot time-dependent probability distribution for superstates
#    analyse_trajs_obj.calc_timedepdistribn()
#    analyse_trajs_obj.plot_probdistribn(nxticks)
    # plot representative trajectories
    analyse_trajs_obj.plot_trajs(nxticks,nyticks,op_fname,op_minval,op_maxval)
