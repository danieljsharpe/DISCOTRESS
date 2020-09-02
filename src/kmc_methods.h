/*
File containing classes and functions for handling enhanced kinetic Monte Carlo simulations and propagating the trajectories

This file is a part of DISCOTRESS, a software package to simulate the dynamics on arbitrary continuous- and discrete-time Markov chains (CTMCs and DTMCs).
Copyright (C) 2020 Daniel J. Sharpe

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef __KMC_METHODS_H_INCLUDED__
#define __KMC_METHODS_H_INCLUDED__

#include "ktn.h"
#include <limits>
#include <utility>
#include <unordered_map>
#include <string>
#include <typeinfo>
#include <iomanip>
#include <fstream>
#include <omp.h>

using namespace std;

class Traj_Method;

/* data structure for a single trajectory (walker) on the transition network */
struct Walker {

    public:

    explicit Walker()=default;
    ~Walker();
    void dump_walker_info(bool,long double,const Node*,bool=false); // write trajectory data to file
    void dump_tp_distribn(); // append transition path quantities to distributions file
    void reset_walker_info();

    int walker_id; // ID of walker in set of trajectories
    int path_no; // the trajectory iteration for this walker ID
    int comm_curr, comm_prev; // for WE-kMC
    unsigned long long int k; // path activity
    bool active; // walker is currently a member of the set of active trajectories being propagated
    bool accumprobs=false; // indicates if the Network on which the Walker is active has accumulated transition probs
    long double p; // (log) path probability
    long double t; // path time (stochastically sampled)
    long double s; // entropy flow along path
    const Node *prev_node, *curr_node; // pointers to nodes previously and currently occupied by the walker
    vector<bool> visited;  // element is true when the corresponding bin has been visited along the trajectory
};

/* abstract class for wrapper (trajectory handling) enhanced sampling methods */
class Wrapper_Method {

    protected:

    int maxn_abpaths;           // algorithm terminates when this number of A<-B transition paths have been successfully sampled
    int n_ab;                   // number of successfully simulated A<-B transition paths
    int n_traj;                 // total number of B<-B or A<-B paths simulated
    double tintvl;              // time interval for dumping trajectory data
    vector<int> ab_successes;   // vector of counts of bin appearances along A<-B transition paths
    vector<int> ab_failures;    // vector of counts of bin appearances along B<-B unreactive paths
    vector<double> tp_densities; // probability that a bin is visited along an A<-B transition path
    vector<double> committors;  // forward (A<-B) committor probabilities for bins
    bool adaptivecomms; // comunities are defined adaptively (or else are set prior to the simulation)
    int maxit;                  // another termination condition; the maximum number of iterations of the enhanced kMC method
    int seed;                   // seed for random number generator
    bool debug;                 // debug printing on/off
    vector<Walker> walkers;     // list of active trajectories (walkers) on the network
    void (*kmc_func)(Walker&);  // function pointer to kMC algorithm for propagating the trajectory   

    public:

    Wrapper_Method();
    virtual ~Wrapper_Method();
    void setup_wrapper_method(int,double,int,int,bool); // set the protected members of the Wrapper_Method class
    virtual void run_enhanced_kmc(const Network&,Traj_Method*)=0; // pure virtual function
    static const Node *get_initial_node(const Network&, Walker&,int); // sample an initial node
    void set_standard_kmc(void(*)(Walker&)); // function to set the kmc_std_method
    static vector<int> find_comm_onthefly(const Network&,const Node*,double,int); // find a community on-the-fly based on max allowed rate and size
    void update_tp_stats(Walker&,bool,bool); // update the transition path statistics, depends on if the path is a transn path or is unreactive
    void calc_tp_stats(int);    // calculate the transition path statistics from the observed counts
    void write_tp_stats(int);   // write transition path statistics to file
    static long double rand_unif_met(int=19); // draw uniform random number between 0 and 1

    template <typename T>
    static void write_vec(const vector<T>& vec, string fname, int precision=30) {
        ofstream vec_f; vec_f.open(fname);
        if (typeid(T)==typeid(double) || typeid(T)==typeid(long double)) {
            vec_f.precision(precision); vec_f.setf(ios::scientific,ios::floatfield); }
        for (const T elem: vec) vec_f << setw(32) << elem << endl;
    }
};

/* class to handle the simulation of many short nonequilibrium trajectories, used to obtain data required for coarse-graining */
class DIMREDN : public Wrapper_Method {

    private:

    vector<int> ntrajsvec; // vector containing number of trajectories to simulate starting from each community in turn
    long double dt;        // length in time of each trajectory

    public:
    
    DIMREDN(const Network&,vector<int>,long double);
    ~DIMREDN();
    void run_enhanced_kmc(const Network&,Traj_Method*);
};

/* no wrapper enhanced sampling class, simply propagates the dynamics of trajectories using the chosen method */
class STD_KMC : public Wrapper_Method {

    private:

    Walker walker={walker_id:0,path_no:0,comm_curr:0,comm_prev:0,k:0,active:true,accumprobs:false,\
                   p:-numeric_limits<long double>::infinity(),t:0.,s:0.}; // method uses only a single walker

    public:

    STD_KMC(const Network&,bool);
    ~STD_KMC();
    void run_enhanced_kmc(const Network&,Traj_Method*);
};

/* Weighted ensemble kMC */
class WE_KMC : public Wrapper_Method {

    private:

    int nwalkers;
    double taure; // time interval between checking communities and resampling trajectories
    double adaptminrate;

    void we_resampling();

    public:

    WE_KMC(const Network&,double,bool);
    ~WE_KMC();
    void run_enhanced_kmc(const Network&,Traj_Method*);
};

/* Forward flux sampling kMC */
class FFS_KMC : public Wrapper_Method {

    private:


    public:

    FFS_KMC(const Network&);
    ~FFS_KMC();
    void run_enhanced_kmc(const Network&,Traj_Method*);
};

/* non-equilibrium umbrella sampling kMC */
class NEUS_KMC : public Wrapper_Method {

    public:

    NEUS_KMC(const Network&);
    ~NEUS_KMC();
    void run_enhanced_kmc(const Network&,Traj_Method*);
};

/* milestoning kMC */
class MILES_KMC : public Wrapper_Method {

    public:

    MILES_KMC(const Network&);
    ~MILES_KMC();
    void run_enhanced_kmc(const Network&,Traj_Method*);
};

/* abstract class for methods to propagate individual trajectories */
class Traj_Method {

    protected:

    bool discretetime;          // transition probabilities represent a discrete-time Markov chain
    double tintvl;              // time interval for dumping trajectory data
    double next_tintvl;         // next time for dumping trajectory data
    bool dumpintvls;            // specifies that trajectory data is to be dumped at the time intervals
    int seed;
    bool debug;

    public:

    Traj_Method();
    virtual ~Traj_Method();
    Traj_Method(const Traj_Method&);
    virtual Traj_Method* clone() {}
    void setup_traj_method(double,bool,bool,int,bool) ; // set the protected members of the Traj_Method class
    void dump_traj(Walker&,bool,bool,long double=numeric_limits<long double>::infinity()); // call function to dump walker info and then update next_tintvl;
    virtual void kmc_iteration(const Network&,Walker&)=0;
    virtual void do_bkl_steps(const Network&,Walker&,long double=numeric_limits<long double>::infinity()) {} // dummy function overridden in KPS and MCAMC to do BKL steps after a basin escape
    virtual void reset_nodeptrs() {} // dummy function overridden in KPS and MCAMC to reset basin and absorbing node pointers when A is hit
    bool statereduction=false;    // purpose of the computation is to perform a state reduction procedure, not a simulation
};

/* rejection-free algorithm of Bortz, Kalos and Lebowitz (aka n-fold way algorithm) */
class BKL : public Traj_Method {

    public:

    BKL(const Network&,bool);
    ~BKL();
    BKL(const BKL&);
    BKL* clone() { return new BKL(*this); } // NB this calls copy constructor for BKL
    void kmc_iteration(const Network&,Walker&);
    static void bkl(Walker&,bool,int);
};

/* kinetic path sampling (kPS)
   Note that the number of kMC self-hops/transition hops are stored in the Node and Edge data structures,
   respectively, of the subnetwork stored via the ktn_kps pointer. */
class KPS : public Traj_Method {

    private:

    Network *ktn_kps=nullptr; // pointer to the subnetwork of the TN that kPS internally uses and transforms
    Network *ktn_kps_orig=nullptr; // pointer to the original subnetwork of the TN
    Network *ktn_kps_gt=nullptr; // pointer to the graph-transformed subnetwork (used if recycling GT of a basin)
    Network *ktn_l=nullptr, *ktn_u=nullptr; // pointers to arrays used in LU-style decomposition of transition matrix
    vector<int> basin_ids; // used to indicate the set to which each node belongs for the current kPS iteration
        // (eliminated=1, transient noneliminated=2, absorbing boundary=3, absorbing nonboundary=0)
    vector<int> eliminated_nodes; // vector of IDs of eliminated nodes (in order)
    unordered_map<int,int> nodemap; // map of node IDs from original network to subnetwork
    long double tau;     // lag time at which transition matrix is evaluated
    int nelim;      // maximum number of nodes of a trapping basin to be eliminated
    int N_c;        // number of nodes connected to the eliminated states of the current trapping basin
    int N, N_B;     // number of eliminated nodes / total number of nodes for the currently active trapping basin
    int N_e;        // number of edges in the subnetwork
    const Node *alpha=nullptr, *epsilon=nullptr; // final and initial microstates of current escape trajectory
        // NB these pointers point to nodes in the original network, passed as the arg to kmc_iteration()
    bool adaptivecomms;
    double adaptminrate; // maximum allowed rate in finding a community on-the-fly
    int kpskmcsteps; // number of kMC steps to run after each kPS trapping basin escape trajectory sampled
    // series of truth values specify which state reduction procedures to perform
    bool committor=false, absorption=false, fundamentalred=false, fundamentalirred=false, mfpt=false, gth=false;
    vector<long double> mfpt_vals; // vector of MFPTs (elem is non-zero for non-absorbing nodes)
    long double mu; // sum of (unnormalised) stationary probabilities in GTH algorithm

    void setup_basin_sets(const Network&,Walker&,bool);
    long double iterative_reverse_randomisation();
    Node *sample_absorbing_node();
    void graph_transformation(const Network&);
    void gt_iteration(Node*);
    vector<pair<Node*,Edge*>> undo_gt_iteration(Node*);
    void update_path_quantities(Walker&,long double,const Node*);
    Network *get_subnetwork(const Network&,bool);
    void do_bkl_steps(const Network&,Walker&,long double=numeric_limits<long double>::infinity());
    void reset_nodeptrs();
    void calc_committor(const Network&);
    void calc_absprobs(); void calc_mfpt(); void calc_gth();
    void calc_fundamentalred(const Network&);
    void write_renormalised_probs(string);
    void rewrite_stat_probs(const Network&);
    static long double committor_boundary_node(const Network&,int,const vector<long double>,int);

    public:

    KPS(const Network&,bool,int,long double,int,bool,double);
    ~KPS();
    KPS(const KPS&);
    KPS* clone() { return new KPS(*this); }
    void set_statereduction_procs(bool,bool,bool,bool,bool,bool);
    void kmc_iteration(const Network&,Walker&);
    static long double calc_gt_factor(Node*);
    static void reset_kmc_hop_counts(Network&);
    static long double gamma_distribn(unsigned long long int,long double,int);
    static unsigned long long int binomial_distribn(unsigned long long int,long double,int);
    static unsigned long long int negbinomial_distribn(unsigned long long int,long double,int);
    static long double exp_distribn(long double,int);
    static void test_ktn(const Network&);
};

/* Monte Carlo with absorbing Markov chains (MCAMC) */
class MCAMC : public Traj_Method {

    private:

    int kpskmcsteps; // number of kMC steps to run after each MCAMC trapping basin escape trajectory sampled
    bool meanrate; // if True, use (approximate) mean rate method, else use (exact) FPTA method
    const Node *alpha=nullptr, *epsilon=nullptr; // final and initial microstates of current escape trajectory

    public:

    MCAMC(const Network&,bool,int,bool);
    ~MCAMC();
    MCAMC(const MCAMC&);
    MCAMC* clone() { return new MCAMC(*this); }
    void kmc_iteration(const Network&,Walker&);
    void do_bkl_steps(const Network&,Walker&,long double=numeric_limits<long double>::infinity());
    void reset_nodeptrs();
};

#endif
