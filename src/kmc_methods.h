/*
Classes and functions for handling enhanced kinetic Monte Carlo simulations and propagating the trajectories
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

using namespace std;

/* data structure for a single trajectory (walker) on the transition network */
struct Walker {

    public:

    explicit Walker()=default;
    ~Walker();
    void dump_walker_info(int,bool,bool,bool=true); // write path quantities to files
    void reset_walker_info();

    int walker_id; // ID of walker in set of trajectories
    int comm_curr, comm_prev; // for WE-kMC
    unsigned long long int k; // path activity
    bool active; // walker is currently a member of the set of active trajectories being propagated
    bool accumprobs=false; // indicates if the Network on which the Walker is active has accumulated transition probs
    long double p; // (log) path probability
    long double t; // path time (stochastically sampled)
    long double s; // entropy flow along path
    const Node *curr_node; // pointer to node currently occupied by the walker
    vector<bool> visited;  // element is true when the corresponding bin has been visited along the trajectory
};

/* abstract class containing functions to handle enhanced sampling kMC methods */
class KMC_Enhanced_Methods {

    protected:

    int maxn_abpaths;           // algorithm terminates when this number of A-B paths have been successfully sampled
    int n_ab;                   // number of successfully simulated A<-B transition paths
    int n_traj;                 // total number of B<-B or A<-B paths simulated
    double tintvl;              // time interval for dumping trajectory data
    vector<int> ab_successes;   // vector of counts of bin appearances along A<-B transition paths
    vector<int> ab_failures;    // vector of counts of bin appearances along B<-B unreactive paths
    vector<double> tp_densities; // probability that a bin is visited along an A<-B transition path
    vector<double> committors;  // forward (A<-B) committor probabilities for bins
    int maxit;                  // another termination condition; the maximum number of iterations of the enhanced kMC method
    int seed;                   // seed for random number generator
    bool debug;                 // debug printing on/off
    void (*kmc_func)(Walker&);  // function pointer to kMC algorithm for propagating the trajectory   

    public:

    KMC_Enhanced_Methods();
    virtual ~KMC_Enhanced_Methods();
    virtual void run_enhanced_kmc(const Network&)=0; // pure virtual function
    Node *get_initial_node(const Network&, Walker&); // sample an initial node
    void set_standard_kmc(void(*)(Walker&)); // function to set the kmc_std_method
    static vector<int> find_comm_onthefly(const Network&,const Node*,double,int); // find a community on-the-fly based on max allowed rate and size
    void update_tp_stats(Walker&,bool,bool); // update the transition path statistics, depends on if the path is a transn path or is unreactive
    void calc_tp_stats(int);    // calculate the transition path statistics from the observed counts
    void write_tp_stats(int);   // write transition path statistics to file

    template <typename T>
    static void write_vec(const vector<T>& vec, string fname) {
        ofstream vec_f; vec_f.open(fname);
        if (typeid(T)==typeid(double) || typeid(T)==typeid(long double)) {
            vec_f.precision(30); vec_f.setf(ios::scientific,ios::floatfield); }
        for (const T elem: vec) vec_f << setw(30) << elem << endl;
    }
};

/* Standard kMC, simply propagates the dynamics of a single trajectory using the chosen standard method */
class STD_KMC : public KMC_Enhanced_Methods {

    private:

    Walker walker={walker_id:0,comm_curr:0,comm_prev:0,k:0,active:true,accumprobs:false,\
                   p:-numeric_limits<long double>::infinity(),t:0.,s:0.}; // method uses only a single walker

    public:

    STD_KMC(const Network&,int,int,double,int);
    ~STD_KMC();
    void run_enhanced_kmc(const Network&);
};

/* Weighted ensemble kMC */
class WE_KMC : public KMC_Enhanced_Methods {

    private:

    vector<Walker> walkers; // list of active trajectories (walkers) on the network
    int nwalkers;
    double tau; // time interval between checking communities and resampling trajectories
    bool adaptivecomms;
    double adaptminrate;

    void we_resampling();

    public:

    WE_KMC(const Network&,int,int,long double,double,bool,int,bool);
    ~WE_KMC();
    void run_enhanced_kmc(const Network&);
};

/* Kinetic path sampling (kPS)
   Note that kPS only uses a single walker due to its memory requirements.
   Also note that the number of kMC self-hops/transition hops are stored in the Node and Edge data structures,
   respectively, of the subnetwork stored via the ktn_kps pointer. */
class KPS : public KMC_Enhanced_Methods {

    private:

    Network *ktn_kps=nullptr; // pointer to the subnetwork of the TN that kPS internally uses and transforms
    Network *ktn_kps_orig=nullptr; // pointer to the original subnetwork of the TN
    Network *ktn_kps_gt=nullptr; // pointer to the graph-transformed subnetwork (used if recycling GT of a basin)
    Network *ktn_l=nullptr, *ktn_u=nullptr; // pointers to arrays used in LU-style decomposition of transition matrix
    Walker walker={walker_id:0,comm_curr:0,comm_prev:0,k:0,active:true,accumprobs:false,\
                   p:-numeric_limits<long double>::infinity(),t:0.,s:0.};
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
        // NB these pointers point to nodes in the original network, passed as the arg to run_enhanced_kmc()
    bool adaptivecomms; // comunities are defined adaptively (or else are set prior to the simulation)
    double adaptminrate; // maximum allowed rate in finding a community on-the-fly
    bool pfold;      // calculate the committor functions instead of performing a kPS simulation
    int kpskmcsteps; // number of kMC steps to run after each kPS trapping basin escape trajectory sampled
    double next_tintvl; // next time interval for dumping trajectory data

    void setup_basin_sets(const Network&,bool);
    long double iterative_reverse_randomisation();
    Node *sample_absorbing_node();
    void graph_transformation(const Network&);
    void gt_iteration(Node*);
    vector<pair<Node*,Edge*>> undo_gt_iteration(Node*);
    void update_path_quantities(long double,const Node*);
    Network *get_subnetwork(const Network&,bool);
    void calc_pfold(const Network&);

    public:

    KPS(const Network&,int,int,int,long double,double,int,bool,double,bool,int,bool);
    ~KPS();
    void run_enhanced_kmc(const Network&);
    static long double calc_gt_factor(Node*);
    static void reset_kmc_hop_counts(Network&);
    static long double gamma_distribn(unsigned long long int,long double,int);
    static unsigned long long int binomial_distribn(unsigned long long int,long double,int);
    static unsigned long long int negbinomial_distribn(unsigned long long int,long double,int);
    static long double exp_distribn(long double,int);
    static void test_ktn(const Network&);
};

/* Forward flux sampling kMC */
class FFS_KMC : public KMC_Enhanced_Methods {

    private:

    vector<Walker> walkers;

    public:

    FFS_KMC(const Network&);
    ~FFS_KMC();
    void run_enhanced_kmc(const Network&);
};

/* Monte Carlo with absorbing Markov chains (MCAMC) */
class MCAMC : public KMC_Enhanced_Methods {

    public:

    MCAMC(const Network&);
    ~MCAMC();
    void run_enhanced_kmc(const Network&);
};

/* non-equilibrium umbrella sampling kMC */
class NEUS_KMC : public KMC_Enhanced_Methods {

    public:

    NEUS_KMC(const Network&);
    ~NEUS_KMC();
    void run_enhanced_kmc(const Network&);
};

/* milestoning kMC */
class MILES_KMC : public KMC_Enhanced_Methods {

    public:

    MILES_KMC(const Network&);
    ~MILES_KMC();
    void run_enhanced_kmc(const Network&);
};

/* transition path sampling kMC */
class TPS_KMC : public KMC_Enhanced_Methods {

    public:

    TPS_KMC(const Network&);
    ~TPS_KMC();
    void run_enhanced_kmc(const Network&);
};

/* class containing functions to propagate KMC trajectories */
class KMC_Standard_Methods {

    public:

    KMC_Standard_Methods();
    ~KMC_Standard_Methods();
    static void bkl(Walker&); // rejection-free algorithm of Bortz, Kalos and Lebowitz (aka n-fold way algorithm)
    static void rejection_kmc(Walker&); // kMC algorithm where some moves are rejected
    static void leapfrog(Walker&); // leapfrog algorithm of Trygubenko & Wales
    static long double rand_unif_met(int=19); // draw random number between 0 and 1
    static vector<double> calc_committors(const Network&,int); // calculate committor functions from counts in Node structures of Network obj
    static vector<double> calc_tp_density(const Network&,int); // calculate transn path density from counts in Node structures of Network obj
};

#endif
