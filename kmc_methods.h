/*
Classes and functions for handling enhanced kinetic Monte Carlo simulations and propagating the trajectories
*/

#ifndef __KMC_METHODS_H_INCLUDED__
#define __KMC_METHODS_H_INCLUDED__

#include "ktn.h"
#include <map>

using namespace std;

/* data structure for a single trajectory (walker) on the transition network */
struct Walker {

    public:

    explicit Walker()=default;
    ~Walker();
    int walker_id; // ID of walker in set of trajectories
    int bin_curr, bin_prev; // for WE-kMC
    int k; // path activity
    bool active; // walker is currently a member of the set of active trajectories being propagated
    double p; // path probability
    double t; // path time (stochastically sampled)
    double s; // entropy flow along path
};

/* abstract class containing functions to handle enhanced sampling kMC methods */
class KMC_Enhanced_Methods {

    protected:

    int n_abpaths;  // algorithm terminates when this number of A-B paths have been successfully sampled
    int seed; // seed for random number generator

    public:

    KMC_Enhanced_Methods();
    virtual ~KMC_Enhanced_Methods();
    virtual void run_enhanced_kmc(const Network&)=0; // pure virtual function
    void find_bin_onthefly();
};

/* Weighted ensemble kMC */
class WE_KMC : public KMC_Enhanced_Methods {

    private:

    vector<Walker> walkers; // list of active trajectories (walkers) on the network

    public:

    WE_KMC(const Network&,double,int,bool);
    ~WE_KMC();
    void run_enhanced_kmc(const Network&);

    int nbins;
    int nwalkers;
    double tau; // time interval between checking bins
};

/* Kinetic path sampling (kPS) */
class KPS : public KMC_Enhanced_Methods {

    private:

    Network *ktn_kps=nullptr; // pointer to the subnetwork of the TN that kPS internally uses and transforms
    Network *ktn_kps_orig=nullptr; // pointer to the original subnetwork of the TN
    Network *ktn_l=nullptr, *ktn_u=nullptr; // pointers to arrays used in LU-style decomposition of transition matrix
    Walker walker={walker_id:0,bin_curr:0,bin_prev:0,k:0,active:true,p:0.,t:0.,s:0.}; // kPS only uses a single walker due to its memory requirements
    // array<array<int>> H; // hopping matrix (sparse format)
    vector<int> h;  // flicker vector
    vector<int> basin_ids; // used to indicate the set to which each node belongs for the current kPS iteration
        // (eliminated=1, transient noneliminated=2, absorbing boundary=3, absorbing nonboundary=0)
    vector<Node*> eliminated_nodes; // vector of eliminated nodes (in order)
    map<int,int> nodemap; // map of ID's for full network to subnetwork
    double tau;     // lag time at which transition matrix is evaluated
    int nelim;      // maximum number of nodes of a trapping basin to be eliminated
    int nbins;      // number of bins specified (if not adaptivebins)
    int N_c;        // number of nodes connected to the eliminated states of the current trapping basin
    int N, N_B;     // number of eliminated nodes / total number of nodes for the currently active trapping basin
    int N_e;        // number of edges in the subnetwork
    const Node *alpha=nullptr, *epsilon=nullptr; // final and initial microstates of current escape trajectory
        // NB these pointers point to nodes in the original network, passed as the arg to run_enhanced_kmc()
    bool adaptivebins; // bins are defined adaptively (or else are set prior to the simulation)
    bool initcond;  // an initial condition to sample the starting node has been specified (Y/N)
    int kpskmcsteps; // number of kMC steps to run after each kPS trapping basin escape trajectory sampled
    int n_kpsmaxit; // algorithm terminates when maximum number of kPS trapping basin escapes have been simulated
    bool debug;

    void setup_basin_sets(const Network&);
    double iterative_reverse_randomisation();
    Node *sample_absorbing_node();
    void graph_transformation(const Network&);
    void gt_iteration(Node*);
    void undo_gt_iteration(Node*);
    void update_path_quantities(double);

    Network *get_subnetwork(const Network&);

    public:

    KPS(const Network&,int,int,int,double,int,int,bool,bool,int,bool);
    ~KPS();
    void run_enhanced_kmc(const Network&);
    static double calc_gt_factor(Node*);
    static double gamma_distribn(int,double,int);
    static int binomial_distribn(int,double,int);
    static int negbinomial_distribn(int,double,int);
    static double exp_distribn(double,int);
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

/* accelerated superbasin kMC */
class AS_KMC : public KMC_Enhanced_Methods {

    public:

    AS_KMC(const Network&);
    ~AS_KMC();
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
    static void leapfrog(Walker&); // leapfrog algorithm of Wales
    static void rejection_kmc(Walker&); // kMC algorithm where some moves are rejected
    vector<int> setup_move_probs(const Network&); // calculate the lists of move probabilities for each node
    static double rand_unif_met(int);
};

#endif
