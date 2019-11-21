/*
Classes and functions for handling enhanced kinetic Monte Carlo simulations and propagating the trajectories
*/

#ifndef __KMC_METHODS_H_INCLUDED__
#define __KMC_METHODS_H_INCLUDED__

#include "ktn.h"

using namespace std;

/* data structure for a single trajectory (walker) on the transition network */
struct Walker {

    public:

    Walker();
    ~Walker();
    int walker_id;
    int bin_curr, bin_prev; // for WE-kMC
    bool active; // walker is currently a member of the set of active trajectories being propagated
};

/* abstract class containing functions to handle enhanced sampling kMC methods */
class KMC_Enhanced_Methods {

    private:

    int n_abpaths;  // algorithm terminates when this number of A-B paths have been successfully sampled

    public:

    KMC_Enhanced_Methods();
    ~KMC_Enhanced_Methods();
    virtual void run_enhanced_kmc(Network&)=0; // pure virtual function
    void find_bin_onthefly();
};

/* Weighted ensemble kMC */
class WE_KMC : public KMC_Enhanced_Methods {

    private:

    vector<Walker> walkers; // list of active trajectories (walkers) on the network

    public:

    WE_KMC(Network&,double,int,bool);
    ~WE_KMC();
    void run_enhanced_kmc(Network&);

    int nbins;
    int nwalkers;
    double tau; // time interval between checking bins
};

/* Kinetic path sampling (kPS) */
class KPS : public KMC_Enhanced_Methods {

    private:

    // array<array<int>> H; // hopping matrix (sparse format)
    vector<int> h;  // flicker vector
    vector<int> basin_ids; // used to indicate the set to which each node belongs for the current kPS iteration
        // (eliminated=1, transient noneliminated boundary=2, transient noneliminated nonboundary=3,
        //  absorbing boundary=4, absorbing nonboundary=0)
    double tau;     // lag time at which transition matrix is evaluated
    int nelim;      // maximum number of nodes of a trapping basin to be eliminated
    int nbins;      // number of bins specified (if not adaptivebins)
    int N_c;        // number of nodes connected to the eliminated states of the current trapping basin
    int N;          // number of eliminated nodes for the currently active trapping basin
    Node *alpha=nullptr, *epsilon=nullptr; // final and initial microstates of current escape trajectory
    bool adaptivebins; // bins are defined adaptively (or else are set prior to the simulation)
    bool initcond;  // an initial condition to sample the starting node has been specified (Y/N)
    int kpskmcsteps; // number of kMC steps to run after each kPS trapping basin escape trajectory sampled
    int n_kpsmaxit; // algorithm terminates when maximum number of kPS trapping basin escapes have been simulated

    void setup_basin_sets(Network&);
    double iterative_reverse_randomisation();
    Node *sample_absorbing_node(Network&);
    void graph_transformation(Network&);

    public:

    KPS(Network&,int,int,int,double,int,int,bool,bool);
    ~KPS();
    void run_enhanced_kmc(Network&);
    static double gamma_distribn(int,double);
    static int binomial_distribn(int,double);
    static int negbinomial_distribn(int,double);
    static double exp_distribn(double);
    void quack();
};

/* Forward flux sampling kMC */
class FFS_KMC : public KMC_Enhanced_Methods {

    private:

    vector<Walker> walkers;

    public:

    FFS_KMC();
    ~FFS_KMC();
    void run_enhanced_kmc(Network&);
};

/* accelerated superbasin kMC */
class AS_KMC : public KMC_Enhanced_Methods {

    public:

    AS_KMC();
    ~AS_KMC();
    void run_enhanced_kmc(Network&);
};

/* non-equilibrium umbrella sampling kMC */
class NEUS_KMC : public KMC_Enhanced_Methods {

    public:

    NEUS_KMC();
    ~NEUS_KMC();
    void run_enhanced_kmc(Network&);
};

/* class containing functions to propagate KMC trajectories */
class KMC_Standard_Methods {

    public:

    KMC_Standard_Methods();
    ~KMC_Standard_Methods();
    static void bkl(Walker&); // rejection-free algorithm of Bortz, Kalos and Lebowitz (aka n-fold way algorithm)
    static void leapfrog(Walker&); // leapfrog algorithm of Wales
    static void rejection_kmc(Walker&); // kMC algorithm where some moves are rejected
    vector<int> setup_move_probs(Network&); // calculate the lists of move probabilities for each node
    static double rand_unif();
};

#endif
