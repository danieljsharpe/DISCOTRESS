/*
Classes and functions for handling enhanced kinetic Monte Carlo simulations and propagating the trajectories
*/

#ifndef __KMC_METHODS_H_INCLUDED__
#define __KMC_METHODS_H_INCLUDED__

#include <vector>
#include <array>
#include <set>

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

/* class containing functions to handle enhanced sampling kMC methods */
class KMC_Enhanced_Methods {

    public:

    KMC_Enhanced_Methods();
    ~KMC_Enhanced_Methods();
    void find_bin_onthefly();
};

/* Weighted ensemble kMC */
class WE_KMC : KMC_Enhanced_Methods {

    private:

    vector<Walker> walkers; // list of active trajectories (walkers) on the network

    public:

    WE_KMC();
    ~WE_KMC();

    int nbins;
    int nwalkers;
    double tau; // time interval between checking bins
};

/* Kinetic path sampling (kPS) */
class KPS : KMC_Enhanced_Methods {

    private:

    // array<array<double>> H; // hopping matrix
    // array<double> h; // flicker vector
    double tau; // lag time at which transition matrix is evaluated

    public:

    KPS();
    ~KPS();
};

/* Forward flux sampling kMC */
class FFS_KMC : KMC_Enhanced_Methods {

    private:

    vector<Walker> walkers;

    public:

    FFS_KMC();
    ~FFS_KMC();
};

/* accelerated superbasin kMC */
class AS_KMC : KMC_Enhanced_Methods {

    public:

    AS_KMC();
    ~AS_KMC();
};

/* non-equilibrium umbrella sampling kMC */
class NEUS_KMC : KMC_Enhanced_Methods {

    public:

    NEUS_KMC();
    ~NEUS_KMC();
};

/* class containing functions to propagate KMC trajectories */
class KMC_Standard_Methods {

    public:

    KMC_Standard_Methods();
    ~KMC_Standard_Methods();
    static void bkl(Walker &walker); // rejection-free algorithm of Bortz, Kalos and Lebowitz (aka n-fold way algorithm)
    static void leapfrog(Walker &walker); // leapfrog algorithm of Wales
    static void rejection_kmc(Walker &walker); // kMC algorithm where some moves are rejected
    vector<int> setup_move_probs(); // calculate the lists of move probabilities for each node

    set<int> minA, minB; // A and B endpoint nodes (A<-B)
};

#endif
