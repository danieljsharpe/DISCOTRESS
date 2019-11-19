/*
C++ code providing a suite for kinetic Monte Carlo simulations, including various advanced sampling strategies,
applicable to any generic kinetic transition network

Compile with:
g++ -std=c++11 kmc_suite.cpp kmc_methods.cpp we_kmc.cpp kps.cpp ffs_kmc.cpp as_kmc.cpp neus_kmc.cpp keywords.cpp ktn.cpp -o kmc_suite

Daniel J. Sharpe
May 2019
*/

#include "kmc_methods.h"
#include "keywords.h"
#include "debug_tests.h"
#include <omp.h>
#include <vector>
#include <iostream>
#include <cstring>

using namespace std;

/* main class to set up and drive a standard or enhanced kinetic Monte Carlo simulation */
class KMC_Suite {

    private:

    // function pointer to KMC algorithm for propagating the trajectory

    public:

    KMC_Suite();
    ~KMC_Suite();

    Network *ktn; // the network to be simulated
    // object to handle enhanced sampling
};

KMC_Suite::KMC_Suite () {

    const char *inpfname = "input.kmc"; // file containing input keywords
    cout << ">>>>> reading keywords..." << endl;
    Keywords my_kws = read_keywords(inpfname);
    cout << ">>>>> setting up transition network..." << endl;
    const char *conns_fname="ts_conns.dat", *wts_fname="ts_weights.dat", \
               *stat_probs_fname = "stat_prob.dat";
    vector<pair<int,int>> ts_conns = Read_files::read_two_col<int>(conns_fname);
    vector<double> ts_wts = Read_files::read_one_col<double>(wts_fname);
    vector<double> stat_probs = Read_files::read_one_col<double>(stat_probs_fname);
    vector<int> communities;
    if (strlen(my_kws.binfile)>0) communities = Read_files::read_one_col<int>(my_kws.binfile);
    vector<int> nodesA = Read_files::read_one_col<int>(my_kws.minafile.c_str());
    vector<int> nodesB = Read_files::read_one_col<int>(my_kws.minbfile.c_str());
    if (!((nodesA.size()==my_kws.nA) || (nodesB.size()==my_kws.nB))) throw exception();
    cout << ">>>>> setting up the transition network data structure" << endl;
    Network ktn_obj(my_kws.n_nodes,my_kws.n_edges);
    if (strlen(my_kws.binfile)>0) {
        Network::setup_network(ktn_obj,ts_conns,ts_wts,stat_probs,nodesA,nodesB,communities);
    } else {
        Network::setup_network(ktn_obj,ts_conns,ts_wts,stat_probs,nodesA,nodesB);
    }
    ktn = &ktn_obj;
    cout << ">>>>> setting up simulation..." << endl;
    // set up enhanced sampling class, if any chosen
    if (my_kws.enh_method==1) {        // WE simulation
        // do stuff to set up WE-kMC simulation
        // we_kmc_obj = WE_KMC()
        // we_kmc_obj.walkers.reserve(my_kws.nwalkers);
    } else if (my_kws.enh_method==2) { // kPS simulation
        KPS kps_obj(*ktn,my_kws.nelim,my_kws.tau);
    } else if (my_kws.enh_method==3) { // FFS simulation
    } else if (my_kws.enh_method==4) { // AS-kMC simulation
    } else if (my_kws.enh_method==5) { // NEUS-kMC simulation
    }
    // set up method to propagate kMC trajectories
    if (my_kws.enh_method==2) {        // kPS algorithm (does not use explicit kMC simulation)
    } else if (my_kws.kmc_method==1) { // BKL algorithm
    } else if (my_kws.kmc_method==2) { // rejection algorithm
    } else if (my_kws.kmc_method==3) { // leapfrog algorithm
    }

    if (my_kws.debug) run_debug_tests(ktn_obj);
}

KMC_Suite::~KMC_Suite() {

}


int main(int argc, char** argv) {

    KMC_Suite kmc_suite_obj;

    return 0;
}
