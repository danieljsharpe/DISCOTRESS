/*
C++ code providing a suite for kinetic Monte Carlo simulations, including various advanced sampling strategies,
applicable to any generic kinetic transition network

Compile with:
g++ -std=c++17 kmc_suite.cpp kmc_methods.cpp we_kmc.cpp kps.cpp ffs_kmc.cpp as_kmc.cpp neus_kmc.cpp keywords.cpp ktn.cpp -o kmc_suite

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

    public:

    KMC_Suite();
    ~KMC_Suite();

    Network *ktn; // the network to be simulated
    KMC_Enhanced_Methods *enh_method=nullptr; // object to handle enhanced sampling
    void *kmc_func; // function pointer to KMC algorithm for propagating the trajectory
    bool debug=false;
};

KMC_Suite::KMC_Suite () {

    const char *inpfname = "input.kmc"; // file containing input keywords
    cout << "kmc_suite> reading keywords..." << endl;
    Keywords my_kws = read_keywords(inpfname);
    cout << "kmc_suite> reading input data files..." << endl;
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

    cout << "kmc_suite> setting up the transition network data structure..." << endl;
    ktn = new Network(my_kws.n_nodes,my_kws.n_edges);
    if (strlen(my_kws.binfile)>0) {
        Network::setup_network(*ktn,ts_conns,ts_wts,stat_probs,nodesA,nodesB,communities);
    } else {
        Network::setup_network(*ktn,ts_conns,ts_wts,stat_probs,nodesA,nodesB);
    }

    cout << "kmc_suite> setting up simulation..." << endl;
    cout << "kmc_suite> max. no. of A-B paths: " << my_kws.nabpaths << " \tmax. iterations: " << my_kws.maxit << "\n" << endl;
    // set up enhanced sampling class, if any chosen
    if (my_kws.enh_method==1) {        // WE simulation
        WE_KMC *we_kmc_ptr = new WE_KMC(*ktn,my_kws.tau,my_kws.nbins,my_kws.adaptivebins);
        enh_method = we_kmc_ptr;
    } else if (my_kws.enh_method==2) { // kPS simulation
        ktn->get_tmtx_lin(my_kws.tau);
        KPS *kps_ptr = new KPS(*ktn,my_kws.nabpaths,my_kws.maxit,my_kws.nelim,my_kws.tau,my_kws.nbins,my_kws.kpskmcsteps, \
                    my_kws.adaptivebins,my_kws.initcond,my_kws.seed,my_kws.debug);
        enh_method = kps_ptr;
    } else if (my_kws.enh_method==3) { // FFS simulation
    } else if (my_kws.enh_method==4) { // AS-kMC simulation
    } else if (my_kws.enh_method==5) { // NEUS-kMC simulation
    } else {
        enh_method=nullptr;
    }
    // set up method to propagate kMC trajectories
    if (my_kws.enh_method==2) {        // kPS algorithm (does not use explicit kMC simulation)
        kmc_func=nullptr;
    } else if (my_kws.kmc_method==1) { // BKL algorithm
    } else if (my_kws.kmc_method==2) { // rejection algorithm
    } else if (my_kws.kmc_method==3) { // leapfrog algorithm
    }
    if (my_kws.debug) debug=true;
}

KMC_Suite::~KMC_Suite() {
    delete ktn;
    if (enh_method) delete enh_method;
}


int main(int argc, char** argv) {

    KMC_Suite kmc_suite_obj;
    if (kmc_suite_obj.debug) run_debug_tests(*kmc_suite_obj.ktn);
    kmc_suite_obj.enh_method->run_enhanced_kmc(*kmc_suite_obj.ktn);
/*
    Node newnode = kmc_suite_obj.ktn->nodes[0];
    cout << newnode.node_id << "   " << newnode.udeg << endl;
*/

    cout << "kmc_suite> finished program, exiting" << endl;
    return 0;
}
