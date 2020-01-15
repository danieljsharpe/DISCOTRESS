/*
File containing functions relating to weighted ensemble sampling kMC (WE-kMC)
*/

#include "kmc_methods.h"
#include <omp.h>
#include <iostream>

using namespace std;

WE_KMC::WE_KMC(const Network &ktn, int maxn_abpaths, int maxit, double tau, int nbins, \
               bool adaptivebins, int seed, bool debug) {

    cout << "wekmc> running WE-kMC with parameters:\n  lag time: " << tau << " \tno. of bins: " \
         << nbins << "\n  adaptive binning (y/n): " << adaptivebins \
         << "\n  random seed: " << seed << " \tdebug printing: " << debug << endl;
    this->tau=tau; this->nbins=nbins; this->adaptivebins=adaptivebins;
    this->maxn_abpaths=maxn_abpaths; this->maxit=maxit;
    this->seed=seed; this->debug=debug;
    if (!adaptivebins) walkers.reserve(nbins); // quack use proper argument
}

WE_KMC::~WE_KMC() {}

void WE_KMC::run_enhanced_kmc(const Network &ktn) {

    cout << "wekmc> beginning WE-kMC simulation" << endl;
    n_ab=0; int n_wekmcit=0;
    double tau_r=tau; // next resampling time
    // setup walkers (set nwalkers here)
    cout << "wekmc> finished initialising walkers" << endl;
    while ((n_ab<maxn_abpaths) and (n_wekmcit<maxit)) { // algorithm terminates when max. no. of iterations of resampling procedure have been performed
        for (auto &walker: walkers) {
            if (!walker.active) continue;
//            while (walker.t<tau_r) // quack simulate
        }
        tau_r += tau;
        break;
    }
    cout << "wekmc> finished WE-kMC simulation" << endl;
}

void WE_KMC::we_resampling() {

    if (debug) cout << "wekmc> resampling trajectories" << endl;
}
