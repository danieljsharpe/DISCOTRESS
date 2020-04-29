/*
File containing functions relating to weighted ensemble sampling kMC (WE-kMC).

WE sampling can be used to sample nonequilibrium and equilibrium transition path ensembles. See:
G. A. Huber and S. Kim, Biophys. J. 70, 97-110 (1996).
R. M. Donovan, A. J. Sedgewick, J. R. Faeder, and D. M. Zuckerman, J. Chem. Phys. 139, 115105 (2013).
D. M. Zuckerman and L. T. Chong, Annu. Rev. Biophys. 46, 43-57 (2017).

This file is a part of DISCOTRESS, a software package to simulate the dynamics on arbitrary continuous time Markov chains (CTMCs).
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

#include "kmc_methods.h"
#include <omp.h>
#include <iostream>

using namespace std;

WE_KMC::WE_KMC(const Network &ktn, int maxn_abpaths, int maxit, long double tau, double tintvl, \
               bool adaptivecomms, int seed, bool debug) {

    cout << "wekmc> running WE-kMC with parameters:\n  lag time: " << tau << " \tno. of communities: " \
         << ktn.ncomms << "\n  adaptive definition of communities (y/n): " << adaptivecomms \
         << "\n  random seed: " << seed << " \tdebug printing: " << debug << endl;
    this->tau=static_cast<double>(tau); this->adaptivecomms=adaptivecomms;
    this->maxn_abpaths=maxn_abpaths; this->maxit=maxit; this->tintvl=tintvl;
    this->seed=seed; this->debug=debug;
    if (!adaptivecomms) walkers.reserve(ktn.ncomms); // quack use proper argument and set accumprobs
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
