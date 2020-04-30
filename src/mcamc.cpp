/*
File containing functions relating to the Monte Carlo with absorbing Markov chains (MCAMC) algorithm. See:
M. A. Novotny, Phys. Rev. Lett. 74, 1-5 (1995).
M. A. Novotny, Comput. Phys. Commun. 147, 659-664 (2002).
B. Puchala, M. L. Falk and K. Garikipati, J. Chem. Phys. 132, 134104 (2010).

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

using namespace std;

MCAMC::MCAMC(const Network &ktn, int maxn_abpaths, int maxit, double tintvl, bool meanrate) {

    if (maxn_abpaths>0) {
        cout << "mcamc> simulating the A<-B transition path ensemble with MCAMC (FPTA) algorithm:\n" \
             << "  max. no. of A<-B paths: " << maxn_abpaths << " \tmax. iterations: " << maxit \
             << "\n  time interval for dumping trajectory data: " << tintvl << endl;
    } else {
        cout << "mcamc> simulating many short trajectories with MCAMC (FPTA) algorithm" << endl;
    }
    cout << "kps> MCAMC parameters:\n  FPTA (0) or mean rate method (1)?: " << meanrate << endl;
    this->meanrate=meanrate;
}

MCAMC::~MCAMC() {}

void MCAMC::run_enhanced_kmc(const Network &ktn) {

    cout << "mcamc> beginning MCAMC simulation to sample the A<-B TPE" << endl;
}

void MCAMC::run_dimreduction(const Network &ktn, vector<int> ntrajsvec) {

    cout << "mcamc> beginning MCAMC simulation to obtain trajectory data for dimensionality reduction" << endl;
}
