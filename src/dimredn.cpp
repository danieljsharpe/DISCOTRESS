/*
File containing functions to handle simulation of many short nonequilibrium trajectories, used to obtain data required for coarse-graining
a transition network

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

DIMREDN::DIMREDN(const Network &ktn, vector<int> ntrajsvec, double dt, int seed) {

    cout << "dimredn> constructing DIMREDN class" << endl;
    this->ntrajsvec=ntrajsvec; this->dt=dt;
    this->seed=seed;
}

DIMREDN::~DIMREDN() {}

void DIMREDN::run_enhanced_kmc(const Network &ktn, Traj_Method &traj_method_obj) {

    cout << "dimredn> beginning simulation to obtain trajectory data for dimensionality reduction" << endl;
    for (int i=0;i<ntrajsvec.size();i++) {
        cout << "i: " << ntrajsvec[i] << endl; }
}
