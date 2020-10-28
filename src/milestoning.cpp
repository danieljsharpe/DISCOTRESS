/*
File containing functions relating to milestoning kMC.

Milestoning is used to simulate the equilibrium transition path ensemble. See:
A. K. Faradjian and R. Elber, J. Chem. Phys. 120, 10880-10889 (2004).
J. M. Bello-Rivas and R. Elber, J. Chem. Phys. 142, 094102 (2015).

This file is a part of DISCOTRESS, a software package to simulate the dynamics on arbitrary continuous- and discrete-time Markov chains (CTMCs and DTMCs).
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

MILES::MILES(const Network &ktn, const Wrapper_args &wrapper_args) : Wrapper_Method(wrapper_args) {}
MILES::~MILES() {}

void MILES::run_enhanced_kmc(const Network &ktn, Traj_Method *traj_method_obj) {

}
