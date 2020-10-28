/*
File containing functions relating to non-equilibrium umbrella sampling kMC (NEUS-kMC)

NEUS is used to simulate the equilibrium transition path ensemble. See:
A. Dickson and A. Warmflash and A. R. Dinner, J. Chem. Phys. 130, 074104 (2009).
A. Dickson and A. Warmflash and A. R. Dinner, J. Chem. Phys. 131, 154104 (2009).
A. Dickson and A. R. Dinner, Annu. Rev. Phys. Chem. 61, 441-459 (2010).
E. Vanden-Eijnden and M. Venturoli, J. Chem. Phys. 131, 044120 (2009).

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

using namespace std;

NEUS::NEUS(const Network &ktn, const Wrapper_args &wrapper_args) : Wrapper_Method(wrapper_args) {}
NEUS::~NEUS() {}

void NEUS::run_enhanced_kmc(const Network &ktn, Traj_Method *traj_method_obj) {

}
