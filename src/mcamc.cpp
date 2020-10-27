/*
File containing functions relating to the Monte Carlo with absorbing Markov chains (MCAMC) algorithm. See:
M. A. Novotny, Phys. Rev. Lett. 74, 1-5 (1995).
M. A. Novotny, Comput. Phys. Commun. 147, 659-664 (2002).
B. Puchala, M. L. Falk and K. Garikipati, J. Chem. Phys. 132, 134104 (2010).

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

MCAMC::MCAMC(const Network &ktn, int kpskmcsteps, bool meanrate, \
             const Traj_args &traj_args) : Traj_Method(traj_args) {

    cout << "kps> MCAMC parameters:\n  FPTA (0) or mean rate method (1)?: " << meanrate \
         << "\n  no. of kMC steps after MCAMC iteration: " << kpskmcsteps << endl;
    this->kpskmcsteps=kpskmcsteps; this->meanrate=meanrate;
}

MCAMC::~MCAMC() {}

MCAMC::MCAMC(const MCAMC &mcamc_obj) : Traj_Method(mcamc_obj) {
    this->kpskmcsteps=mcamc_obj.kpskmcsteps; this->meanrate=mcamc_obj.meanrate;
}

void MCAMC::kmc_iteration(const Network &ktn, Walker &walker) {

    cout << "mcamc> running a single iteration of MCAMC" << endl;
}

/* perform specified number of BKL iterations after a basin escape */
void MCAMC::do_bkl_steps(const Network &ktn, Walker &walker, long double maxtime) {

    cout << "mcamc> doing BKL steps after a basin escape" << endl;
}

void MCAMC::reset_nodeptrs() {
    epsilon=nullptr; alpha=nullptr;
}
