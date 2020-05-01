/*
File containing functions relating to forward flux sampling kMC (FFS-kMC).

FFS is used to simulate nonequilibrium transition path ensembles. See:
R. J. Allen, P. B. Warren, and P. R. ten Wolde, Phys. Rev. Lett. 94, 018104 (2005).
R. J. Allen, D. Frenkel and P. R. ten Wolde, J. Chem. Phys. 124, 024102 (2006).
R. J. Allen, C. Valerani and P. R. ten Wolde, J. Phys.: Condens. Matter 21, 463102 (2009).

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

FFS_KMC::FFS_KMC(const Network &ktn) {}
FFS_KMC::~FFS_KMC() {}

void FFS_KMC::run_enhanced_kmc(const Network &ktn, Traj_Method &traj_method_obj) {

}
