/*
File containing functions relating to transition path sampling kMC

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

TPS_KMC::TPS_KMC(const Network &ktn) {}
TPS_KMC::~TPS_KMC() {}

void TPS_KMC::run_enhanced_kmc(const Network &ktn) {

}
