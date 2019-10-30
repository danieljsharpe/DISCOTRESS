/*
File containing kinetic Monte Carlo simulation algorithms to propagate individual trajectories
*/

#include "kmc_methods.h"

using namespace std;

Walker::Walker() {}
Walker::~Walker() {}

KMC_Enhanced_Methods::KMC_Enhanced_Methods() {}
KMC_Enhanced_Methods::~KMC_Enhanced_Methods() {}

void KMC_Enhanced_Methods::find_bin_onthefly() {

}

KMC_Standard_Methods::KMC_Standard_Methods() {}
KMC_Standard_Methods::~KMC_Standard_Methods() {}

void KMC_Standard_Methods::bkl(Walker &walker) {

}

void KMC_Standard_Methods::leapfrog(Walker &walker) {

}

void KMC_Standard_Methods::rejection_kmc(Walker &walker) {

}

vector<int> KMC_Standard_Methods::setup_move_probs(Network &ktn) {

}
