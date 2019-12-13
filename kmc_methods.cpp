/*
File containing kinetic Monte Carlo simulation algorithms to propagate individual trajectories
*/

#include "kmc_methods.h"
#include <random>

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

/* draw a uniform random number between 0 and 1, used in Metropolis conditions etc. */
double KMC_Standard_Methods::rand_unif_met(int seed) {
    static default_random_engine generator(seed);
    static uniform_real_distribution<double> unif_real_distrib(0.,1.);
    return unif_real_distrib(generator);
}

vector<int> KMC_Standard_Methods::setup_move_probs(const Network &ktn) {

}
