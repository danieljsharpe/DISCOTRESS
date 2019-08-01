/*
C++ code to perform enhanced kinetic Monte Carlo sampling for any generic kinetic transition network / transition rate matrix.

Compile with:
g++ -std=c++11 enhanced_kmc.cpp kmc_methods.cpp keywords.cpp -o enhanced_kmc

Daniel J. Sharpe
May 2019
*/

#include "kmc_methods.h"
#include "keywords.h"
#include <omp.h>
#include <vector>
#include <iostream>

using namespace std;

class Walker;

/* class to set up and drive Kinetic Monte Carlo simulation */
class Enhanced_KMC {

    private:

    void setup_simn(const vector<pair<int,int>>&, const vector<double>&, const vector<int>& = {});
    
    // function pointer to KMC algorithm for propagating the trajectory

    public:

    Enhanced_KMC();
    ~Enhanced_KMC();
    vector<Walker> walkers; // vector of walkers on the network

};

Enhanced_KMC::Enhanced_KMC () {

    const char *inpfname = "input.kmc";
    Keywords my_kws = read_keywords(inpfname);
    if (my_kws.enh_method==1) { // WE simulation
        walkers.reserve(my_kws.nwalkers);
    } else if (my_kws.enh_method==2) { // FFS simulation
    }
    const char *conns_fname="ts_conns.dat", *wts_fname="ts_weights.dat", \
               *comms_fname="communities.dat";
    vector<pair<int,int>> ts_conns = Read_files::read_two_col<int>(conns_fname);
    vector<double> ts_wts = Read_files::read_one_col<double>(wts_fname);
    cout << "Finished reading input files" << endl;
    setup_simn(ts_conns,ts_wts);
}

Enhanced_KMC::~Enhanced_KMC() {

}

/* Set up the data for the kinetic Monte Carlo simulation */
void Enhanced_KMC::setup_simn(const vector<pair<int,int>> &ts_conns, const vector<double> &ts_wts, \
                         const vector<int> &comms) {

}

/* class containing data for a single walker */
class Walker {

    public:

    Walker();
    ~Walker();
    int walker_id;
    int bin_curr, bin_prev; // for weighted ensemble KMC
    bool active; // walker is currently a member of the active set of trajectories being propagated

    void simulate(double);
};

Walker::Walker() {}
Walker::~Walker() {}

int main(int argc, char** argv) {

    Enhanced_KMC ekmc_obj;

    return 0;
}
