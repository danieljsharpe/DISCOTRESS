/*

DISCOTRESS
DIscrete State COntinuous Time Rare Event Simulation Suite
Author: Daniel J. Sharpe (daniel.j.sharpe@gmail.com; github.com/danieljsharpe)

DISCOTRESS is a software package to simulate the dynamics on arbitrary continuous time Markov chains (CTMCs).
DISCOTRESS is designed to enable simulation of the dynamics even for CTMCs that exhibit strong metastability
(ie rare event dynamics), where standard simulation methods fail.

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
#include "keywords.h"
#include "debug_tests.h"
#include <vector>
#include <iostream>

using namespace std;

/* main class to set up and drive a standard or enhanced kinetic Monte Carlo simulation */
class Discotress {

    private:

    public:

    Discotress();
    ~Discotress();

    Network *ktn; // the network to be simulated
    Wrapper_Method *wrapper_method_obj=nullptr; // object to handle enhanced sampling
    Traj_Method *traj_method_obj=nullptr; // object to handle individual trajectories
    bool debug=false;
};

Discotress::Discotress () {

    const char *inpfname = "input.kmc"; // file containing input keywords
    cout << "discotress> reading keywords..." << endl;
    Keywords my_kws = read_keywords(inpfname);
    cout << "discotress> reading input data files..." << endl;
    const char *conns_fname="ts_conns.dat", *wts_fname="ts_weights.dat", \
               *stat_probs_fname = "stat_prob.dat";
    vector<pair<int,int>> ts_conns = Read_files::read_two_col<int>(conns_fname);
    vector<long double> ts_wts = Read_files::read_one_col<long double>(wts_fname);
    vector<double> stat_probs = Read_files::read_one_col<double>(stat_probs_fname);
    vector<int> communities, bins;
    if (my_kws.commsfile!=nullptr) {
        communities = Read_files::read_one_col<int>(my_kws.commsfile);
        if (my_kws.binfile!=nullptr) { bins = Read_files::read_one_col<int>(my_kws.binfile);
        } else { bins = communities; } // copy community vector to bin vector
    }
    vector<int> nodesA, nodesB;
    vector<int> ntrajsvec;
    if (my_kws.ntrajsfile==nullptr) { // simulating the A<-B TPE, read in info on A and B sets
        nodesA = Read_files::read_one_col<int>(my_kws.nodesafile.c_str());
        nodesB = Read_files::read_one_col<int>(my_kws.nodesbfile.c_str());
        if (!((nodesA.size()==my_kws.nA) || (nodesB.size()==my_kws.nB))) throw exception();
    } else { // simulating trajectories to obtain data for coarse-graining, read in info on number of trajs for each comm
        ntrajsvec = Read_files::read_one_col<int>(my_kws.ntrajsfile);
        if (ntrajsvec.size()!=my_kws.ncomms) throw exception();
    }
    vector<double> init_probs;
    if (my_kws.initcond) init_probs = Read_files::read_one_col<double>(my_kws.initcondfile);
    cout << "discotress> setting up the transition network data structure..." << endl;
    ktn = new Network(my_kws.n_nodes,my_kws.n_edges);
    if (my_kws.commsfile!=nullptr) {
        Network::setup_network(*ktn,ts_conns,ts_wts,stat_probs,nodesA,nodesB,my_kws.transnprobs, \
            my_kws.tau,my_kws.ncomms,communities,bins);
    } else {
        Network::setup_network(*ktn,ts_conns,ts_wts,stat_probs,nodesA,nodesB,my_kws.transnprobs,my_kws.tau,my_kws.ncomms);
    }
    if (my_kws.dumpwaittimes) ktn->dumpwaittimes();
    if (my_kws.initcond) ktn->set_initcond(init_probs);
    cout << "discotress> setting up the simulator object..." << endl;
    // set up wrapper enhanced sampling class
    if (my_kws.wrapper_method==0) {        // special wrapper to simulate many short nonequilibrium trajectories for dimensionality reduction
        DIMREDN *dimredn_ptr = new DIMREDN(*ktn,ntrajsvec,my_kws.dt,my_kws.seed);
        wrapper_method_obj = dimredn_ptr;
    } else if (my_kws.wrapper_method==1) {        // standard kMC simulation, no enhanced sampling
        STD_KMC *std_kmc_ptr = new STD_KMC(*ktn,my_kws.nabpaths,my_kws.maxit,my_kws.tintvl,my_kws.adaptivecomms,my_kws.seed);
        wrapper_method_obj = std_kmc_ptr;
    } else if (my_kws.wrapper_method==2) { // WE simulation
        WE_KMC *we_kmc_ptr = new WE_KMC(*ktn,my_kws.nabpaths,my_kws.maxit,my_kws.tau,my_kws.tintvl,my_kws.adaptivecomms, \
                    my_kws.seed,my_kws.debug);
        wrapper_method_obj = we_kmc_ptr;
    } else if (my_kws.wrapper_method==3) { // FFS simulation
    } else if (my_kws.wrapper_method==4) { // NEUS-kMC simulation
    } else if (my_kws.wrapper_method==5) { // milestoning simulation
    } else {
        throw exception(); // a wrapper method object must be set
    }
    // set up method to propagate individual trajectories
    if (my_kws.traj_method==1) {            // BKL algorithm
        ktn->get_cum_branchprobs();
        wrapper_method_obj->set_standard_kmc(&BKL::bkl);
    } else if (my_kws.traj_method==2) {     // kPS algorithm
        if (my_kws.branchprobs) { ktn->get_tmtx_branch();
        } else if (!my_kws.transnprobs) { ktn->get_tmtx_lin(my_kws.tau); }
        KPS *kps_ptr = new KPS(*ktn,my_kws.nelim,my_kws.tau,my_kws.kpskmcsteps, \
                    my_kws.adaptivecomms,my_kws.adaptminrate,my_kws.pfold,my_kws.tintvl,my_kws.seed,my_kws.debug);
        traj_method_obj = kps_ptr;
    } else if (my_kws.traj_method==3) {     // MCAMC algorithm
        MCAMC *mcamc_ptr = new MCAMC(*ktn,my_kws.kpskmcsteps,my_kws.meanrate);
        traj_method_obj = mcamc_ptr;
    }
    if (my_kws.debug) debug=true;
    cout << "discotress> finished setting up simulation" << endl;
}

Discotress::~Discotress() {
    if (wrapper_method_obj) delete wrapper_method_obj;
    if (traj_method_obj) delete traj_method_obj;
    delete ktn;
}


int main(int argc, char** argv) {

    Discotress discotress_obj;
    if (discotress_obj.debug) run_debug_tests(*discotress_obj.ktn);
    discotress_obj.wrapper_method_obj->run_enhanced_kmc(*discotress_obj.ktn,*discotress_obj.traj_method_obj);

/*
    Node newnode = discotress_obj.ktn->nodes[0];
    cout << newnode.node_id << "   " << newnode.udeg << endl;
*/

    cout << "discotress> finished, exiting program normally" << endl;
    return 0;
}
