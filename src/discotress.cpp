/*

DISCOTRESS
DIscrete State COntinuous Time Rare Event Simulation Suite
Author: Daniel J. Sharpe (daniel.j.sharpe@gmail.com; github.com/danieljsharpe)

DISCOTRESS is a software package to simulate the dynamics on arbitrary continuous- and discrete-time Markov chains (CTMCs and DTMCs).
DISCOTRESS is designed to enable simulation of the dynamics even for Markov chains that exhibit strong metastability
(i.e. rare event dynamics), where standard simulation methods fail.

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

    void print_discotress_begin();
    void print_discotress_end();

    public:

    Discotress();
    ~Discotress();

    Network *ktn; // the network to be simulated
    Wrapper_Method *wrapper_method_obj=nullptr; // object to handle enhanced sampling
    Traj_Method *traj_method_obj=nullptr; // object to handle individual trajectories
    bool debug=false;
};

Discotress::Discotress () {

    print_discotress_begin();
    const char *inpfname = "input.kmc"; // file containing input keywords
    cout << "discotress> reading keywords..." << endl;
    Keywords my_kws = read_keywords(inpfname);
    if (my_kws.discretetime) {
        cout << "discotress> simulating a discrete-time Markov chain" << endl;
    } else {
        cout << "discotress> simulating a continuous-time Markov chain" << endl;
    }
    omp_set_num_threads(my_kws.nthreads);
    cout << "discotress> simulation will use max of " << my_kws.nthreads << " threads" << endl;
    long double dummy_randno = Wrapper_Method::rand_unif_met(my_kws.seed); // seed this generator
    if (my_kws.debug) debug=true;

    // read input files
    cout << "discotress> reading input data files..." << endl;
    const char *conns_fname="edge_conns.dat", *wts_fname="edge_weights.dat", \
               *stat_probs_fname = "stat_prob.dat";
    vector<pair<int,int>> conns = Read_files::read_two_col<int>(conns_fname);
    vector<pair<long double,long double>> weights = Read_files::read_two_col<long double>(wts_fname);
    vector<long double> stat_probs = Read_files::read_one_col<long double>(stat_probs_fname);
    vector<int> communities, bins;
    if (my_kws.commsfile!=nullptr) {
        communities = Read_files::read_one_col<int>(my_kws.commsfile);
        if (my_kws.binsfile!=nullptr) { bins = Read_files::read_one_col<int>(my_kws.binsfile);
        } else { bins = communities; } // copy community vector to bin vector
    }
    vector<int> nodesAvec, nodesBvec;
    vector<int> ntrajsvec;
    if (my_kws.wrapper_method!=2) { // simulating the A<-B TPE, read in info on A and B sets
        nodesAvec = Read_files::read_one_col<int>(my_kws.nodesafile.c_str());
        nodesBvec = Read_files::read_one_col<int>(my_kws.nodesbfile.c_str());
        if ((nodesAvec.size()!=my_kws.nA) || (nodesBvec.size()!=my_kws.nB)) {
            cout << "discotress> error: expected numbers of A and B nodes not consistent with lists of nodes in files" << endl;
            throw exception();
        }
        cout << "discotress> simulating " << my_kws.nabpaths << " transition paths. Max. no. of iterations: " << my_kws.maxit << endl;
    } else { // simulating trajectories to obtain data for coarse-graining, read in info on number of trajs for each comm
        ntrajsvec = Read_files::read_one_col<int>(my_kws.ntrajsfile);
        if (ntrajsvec.size()!=my_kws.ncomms) throw exception();
        cout << "discotress> simulating trajectories of max time: " << my_kws.trajt << "   for dimensionality reduction" << endl;
    }
    vector<double> init_probs;
    if (my_kws.initcond) init_probs = Read_files::read_one_col<double>(my_kws.initcondfile);

    // set up the Markov chain (Network) data structure
    cout << "discotress> setting up the Markovian network data object..." << endl;
    ktn = new Network(my_kws.n_nodes,my_kws.n_edges);
    if (my_kws.commsfile!=nullptr) {
        Network::setup_network(*ktn,conns,weights,stat_probs,nodesAvec,nodesBvec,my_kws.discretetime,my_kws.noloop, \
            my_kws.branchprobs,my_kws.tau,my_kws.ncomms,communities,bins);
    } else {
        Network::setup_network(*ktn,conns,weights,stat_probs,nodesAvec,nodesBvec,my_kws.discretetime,my_kws.noloop, \
            my_kws.branchprobs,my_kws.tau,my_kws.ncomms);
    }
    cout << "discotress> no. of nodes: " << ktn->n_nodes << "   in A: " << ktn->nodesA.size() << "   in B: " << ktn->nodesB.size() << endl;
    cout << "discotress> no. of edges: " << ktn->n_edges << "      no. of communities: " << ktn->ncomms << endl;
    if (my_kws.dumpwaittimes) ktn->dumpwaittimes();
    if (my_kws.initcond) ktn->set_initcond(init_probs);
    if (my_kws.statereduction && my_kws.pathlengths) { // override mean waiting times to represent mean number of steps to exit
        for (vector<Node>::iterator it_nodevec=ktn->nodes.begin();it_nodevec!=ktn->nodes.end();++it_nodevec) {
            it_nodevec->t_esc=1.L; }
    }

    // set up Traj_Method object (method to propagate trajectories associated with Walker objects)
    cout << "discotress> setting up the object to propagate individual trajectories..." << endl;
    Traj_args traj_args{my_kws.discretetime,my_kws.statereduction,my_kws.tintvl,my_kws.dumpintvls, \
                        my_kws.seed,my_kws.debug};
    if (my_kws.traj_method==1) {            // BKL algorithm
        if (my_kws.accumprobs) ktn->set_accumprobs();
        BKL *bkl_ptr = new BKL(*ktn,traj_args);
        traj_method_obj = bkl_ptr;
    } else if (my_kws.traj_method==2) {     // KPS algorithm
        KPS *kps_ptr = new KPS(*ktn,my_kws.nelim,my_kws.kpskmcsteps,my_kws.adaptivecomms,my_kws.adaptminrate,traj_args);
        if (my_kws.statereduction) {
            SR_args sr_args{my_kws.absorption,my_kws.committor,my_kws.fundamentalirred,my_kws.fundamentalred, \
                            my_kws.gth,my_kws.mfpt};
            kps_ptr->set_statereduction_procs(sr_args);
        }
        traj_method_obj = kps_ptr;
    } else if (my_kws.traj_method==3) {     // MCAMC algorithm
        MCAMC *mcamc_ptr = new MCAMC(*ktn,my_kws.kpskmcsteps,my_kws.meanrate,traj_args);
        traj_method_obj = mcamc_ptr;
    } else {
        throw exception(); // a trajectory method object must be set
    }

    // set up Wrapper_Method object (enhanced sampling method to handle set of Walker objects)
    cout << "discotress> setting up the enhanced sampling wrapper object..." << endl;
    bool indepcomms=false; // walkers correspond to independent communities or milestones
    if (my_kws.wrapper_method==2 || my_kws.wrapper_method==6) indepcomms=true;
    Wrapper_args wrapper_args{my_kws.nwalkers,ktn->nbins,my_kws.nabpaths,my_kws.tintvl,my_kws.maxit,indepcomms, \
                              my_kws.adaptivecomms,my_kws.seed,my_kws.debug};
    if (my_kws.wrapper_method==0) {        // standard simulation of A<-B paths, no enhanced sampling
        wrapper_args.nwalkers=my_kws.nthreads;
        BTOA *btoa_ptr = new BTOA(*ktn,wrapper_args);
        wrapper_method_obj = btoa_ptr;
    } else if (my_kws.wrapper_method==1) { // standard simulation of paths of fixed total time, no enhanced sampling
        if (my_kws.steadystate) wrapper_args.nwalkers=my_kws.nthreads;
        FIXEDT *fixedt_ptr = new FIXEDT(*ktn,my_kws.trajt,my_kws.steadystate,my_kws.ssrec,wrapper_args);
        wrapper_method_obj = fixedt_ptr;
    } else if (my_kws.wrapper_method==2) { // special wrapper to simulate many short nonequilibrium trajectories for dimensionality reduction
        wrapper_args.nwalkers=my_kws.nthreads;
        DIMREDN *dimredn_ptr = new DIMREDN(*ktn,ntrajsvec,my_kws.trajt,wrapper_args);
        wrapper_method_obj = dimredn_ptr;
    } else if (my_kws.wrapper_method==3) { // WE simulation
        WE *we_ptr = new WE(*ktn,my_kws.taure,wrapper_args);
        wrapper_method_obj = we_ptr;
    } else if (my_kws.wrapper_method==4) { // FFS simulation
    } else if (my_kws.wrapper_method==5) { // NEUS-kMC simulation
    } else if (my_kws.wrapper_method==6) { // milestoning simulation
    } else if (my_kws.wrapper_method==7) { // recursive enumeration algorithm for k shortest paths problem
        wrapper_args.nwalkers=0; // REA class does not store paths in walkers vector, instead has its own arrays
        REA *rea_ptr = new REA(*ktn,my_kws.discretetime,my_kws.writerea,wrapper_args);
        wrapper_method_obj = rea_ptr;
    } else {
        throw exception(); // a wrapper method object must be set
    }

    cout << "discotress> finished setting up simulation" << endl;
}

Discotress::~Discotress() {
    if (wrapper_method_obj) delete wrapper_method_obj;
    if (traj_method_obj) delete traj_method_obj;
    delete ktn;
    print_discotress_end();
}

void Discotress::print_discotress_begin() {
    // this ASCII art was produced using the FIGlet JS app (github.com/patorjk/figlet.js)
    cout << "\n\n\n\n" \
         << "          ::::::::: ::::::::::: ::::::::   ::::::::   :::::::: ::::::::::: :::::::::  :::::::::: ::::::::   :::::::: \n" \
         << "          :+:    :+:    :+:    :+:    :+: :+:    :+: :+:    :+:    :+:     :+:    :+: :+:       :+:    :+: :+:    :+:\n" \
         << "          +:+    +:+    +:+    +:+        +:+        +:+    +:+    +:+     +:+    +:+ +:+       +:+        +:+       \n" \
         << "          +#+    +:+    +#+    +#++:++#++ +#+        +#+    +:+    +#+     +#++:++#:  +#++:++#  +#++:++#++ +#++:++#++\n" \
         << "          +#+    +#+    +#+           +#+ +#+        +#+    +#+    +#+     +#+    +#+ +#+              +#+        +#+\n" \
         << "          #+#    #+#    #+#    #+#    #+# #+#    #+# #+#    #+#    #+#     #+#    #+# #+#       #+#    #+# #+#    #+#\n" \
         << "          ######### ########### ########   ########   ########     ###     ###    ### ########## ########   ######## \n" \
         << endl;
    cout << "\n\n                                                (C) Daniel J. Sharpe 2020\n\n\n\n" << endl;
}

void Discotress::print_discotress_end() {
    cout << "\n\n\n\n" \
         << "                        ########==                                               \n" \
         << "                    ####==++++++====          ________________________________   \n" \
         << "                  ====++++XXXXXX++====       |                                |_ \n" \
         << "              ==##++++++XXXX #  XX  XX##     |           Thank you for          |\n" \
         << "              ====++++XXXXXX  XX++    ##     |                                  |\n" \
         << "              ==++++++++XXXXXX++##==XX##     |               using              |\n" \
         << "            ##==####++++++++++####==XX##   __|                                  |\n" \
         << "            ##++++==##++++++####  ####  __|           D I S C O T R E S S       |\n" \
         << "            ##++++++##++++++==##  ==    \\______                              __|\n" \
         << "          ##++++++++==##++++####               |_____________________________|   \n" \
         << "          ##++++++++XX##++++==##                                                 \n" \
         << "          ##++++++XXXX##++++==##                                                 \n" \
         << "          ##++++++XX++##++++==##                                                 \n" \
         << "          ##XX    XX==##++++==##                                                 \n" \
         << "          ##XXXXXX++==##++==##                                                   \n" \
         << "          ##XX++==##==##++==##                                                   \n" \
         << "          ##=====+##==##++==##                                                   \n" \
         << "          ##====++######++####                                                   \n" \
         << "          ##==######++++____++++_________________________________________________\n" \
         << "           ##==##   ##++#   ##++#                                                \n" \
         << "           ##==##     ##      ##                                                 \n" \
         << "           ##=#                                                                  \n" \
         << "          #=#                                                                    \n" \
         << "          ##                                                                     \n" \
         << "          #                                                                      " << endl;
    cout << "\n\n" << endl;
}

int main(int argc, char** argv) {

    Discotress discotress_obj;
    if (discotress_obj.debug) run_debug_tests(*discotress_obj.ktn);
    discotress_obj.wrapper_method_obj->run_enhanced_kmc(*discotress_obj.ktn,discotress_obj.traj_method_obj);

    cout << "discotress> finished, exiting program normally" << endl;

    return 0;
}
