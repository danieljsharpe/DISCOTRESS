/*
functions to read in keywords

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

#include "keywords.h"
#include <omp.h>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <iostream>
#include <assert.h>

using namespace std;

Keywords read_keywords(const char *kw_file) {

    Keywords my_kws;
    string line, token;
    char delim[] = {" "};
    char *dummy;
    ifstream kw_f(kw_file);
    while (getline(kw_f,line)) {
        vector<string> vecstr;
        stringstream ss(line);
        while (getline(ss,token,*delim)) { // collect keyword + args in vector
            vecstr.emplace_back(token); }
        // PROCESS KEYWORDS
        if (vecstr.empty() || vecstr[0][0]=='!') continue; // blank or comment line
        // main keywords
        if (vecstr[0]=="NNODES") {
            my_kws.n_nodes=stoi(vecstr[1]);
        } else if (vecstr[0]=="NEDGES") {
            my_kws.n_edges=stoi(vecstr[1]);
        } else if (vecstr[0]=="WRAPPER") {
            if (vecstr[1]=="BTOA") {
                my_kws.wrapper_method=0;
            } else if (vecstr[1]=="FIXEDT") {
                my_kws.wrapper_method=1;
            } else if (vecstr[1]=="DIMREDN") {
                my_kws.wrapper_method=2;
            } else if (vecstr[1]=="WE") {
                my_kws.wrapper_method=3; 
            } else if (vecstr[1]=="FFS") {
                my_kws.wrapper_method=4;
            } else if (vecstr[1]=="NEUS") {
                my_kws.wrapper_method=5;
            } else if (vecstr[1]=="MILES") {
                my_kws.wrapper_method=6;
            } else if (vecstr[1]=="REA") {
                my_kws.wrapper_method=7;
            } else { cout << "unrecognised WRAPPER option" << endl; exit(EXIT_FAILURE); }
        } else if (vecstr[0]=="TRAJ") {
            if (vecstr[1]=="BKL") {
                my_kws.traj_method=1;
            } else if (vecstr[1]=="KPS") {
                my_kws.traj_method=2;
            } else if (vecstr[1]=="MCAMC") {
                my_kws.traj_method=3;
            } else { cout << "unrecognised TRAJ option" << endl; exit(EXIT_FAILURE); }
        } else if (vecstr[0]=="NODESAFILE") {
            my_kws.nodesafile=vecstr[1];
            my_kws.nA=stoi(vecstr[2]);
        } else if (vecstr[0]=="NODESBFILE") {
            my_kws.nodesbfile=vecstr[1];
            my_kws.nB=stoi(vecstr[2]);
        // optional keywords relating to simulation parameters and output
        } else if (vecstr[0]=="BINSFILE") {
            my_kws.binsfile = new char[vecstr[1].size()+1];
            copy(vecstr[1].begin(),vecstr[1].end(),my_kws.binsfile);
            my_kws.binsfile[vecstr[1].size()]='\0';
            my_kws.nbins=stoi(vecstr[2]);
        } else if (vecstr[0]=="COMMSFILE") {
            my_kws.commsfile = new char[vecstr[1].size()+1];
            copy(vecstr[1].begin(),vecstr[1].end(),my_kws.commsfile);
            my_kws.commsfile[vecstr[1].size()]='\0'; // trailing character
            my_kws.ncomms=stoi(vecstr[2]);
        } else if (vecstr[0]=="DUMPINTVLS") {
            my_kws.dumpintvls=true;
        } else if (vecstr[0]=="INITCONDFILE") {
            my_kws.initcondfile = new char[vecstr[1].size()+1];
            copy(vecstr[1].begin(),vecstr[1].end(),my_kws.initcondfile);
            my_kws.initcondfile[vecstr[1].size()]='\0';
            my_kws.initcond=true;
        } else if (vecstr[0]=="MAXIT") {
            my_kws.maxit=stoi(vecstr[1]);
        } else if (vecstr[0]=="NABPATHS") {
            my_kws.nabpaths=stoi(vecstr[1]);
        } else if (vecstr[0]=="TINTVL") {
            my_kws.tintvl=stod(vecstr[1]);
        // optional keywords relating to enhanced sampling methods
        } else if (vecstr[0]=="ADAPTIVECOMMS") {
            my_kws.adaptivecomms=true;
            my_kws.adaptminrate=stod(vecstr[1]);
        } else if (vecstr[0]=="COMMSTARGFILE") {
            my_kws.commstargfile = new char[vecstr[1].size()+1];
            copy(vecstr[1].begin(),vecstr[1].end(),my_kws.commstargfile);
            my_kws.commstargfile[vecstr[1].size()]='\0';
        } else if (vecstr[0]=="DIMREDUCTION") {
            my_kws.ntrajsfile = new char[vecstr[1].size()+1];
            copy(vecstr[1].begin(),vecstr[1].end(),my_kws.ntrajsfile);
            my_kws.ntrajsfile[vecstr[1].size()]='\0';
        } else if (vecstr[0]=="KPSKMCSTEPS") {
            my_kws.kpskmcsteps=stoi(vecstr[1]);
        } else if (vecstr[0]=="MEANRATE") {
            my_kws.meanrate=true;
        } else if (vecstr[0]=="NELIM") {
            my_kws.nelim=stoi(vecstr[1]);
        } else if (vecstr[0]=="NWALKERS") {
            my_kws.nwalkers=stoi(vecstr[1]);
        } else if (vecstr[0]=="TAURE") {
            my_kws.taure=stod(vecstr[1]);
        } else if (vecstr[0]=="TRAJT") {
            my_kws.trajt=stold(vecstr[1]);
        // keywords for state reduction procedures
        } else if (vecstr[0]=="ABSORPTION") {
            my_kws.absorption=true;
        } else if (vecstr[0]=="COMMITTOR") {
            my_kws.committor=true;
        } else if (vecstr[0]=="FUNDAMENTALIRRED") {
            my_kws.fundamentalirred=true;
        } else if (vecstr[0]=="FUNDAMENTALRED") {
            my_kws.fundamentalred=true;
        } else if (vecstr[0]=="GTH") {
            my_kws.gth=true;
        } else if (vecstr[0]=="MFPT") {
            my_kws.mfpt=true;
        } else if (vecstr[0]=="PATHLENGTHS") {
            my_kws.pathlengths=true;
        // other optional keywords
        } else if (vecstr[0]=="ACCUMPROBS") {
            my_kws.accumprobs=true;
        } else if (vecstr[0]=="BRANCHPROBS") {
            my_kws.branchprobs=true;
        } else if (vecstr[0]=="DEBUG") {
            my_kws.debug=true;
        } else if (vecstr[0]=="DISCRETETIME") {
            my_kws.discretetime=true;
        } else if (vecstr[0]=="DUMPWAITTIMES") {
            my_kws.dumpwaittimes=true;
        } else if (vecstr[0]=="NTHREADS") {
            my_kws.nthreads=stoi(vecstr[1]);
            assert((my_kws.nthreads>0 && my_kws.nthreads<=omp_get_max_threads()));
        } else if (vecstr[0]=="SEED") {
            my_kws.seed=stoi(vecstr[1]);
        } else if (vecstr[0]=="TAU") {
            my_kws.tau=stold(vecstr[1]);
        } else {
            cout << "keywords> error: unrecognised keyword: " << vecstr[0] << endl;
            exit(EXIT_FAILURE);
        }
    }
    kw_f.close();
    cout << "keywords> finished reading keywords" << endl;

    // check necessary keywords and compatability
    if (my_kws.n_nodes<=0 || my_kws.n_edges<=0 || ((my_kws.nA<=0 || my_kws.nB<=0) && my_kws.wrapper_method!=2)) {
        cout << "keywords> error: network parameters not set correctly" << endl; exit(EXIT_FAILURE); }
    if ((my_kws.nabpaths<=0 && my_kws.wrapper_method!=2) || my_kws.maxit<=0) {
        cout << "keywords> error: termination condition not specified correctly" << endl; exit(EXIT_FAILURE); }
    if (my_kws.commsfile!=nullptr && my_kws.ncomms<=1) {
        cout << "keywords> error: there must be at least two communities in the specified partitioning" << endl; exit(EXIT_FAILURE); }
    if (my_kws.dumpintvls && my_kws.tintvl<=0.) {
        cout << "keywords> error: invalid time interval for dumping trajectory data" << endl; exit(EXIT_FAILURE); }
    if (my_kws.traj_method<=0 || my_kws.wrapper_method<0) {
        cout << "keywords> error: must specify both a wrapper method and a trajectory method" << endl; exit(EXIT_FAILURE); }
    if ((my_kws.discretetime || !my_kws.branchprobs) && my_kws.tau<=0.) {
        cout << "keywords> error: if reading in transition probs for DTMC or otherwise not using branching probs, must specify tau as lag time" << endl;
        exit(EXIT_FAILURE); }
    // set the purpose of the computation to be a state reduction procedure and not a dynamical simulation, if appropriate
    if (my_kws.committor || my_kws.absorption || my_kws.fundamentalred || my_kws.fundamentalirred || my_kws.mfpt || my_kws.gth) {
        assert(my_kws.nabpaths==1);
        if (my_kws.wrapper_method!=0 || my_kws.traj_method!=2) {
            cout << "keywords> error: to perform a state reduction computation, must set WRAPPER BTOA and TRAJ KPS" << endl; exit(EXIT_FAILURE); }
        if (my_kws.ncomms!=2) {
            cout << "keywords> error: a state reduction computation uses only two communities (namely, not A and A)" << endl; exit(EXIT_FAILURE); }
        if ((my_kws.gth || my_kws.fundamentalirred) && my_kws.nA!=1) {
            cout << "keywords> error: the GTH and FUND algorithms can be ran only when there is a single node in A" << endl; exit(EXIT_FAILURE); }
        if (my_kws.fundamentalred && (my_kws.committor || my_kws.absorption || my_kws.fundamentalirred || my_kws.mfpt || my_kws.gth)) {
            cout << "keywords> error: computation of the fundamental matrix for a reducible Markov chain is standalone" << endl; exit(EXIT_FAILURE); }
        if (my_kws.n_nodes-my_kws.nA>my_kws.nelim) {
            cout << "keywords> error: for state reduction must set NELIM to ensure that all nodes not in A are eliminated" << endl; exit(EXIT_FAILURE); }
        my_kws.statereduction=true;
        my_kws.nthreads=1; // use only a single thread for a state reduction computation
    }
    // check specification of wrapper method is valid
    if (my_kws.wrapper_method==0) { // standard simulation of paths initialised in state B and terminating when state A is hit
        // ...
    } else if (my_kws.wrapper_method==1) { // standard simulation of paths with fixed total time
        assert(my_kws.trajt>my_kws.ssrec);
        if (my_kws.trajt<=0. || my_kws.ssrec<0.) {
            cout << "keywords> error: simulation of fixed-time paths not specified correctly" << endl; exit(EXIT_FAILURE); }
    } else if (my_kws.wrapper_method==2) { // special wrapper method to propagate trajectories required for dimensionality reduction
        if (my_kws.ntrajsfile==nullptr || my_kws.trajt<=0. || my_kws.commsfile==nullptr || my_kws.meanrate || my_kws.initcondfile || \
            my_kws.traj_method==1 || my_kws.nA!=0 || my_kws.nB!=0 || !my_kws.dumpintvls) {
            cout << "keywords> error: dimensionality reduction simulation not specified correctly" << endl; exit(EXIT_FAILURE); }
    } else if (my_kws.wrapper_method==3) { // WE simulation
        if (my_kws.taure<=0. || (my_kws.commsfile!=nullptr && !my_kws.adaptivecomms) || my_kws.nwalkers<1 || \
            (my_kws.commstargfile!=nullptr && !my_kws.adaptivecomms) || my_kws.traj_method!=1) {
            cout << "keywords> error: WE simulation not specified correctly" << endl; exit(EXIT_FAILURE); }
    } else if (my_kws.wrapper_method==4) { // FFS simulation
        if (my_kws.commsfile==nullptr) {
            cout << "keywords> error: FFS simulation not specified correctly" << endl; exit(EXIT_FAILURE); }
    } else if (my_kws.wrapper_method==5 || my_kws.nwalkers<1) { // NEUS simulation
        if (my_kws.commsfile==nullptr || my_kws.traj_method!=1) {
            cout << "keywords> error: NEUS simulation not specified correctly" << endl; exit(EXIT_FAILURE); }
    } else if (my_kws.wrapper_method==6 || my_kws.nwalkers<1) { // milestoning simulation
        if (my_kws.commsfile==nullptr) {
            cout << "keywords> error: milestoning simulation not specified correctly" << endl; exit(EXIT_FAILURE); }
    } else if (my_kws.wrapper_method==7) { // recursive enumeration algorithm for k shortest paths
        if (my_kws.nA!=1 || my_kws.nB!=1 || my_kws.nabpaths<1) {
            cout << "keywords> error: REA k shortest paths computation not specified correctly" << endl; exit(EXIT_FAILURE); }
    }
    // check specification of trajectory method is valid
    if (my_kws.traj_method==1) { // BKL algorithm
        // ...
    } else if (my_kws.traj_method==2) { // kPS algorithm
        if ((my_kws.commsfile==nullptr && !my_kws.adaptivecomms) || my_kws.nelim<=0) {
            cout << "keywords> error: kPS algorithm not specified correctly" << endl; exit(EXIT_FAILURE); }
    } else if (my_kws.traj_method==3) { // MCAMC algorithm
        if (my_kws.branchprobs) {
            cout << "keywords> error: MCAMC algorithm not specified correctly" << endl; exit(EXIT_FAILURE); }
    }

    return my_kws;
}
