/*
functions to read in keywords

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

#include "keywords.h"
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
        if (vecstr[0]=="TRAJ") {
            if (vecstr[1]=="BKL") {
                my_kws.traj_method=1;
            } else if (vecstr[1]=="KPS") {
                my_kws.traj_method=2;
            } else if (vecstr[2]=="MCAMC") {
                my_kws.traj_method=3;
            } else { exit(EXIT_FAILURE); }
        } else if (vecstr[0]=="WRAPPER") {
            if (vecstr[1]=="DIMREDN") {
                my_kws.wrapper_method=0;
            } else if (vecstr[1]=="NONE") {
                my_kws.wrapper_method=1;
            } else if (vecstr[1]=="WE") {
                my_kws.wrapper_method=2;
            } else if (vecstr[2]=="FFS") {
                my_kws.wrapper_method=3;
            } else if (vecstr[3]=="NEUS") {
                my_kws.wrapper_method=4;
            } else if (vecstr[4]=="MILES") {
                my_kws.wrapper_method=5;
            } else { exit(EXIT_FAILURE); }
        } else if (vecstr[0]=="NNODES") {
            my_kws.n_nodes=stoi(vecstr[1]);
        } else if (vecstr[0]=="NEDGES") {
            my_kws.n_edges=stoi(vecstr[1]);
        } else if (vecstr[0]=="NABPATHS") {
            my_kws.nabpaths=stoi(vecstr[1]);
        } else if (vecstr[0]=="MAXIT") {
            my_kws.maxit=stoi(vecstr[1]);
        } else if (vecstr[0]=="NODESAFILE") {
            my_kws.nodesafile=vecstr[1];
            my_kws.nA=stoi(vecstr[2]);
        } else if (vecstr[0]=="NODESBFILE") {
            my_kws.nodesbfile=vecstr[1];
            my_kws.nB=stoi(vecstr[2]);
        } else if (vecstr[0]=="DIMREDUCTION") {
            my_kws.ntrajsfile = new char[vecstr[1].size()+1];
            copy(vecstr[1].begin(),vecstr[1].end(),my_kws.ntrajsfile);
            my_kws.ntrajsfile[vecstr[1].size()]='\0';
            my_kws.dt=stod(vecstr[2]);
        } else if (vecstr[0]=="TAU") {
            my_kws.tau=stold(vecstr[1]);
        } else if (vecstr[0]=="TINTVL") {
            my_kws.tintvl=stod(vecstr[1]);
        } else if (vecstr[0]=="INITCONDFILE") {
            my_kws.initcondfile = new char[vecstr[1].size()+1];
            copy(vecstr[1].begin(),vecstr[1].end(),my_kws.initcondfile);
            my_kws.initcondfile[vecstr[1].size()]='\0';
            my_kws.initcond=true;
        } else if (vecstr[0]=="COMMSFILE") {
            my_kws.commsfile = new char[vecstr[1].size()+1];
            copy(vecstr[1].begin(),vecstr[1].end(),my_kws.commsfile);
            my_kws.commsfile[vecstr[1].size()]='\0'; // trailing character
            my_kws.ncomms=stoi(vecstr[2]);
        } else if (vecstr[0]=="COMMSTARGFILE") {
            my_kws.commstargfile = new char[vecstr[1].size()+1];
            copy(vecstr[1].begin(),vecstr[1].end(),my_kws.commstargfile);
            my_kws.commstargfile[vecstr[1].size()]='\0';
        } else if (vecstr[0]=="BINFILE") {
            my_kws.binfile = new char[vecstr[1].size()+1];
            copy(vecstr[1].begin(),vecstr[1].end(),my_kws.binfile);
            my_kws.binfile[vecstr[1].size()]='\0';
            my_kws.nbins=stoi(vecstr[2]);
        } else if (vecstr[0]=="ADAPTIVECOMMS") {
            my_kws.adaptivecomms=true;
            my_kws.adaptminrate=stod(vecstr[1]);
        } else if (vecstr[0]=="KPSKMCSTEPS") {
            my_kws.kpskmcsteps=stoi(vecstr[1]);
        } else if (vecstr[0]=="MEANRATE") {
            my_kws.meanrate=true;
        } else if (vecstr[0]=="NELIM") {
            my_kws.nelim=stoi(vecstr[1]);
        } else if (vecstr[0]=="PFOLD") {
            my_kws.pfold=true;
        } else if (vecstr[0]=="TRANSNPROBS") {
            my_kws.transnprobs=true;
        } else if (vecstr[0]=="BRANCHPROBS") {
            my_kws.branchprobs=true;
        } else if (vecstr[0]=="NTHREADS") {
            my_kws.nthreads=stoi(vecstr[1]);
            assert(my_kws.nthreads>0);
        } else if (vecstr[0]=="DEBUG") {
            my_kws.debug=true;
        } else if (vecstr[0]=="SEED") {
            my_kws.seed=stoi(vecstr[1]);
        } else if (vecstr[0]=="DUMPWAITTIMES") {
            my_kws.dumpwaittimes=true;
        } else {
            cout << "keywords> error: unrecognised keyword: " << vecstr[0] << endl;
            exit(EXIT_FAILURE);
        }
    }
    kw_f.close();
    cout << "keywords> finished reading keywords" << endl;

    // check necessary keywords and compatability
    if (my_kws.n_nodes<=0 || my_kws.n_edges<=0 || ((my_kws.nA<=0 || my_kws.nB<=0) && my_kws.ntrajsfile==nullptr)) {
        cout << "keywords> error: transition network parameters not set correctly" << endl; exit(EXIT_FAILURE); }
    if ((my_kws.nabpaths<=0 && my_kws.ntrajsfile==nullptr) || (my_kws.dt<=0. && my_kws.ntrajsfile!=nullptr) || my_kws.maxit<=0) {
        cout << "keywords> error: termination condition not specified correctly" << endl; exit(EXIT_FAILURE); }
    if (my_kws.commsfile!=nullptr && my_kws.ncomms<=1) {
        cout << "keywords> error: there must be at least two communities in the specified partitioning" << endl; exit(EXIT_FAILURE); }
    if (my_kws.traj_method<=0 || my_kws.wrapper_method<=0) {
        cout << "keywords> error: must specify both a wrapper method and a trajectory method" << endl; exit(EXIT_FAILURE); }
    if (my_kws.transnprobs && my_kws.traj_method!=2) {
        cout << "keywords> error: edge weights must be read in as transition rates if not using kPS" << endl; exit(EXIT_FAILURE); }
    // check specification of wrapper method is valid
    if (my_kws.wrapper_method==0) { // special wrapper method to propagate trajectories required for dimensionality reduction
        if (my_kws.ntrajsfile==nullptr || my_kws.commsfile==nullptr || my_kws.meanrate || my_kws.initcondfile || \
            my_kws.traj_method==1) {
            cout << "keywords> error: dimensionality reduction simulation not specified correctly" << endl; exit(EXIT_FAILURE); }
    } else if (my_kws.wrapper_method==1) { // standard simulation of A<-B transition paths

    } else if (my_kws.wrapper_method==2) { // WE simulation
        if (my_kws.tau<=0. || (my_kws.commsfile!=nullptr && !my_kws.adaptivecomms) || \
            (my_kws.commstargfile!=nullptr && !my_kws.adaptivecomms) || my_kws.traj_method!=1) {
            cout << "keywords> error: WE simulation not specified correctly" << endl; exit(EXIT_FAILURE); }
    } else if (my_kws.wrapper_method==3) { // FFS simulation
        if (my_kws.commsfile==nullptr) {
            cout << "keywords> error: FFS simulation not specified correctly" << endl; exit(EXIT_FAILURE); }
    } else if (my_kws.wrapper_method==4) { // NEUS simulation
        if (my_kws.commsfile==nullptr || my_kws.traj_method!=1) {
            cout << "keywords> error: NEUS simulation not specified correctly" << endl; exit(EXIT_FAILURE); }
    } else if (my_kws.wrapper_method==5) { // milestoning simulation
        if (my_kws.commsfile==nullptr) {
            cout << "keywords> error: milestoning simulation not specified correctly" << endl; exit(EXIT_FAILURE); }
    }
    // check specification of trajectory method is valid
    if (my_kws.traj_method==1) { // BKL algorithm
        if (!my_kws.branchprobs) {
            cout << "keywords> error: BKL algorithm not specified correctly" << endl; exit(EXIT_FAILURE); }
    } else if (my_kws.traj_method==2) { // kPS algorithm
        if ((!(my_kws.tau>0.) && !my_kws.branchprobs) || (my_kws.commsfile==nullptr && !my_kws.adaptivecomms) || my_kws.nelim<=0 || \
            (my_kws.kpskmcsteps>0 && !my_kws.branchprobs) || (my_kws.pfold && my_kws.ncomms!=3)) {
            cout << "keywords> error: kPS algorithm not specified correctly" << endl; exit(EXIT_FAILURE); }
    } else if (my_kws.traj_method==3) { // MCAMC algorithm

    }

    return my_kws;
}
