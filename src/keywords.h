/*
Read keywords from file "input.kmc", also read in information on KTN and communities from files

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

#ifndef __KEYWORDS_H_INCLUDED__
#define __KEYWORDS_H_INCLUDED__

#include <vector>
#include <string>
#include <limits>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <typeinfo>
#include <iostream>

using namespace std;

struct Keywords {

    ~Keywords() {
        if (initcondfile) delete[] initcondfile;
        if (commsfile) delete[] commsfile;
        if (commstargfile) delete[] commstargfile;
        if (binfile) delete[] binfile;
    }

    /* main keywords (see documentation). Here, -1 represents a value that must be set if the parameter is mandatory given
       the combination of chosen keywords. Values of 0 are default values that are valid in any case */
    int n_nodes=0, n_edges=0;
    /* choice of kinetic Monte Carlo method to propagate individual trajectories */
    int kmc_method=-1;
    /* choice of enhanced sampling kinetic Monte Carlo method to accelerate the observation of rare events */
    int enh_method=-1; // note that there is no default, to run kMC with no enhanced sampling this must be set explicitly
    string nodesafile, nodesbfile; // names of the files containing the IDs of the A and B nodes, respectively
    int nA=-1, nB=-1;          // number of nodes in A and B sets, respectively

    // optional arguments pertaining to enhanced sampling methods
    long double tau=-1.;     // "TAU" time interval between checking communities (WE-kMC)
                             //       lag time at which transition probability matrix is evaluated (kPS)
    double tintvl=-1.;       // "TINTVL" time interval for writing trajectory data
    int ncomms=-1;           // number of communities on the network, eg no. of communities for resampling (WE-kMC) or trapping basins (kPS)
    int nbins=-1;            // number of bins on the network, used to calculate transition path statistics
    char *initcondfile=nullptr; // "INITCOND" name of file where nonequilibrium initial probs of nodes in B are specified
    char *commsfile=nullptr; // "COMMSFILE" name of file where communities are defined (WE-kMC, kPS)
    char *commstargfile=nullptr; // "COMMSTARGFILE" name of file where target number of trajectories in each community is defined (WE-kMC)
    char *binfile=nullptr;   // "BINFILE" name of file where bins are defined (for calculating TP statistics)
    bool adaptivecomms=false; // "ADAPTIVECOMMS" communities for resampling (WE-kMC) or trapping basins (kPS) are determined on-the-fly
    int kpskmcsteps=0;       // "KPSKMCSTEPS" number of BKL kMC steps after a trapping basin escape (kPS)
    int nelim=-1;            // "NELIM" maximum number of states to be eliminated from any trapping basin (kPS)
    int nabpaths=-1;         // "NABPATHS" target number of complete A-B paths to simulate
    int maxit=numeric_limits<int>::max(); // "MAXIT" maximum number of iterations of the relevant standard or enhanced kMC algorithm
    double adaptminrate=0.;  // "ADAPTIVECOMMS" minimum transition rate to include in the BFS procedure to define a community on-the-fly

    // other keywords
    bool transnprobs=false;  // "TRANSNPROBS" edge weights are read in as transition probabilities (not as weights)
    bool branchprobs=false;  // "BRANCHPROBS" transition probabilities are calculated as branching probabilities
    bool pfold=false;        // "PFOLD" specifies that a committor function calculation is to be performed instead of a kPS simulation
    bool debug=false;
    int seed=17;
    bool dumpwaittimes=false; // "DUMPWAITTIMES" print waiting times for nodes to file "meanwaitingtimes.dat"

    // implicitly set switches
    bool initcond=false;     // "INITCOND" specifies if a nonequilibrium initial condition for the nodes in set B has been set
};

Keywords read_keywords(const char *);

class Read_files {

    public:

    // read a two-column file
    template <typename T>
    static vector<pair<T,T>> read_two_col(const char *inpfname) {

    string line;
    ifstream inp_f;
    if (!ifstream(inpfname).good()) throw exception(); // check file exists
    inp_f.open(inpfname);
    vector<pair<T,T>> vec_data;
    while (getline(inp_f,line)) {
        vector<string> vecstr;
        istringstream iss(line);
        copy(istream_iterator<string>(iss),istream_iterator<string>(),back_inserter(vecstr));
        if (vecstr.size()!=2) { exit(EXIT_FAILURE); }
        if (typeid(T)==typeid(int)) {
            vec_data.emplace_back(make_pair(stoi(vecstr[0]),stoi(vecstr[1])));
        } else if (typeid(T)==typeid(double)) {
            vec_data.emplace_back(make_pair(stod(vecstr[0]),stod(vecstr[1])));
        }
    }
    inp_f.close();
    return vec_data;
    }

    // read a one-column file
    template <typename T>
    static vector<T> read_one_col(const char *inpfname) {

    string line;
    vector<T> vec_data;
    ifstream inp_f;
    if (!ifstream(inpfname).good()) throw exception(); // check file exists
    inp_f.open(inpfname);
    while (getline(inp_f,line)) {
        if (typeid(T)==typeid(int)) {
            vec_data.emplace_back(stoi(line));
        } else if (typeid(T)==typeid(double)) {
            vec_data.emplace_back(stod(line));
        } else if (typeid(T)==typeid(long double)) {
            vec_data.emplace_back(stold(line));
        }
    }
    inp_f.close();
    return vec_data;
    }

};

#endif
