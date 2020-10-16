/*
Read keywords from file "input.kmc", also read in information on KTN and communities from files

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
#include <omp.h>

using namespace std;

struct Keywords {

    ~Keywords() {
        if (initcondfile) delete[] initcondfile;
        if (commsfile) delete[] commsfile;
        if (commstargfile) delete[] commstargfile;
        if (binsfile) delete[] binsfile;
        if (ntrajsfile) delete[] ntrajsfile;
    }

    /* main keywords (see documentation). Here, -1 represents a value that must be set if the parameter is mandatory given
       the combination of chosen keywords. Values of 0 are default values that are valid in any case */
    int n_nodes=0, n_edges=0;
    /* choice of kinetic Monte Carlo method to propagate individual trajectories */
    int traj_method=-1;
    /* choice of enhanced sampling kinetic Monte Carlo method to wrap simulation of a trajectory ensemble */
    int wrapper_method=-1; // note that there is no default, to run kMC with no enhanced sampling this must be set explicitly
    string nodesafile, nodesbfile; // names of the files containing the IDs of the A and B nodes, respectively
    int nA=0, nB=0;          // number of nodes in A and B sets, respectively

    // optional arguments pertaining to enhanced sampling methods
    long double tau=-1.;     // "TAU" time interval between checking communities (WE-kMC)
                             //       lag time at which transition probability matrix is evaluated (kPS)
    double tintvl=-1.;       // "TINTVL" time interval for writing trajectory data
    long double dt=-1.;      // cf "DIMREDUCTION" time length for all trajectories, initialised from each community in turn
    int ncomms=-1;           // number of communities on the network, eg no. of communities for resampling (WE-kMC) or trapping basins (kPS)
    int nbins=-1;            // number of bins on the network, used to calculate transition path statistics
    int nthreads=omp_get_max_threads(); // number of threads to use in parallel calculations
    char *initcondfile=nullptr; // "INITCOND" name of file where nonequilibrium initial probs of nodes in B are specified
    char *commsfile=nullptr; // "COMMSFILE" name of file where communities are defined (WE-kMC, kPS)
    char *commstargfile=nullptr; // "COMMSTARGFILE" name of file where target number of trajectories in each community is defined (WE-kMC)
    char *binsfile=nullptr;   // "BINSFILE" name of file where bins are defined (for calculating TP statistics)
    char *ntrajsfile=nullptr; // cf "DIMREDUCTION" name of file where number of short trajectories to be ran from each community are defined
    bool adaptivecomms=false; // "ADAPTIVECOMMS" communities for resampling (WE-kMC) or trapping basins (kPS) are determined on-the-fly
    int kpskmcsteps=0;       // "KPSKMCSTEPS" number of BKL kMC steps after a trapping basin escape (kPS)
    int nelim=-1;            // "NELIM" maximum number of states to be eliminated from any trapping basin (kPS)
    double taure=0.;         // "TAURE" time between resampling ensemble of trajectories (WE)
    int nabpaths=-1;         // "NABPATHS" target number of complete A-B paths to simulate
    int maxit=numeric_limits<int>::max(); // "MAXIT" maximum number of iterations of the relevant standard or enhanced kMC algorithm
    double adaptminrate=0.;  // "ADAPTIVECOMMS" minimum transition rate to include in the BFS procedure to define a community on-the-fly

    // keywords for state reduction methods
    bool committor=false;    // "COMMITTOR" specifies that a committor probability calculation is to be performed instead of a kPS simulation
    bool absorption=false;   // "ABSORPTION" specifies that an absorption probability calculation is to be performed
    bool fundamentalred=false; // "FUNDAMENTALRED" specifies that the fundamental matrix of an absorbing (reducible) Markov chain is to be computed
    bool fundamentalirred=false; // "FUNDAMENTALIRRED" specifies that the fundamental matrix of an irreducible Markov chain is to be computed
    bool mfpt=false;         // "MFPT" specifies that the MFPTs for transitions from all non-target nodes are to be computed
    bool gth=false;          // "GTH" specifies that the Grassmann-Taksar-Heyman algorithm for computation of the stationary distribution is to be performed
    bool pathlengths=false;  // "PATHLENGTHS" specifies that mean first passage path lengths (instead of times) are calculated

    // other keywords
    bool transnprobs=false;  // "TRANSNPROBS" edge weights are read in as transition probabilities (not as weights)
    bool discretetime=false; // "DISCRETETIME" edge weights read in as trans probs represent a DTMC at lag time tau
    bool branchprobs=false;  // "BRANCHPROBS" transition probabilities are calculated as branching probabilities
    bool meanrate=false;     // "MEANRATE" use the approximate mean rate method in MCAMC, instead of the exact FPTA method (default)
    bool dumpintvls=false;   // "DUMPINTVLS" trajectory data is dumped at fixed time intervals
    bool debug=false;
    int seed=17;
    bool dumpwaittimes=false; // "DUMPWAITTIMES" print waiting times for nodes to file "meanwaitingtimes.dat"

    // implicitly set switches
    bool initcond=false;     // "INITCOND" specifies if a nonequilibrium initial condition for the nodes in set B has been set
    bool statereduction=false; // is true when the purpose of the computation is to perform a state reduction procedure
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
        } else if (typeid(T)==typeid(long double)) {
            vec_data.emplace_back(make_pair(stold(vecstr[0]),stod(vecstr[1])));
        } else { // inappropriate data type of entries in file
            throw exception();
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
        } else { // inappropriate data type of entries in file
            throw exception();
        }
    }
    inp_f.close();
    return vec_data;
    }

};

#endif
