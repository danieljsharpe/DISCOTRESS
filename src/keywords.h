/*
Read keywords from file "input.kmc". Also read in information on Markov chain (KTN) parameters and topology from input files.

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
    int n_nodes=0, n_edges=0; // "NNODES" / "NEDGES" number of nodes / number of (bidirectional) edges in Markov chain
    /* choice of enhanced sampling method to wrap simulation of a trajectory ensemble */
    int wrapper_method=-1;    // "WRAPPER" enhanced sampling method (note that there is no default)
    /* choice of kinetic Monte Carlo method to propagate individual trajectories */
    int traj_method=-1;       // "TRAJ" method (note that there is no default)
    string nodesafile, nodesbfile; // "NODESAFILE" / "NODESBFILE" names of the files containing the IDs of the A and B nodes, respectively
    int nA=0, nB=0;           // "NODESAFILE" / "NODESBFILE" number of nodes in A and B sets, respectively

    // optional keywords relating to simulation parameters and output
    char *binsfile=nullptr;   // "BINSFILE" name of file where bins are defined (for calculating transition path statistics for bins)
    int nbins=-1;             // "BINSFILE" number of bins on the network, used to calculate transition path statistics
    char *commsfile=nullptr;  // "COMMSFILE" name of file where communities are defined
    int ncomms=-1;            // "COMMSFILE" number of communities on the network, eg no. of communities for resampling (WE-kMC) or trapping basins (kPS)
    bool dumpintvls=false;    // "DUMPINTVLS" trajectory data is dumped at fixed time intervals
    char *initcondfile=nullptr; // "INITCOND" name of file where nonequilibrium initial probs of nodes in B are specified
    int maxit=numeric_limits<int>::max(); // "MAXIT" maximum number of iterations of the relevant standard or enhanced kMC algorithm
    int nabpaths=-1;          // "NABPATHS" target number of complete A-B paths to simulate
    double tintvl=-1.;        // "TINTVL" time interval for writing trajectory data

    // optional keywords pertaining to enhanced sampling methods
    bool adaptivecomms=false; // "ADAPTIVECOMMS" communities for resampling (WE-kMC) or trapping basins (kPS) are determined on-the-fly
    double adaptminrate=0.;   // "ADAPTIVECOMMS" minimum transition rate to include in the BFS procedure to define a community on-the-fly
    char *commstargfile=nullptr; // "COMMSTARGFILE" name of file where target number of trajectories in each community is defined (WE-kMC)
    char *ntrajsfile=nullptr; // "DIMREDUCTION" name of file where number of short trajectories to be ran from each community are defined
    int kpskmcsteps=0;        // "KPSKMCSTEPS" number of BKL kMC steps after a trapping basin escape (kPS or MCAMC)
    bool meanrate=false;      // "MEANRATE" use the approximate mean rate method in MCAMC, instead of the exact FPTA method (default)
    int nelim=-1;             // "NELIM" maximum number of states to be eliminated from any trapping basin (kPS)
    int nwalkers=-1;          // "NWALKERS" for certain enhanced sampling (WRAPPER) methods, number of independent trajectories on the network. For
                              //      certain other enhanced sampling methods, this parameter is ignored and overriden to a default value
    bool reanotirred=false;   // "REANOTIRRED" prevents throwing of errors when a candidate path cannot be found in the REA (expected behaviour for
                              //      reducible, but not irreducible, Markov chains)
    bool steadystate=false;   // "STEADYSTATE" indicates that a small number of trajectories are to be used to estimate steady state dynamical properties
    double ssrec=0.;          // "STEADYSTATE" time interval after which the trajectory is considered to be equilibriated and recording of steady state
                              //      A<-B transition path ensemble statistics begins
    double taure=0.;          // "TAURE" time between resampling ensemble of trajectories (WE)
    long double trajt=0.;     // "TRAJT" max time for trajectories (when simulating trajectories of fixed total time)
    bool writerea=false;      // "WRITEREA" if WRAPPER REA, write trajectory data for the k shortest paths to output files

    // keywords for state reduction methods
    bool absorption=false;    // "ABSORPTION" specifies that an absorption probability calculation is to be performed
    bool committor=false;     // "COMMITTOR" specifies that a committor probability calculation is to be performed instead of a kPS simulation
    bool fundamentalirred=false; // "FUNDAMENTALIRRED" specifies that the fundamental matrix of an irreducible Markov chain is to be computed
    bool fundamentalred=false; // "FUNDAMENTALRED" specifies that the fundamental matrix of an absorbing (reducible) Markov chain is to be computed
    bool gth=false;           // "GTH" specifies that the Grassmann-Taksar-Heyman algorithm for computation of the stationary distribution is to be performed
    bool mfpt=false;          // "MFPT" specifies that the MFPTs for transitions from all non-target nodes are to be computed
    bool pathlengths=false;   // "PATHLENGTHS" specifies that mean first passage path lengths (instead of times) are calculated

    // other keywords
    bool accumprobs=false;    // "ACCUMPROBS" if simulating walkers using the BKL algorithm, optimize efficiency by ordering edges by transition probs
    bool branchprobs=false;   // "BRANCHPROBS" transition probabilities are calculated as branching probabilities
    bool debug=false;         // "DEBUG" turn on extra print statements to aid debugging
    bool discretetime=false;  // "DISCRETETIME" edge weights are read in as transition probabilities (instead of log transition rates). The provided
                              //                edge weights therefore represent a discrete-time Markov chain (DTMC) at lag time tau
    bool dumpwaittimes=false; // "DUMPWAITTIMES" print waiting times for nodes to file "meanwaitingtimes.dat"
    bool noloop=false;        // "NOLOOP" (for a DTMC) renormalize lag times for nodes and outgoing transition probabilities to subsume self-loops
    int nthreads=omp_get_max_threads(); // number of threads to use in parallel calculations
    int seed=17;              // "SEED" seed for random number generators
    long double tau=-1.;      // "TAU" lag time (DTMC) or mean waiting time in linearised transition matrix (CTMC if not using branching probabilities)

    // implicitly set switches
    bool initcond=false;      // "INITCOND" specifies if a nonequilibrium initial condition for the nodes in set B has been set
    bool statereduction=false; // is true when the purpose of the computation is to perform a state reduction procedure

    void check_keywords();    // function to check that keyword specification is appropriate
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
