/*
Read keywords from file "input.kmc", also read in information on KTN and communities from files
*/

#ifndef __KEYWORDS_H_INCLUDED__
#define __KEYWORDS_H_INCLUDED__

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <typeinfo>

#include <iostream>

using namespace std;

struct Keywords {

    ~Keywords() {
        if (binfile) delete[] binfile;
        if (bintargfile) delete[] bintargfile;
        if (bininitfile) delete[] bininitfile;
    }

    // main keywords (see documentation)
    int n_nodes=0, n_edges=0;
    /* choice of kinetic Monte Carlo method to propagate individual trajectories */
    int kmc_method=-1;
    /* choice of enhanced sampling kinetic Monte Carlo method to accelerate the observation of rare events */
    int enh_method=-1; // note that there is no default, to run kMC with no enhanced sampling this must be set explicitly
    string minafile, minbfile; // names of the files containing the IDs of the A and B nodes, respectively
    int nA=0, nB=0;          // number of nodes in A and B sets, respectively

    // optional arguments pertaining to enhanced sampling methods
    double tau=0.;           // "TAU" time interval between checking bins (WE-kMC)
                             //       lag time at which transition matrix is evaluated (kPS)
    int nbins=-1;            // number of bins (WE-kMC) or trapping basins (kPS) on the network
    char *binfile=nullptr;   // "BINFILE" name of file where bins are defined (WE-kMC, kPS)
    char *bintargfile=nullptr; // "BINTARGFILE" name of file where target number of trajectories in each bin is defined (WE-kMC)
    char *bininitfile=nullptr; // "BININITFILE" name of file where initial prob distrib of trajs is defined (WE-kMC)
    bool adaptivebins=false; // "ADAPTIVEBINS" bins (WE-kMC) or trapping basins (kPS) are determined on-the-fly
    int kpskmcsteps=0;       // "KPSKMCSTEPS" number of BKL kMC steps after a trapping basin escape (kPS)
    int nelim=-1;            // "NELIM" maximum number of states to be eliminated from any trapping basin (kPS)

    // other keywords
    bool debug=false;
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
        }
    }
    inp_f.close();
    return vec_data;
    }

};

#endif
