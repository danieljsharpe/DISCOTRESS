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

using namespace std;

struct Keywords {

    /* choice of kinetic Monte Carlo method to propagate trajectories:
       1. "GILLESPIE" Gillespie / stochastic simulation algorithm (SSA)
       2. "TAULEAPING" tau-leaping algorithm
       3. "BKL" Bortz-Kalos-Lebowitz rejection-free algorithm */
    int kmc_method=-1;
    /* choice of enhanced kinetic Monte Carlo method to maintain trajectories:
       0. "NONE" unbiased KMC
       1. "WE" Weighted ensemble KMC
       2. "FFS" Forward-flux sampling KMC */
    int enh_method=-1;

    // optional arguments for WE-KMC
    double tau=0.;           // "TAU" time interval between checking bins
    int nwalkers=0;          // "NWALKERS" number of walkers maintained on the network
    const char *binfile="";  // "BINFILE" (arg 1.) name of file where bins are defined
    int startbin=-1;         // "BINFILE" (arg 2., opt) index of starting bin for initial trajectory
    
};

Keywords read_keywords(const char *);

class Read_files {

    public:

    // read a two-column file
    template <typename T>
    static vector<pair<T,T>> read_two_col(const char *inpfname) {

    string line;
    ifstream inp_f;
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
