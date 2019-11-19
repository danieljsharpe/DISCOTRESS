#include "keywords.h"
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <iostream>

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
                my_kws.kmc_method=1;
            } else if (vecstr[1]=="REJ") {
                my_kws.kmc_method=2;
            } else if (vecstr[2]=="LEAPFROG") {
                my_kws.kmc_method=3;
            } else { exit(EXIT_FAILURE); }
        } else if (vecstr[0]=="ENHANCED") {
            if (vecstr[1]=="NONE") {
                my_kws.enh_method=0;
            } else if (vecstr[1]=="WE") {
                my_kws.enh_method=1;
            } else if (vecstr[1]=="KPS") {
                my_kws.enh_method=2;
            } else if (vecstr[2]=="FFS") {
                my_kws.enh_method=3;
            } else if (vecstr[3]=="ASKMC") {
                my_kws.enh_method=4;
            } else if (vecstr[4]=="NEUS") {
                my_kws.enh_method=5;
            } else { exit(EXIT_FAILURE); }
        } else if (vecstr[0]=="NNODES") {
            my_kws.n_nodes=stoi(vecstr[1]);
        } else if (vecstr[0]=="NEDGES") {
            my_kws.n_edges=stoi(vecstr[1]);
        } else if (vecstr[0]=="MINAFILE") {
            my_kws.minafile=vecstr[1];
            my_kws.nA=stoi(vecstr[2]);
        } else if (vecstr[0]=="MINBFILE") {
            my_kws.minbfile=vecstr[1];
            my_kws.nB=stoi(vecstr[2]);
        } else if (vecstr[0]=="TAU") {
            my_kws.tau=stod(vecstr[1]);
        } else if (vecstr[0]=="BINFILE") {
            my_kws.binfile = new char[vecstr[1].size()+1];
            copy(vecstr[1].begin(),vecstr[1].end(),my_kws.binfile);
            my_kws.binfile[vecstr[1].size()]='\0'; // trailing character
            my_kws.nbins=stoi(vecstr[2]);
        } else if (vecstr[0]=="BINTARGFILE") {
            my_kws.bintargfile = new char[vecstr[1].size()+1];
            copy(vecstr[1].begin(),vecstr[1].end(),my_kws.bintargfile);
        } else if (vecstr[0]=="BININITFILE") {
            my_kws.bininitfile = new char[vecstr[1].size()+1];
            copy(vecstr[1].begin(),vecstr[1].end(),my_kws.bininitfile);
        } else if (vecstr[0]=="ADAPTIVEBINS") {
            my_kws.adaptivebins=true;
        } else if (vecstr[0]=="KPSKMCSTEPS") {
            my_kws.kpskmcsteps=stoi(vecstr[1]);
        } else if (vecstr[0]=="NELIM") {
            my_kws.nelim=stoi(vecstr[1]);
        } else if (vecstr[0]=="DEBUG") {
            my_kws.debug=true;
        } else {
            cout << "Fatal error: Unrecognised keyword: " << vecstr[0] << endl;
            exit(EXIT_FAILURE);
        }
    }
    kw_f.close();
    cout << ">>>>> finished reading keywords" << endl;

    // check necessary keywords and compatability
    if (my_kws.n_nodes<=0 || my_kws.n_edges<=0 || my_kws.nA<=0 || my_kws.nB<=0) {
        cout << "Fatal error: transition network parameters not set correctly" << endl; exit(EXIT_FAILURE); }
    if (my_kws.kmc_method<=0 && my_kws.enh_method!=2) {
        cout << "Fatal error: must choose a kMC trajectory method except with kPS" << endl; exit(EXIT_FAILURE); }
    if (my_kws.enh_method==-1) {
        cout << "Fatal error: Enhanced kMC method not chosen correctly" << endl; exit(EXIT_FAILURE); }
    if (my_kws.enh_method==1) { // WE simulation
        if (!(my_kws.tau>0.) || (!(strlen(my_kws.binfile)>0) && !my_kws.adaptivebins) || \
            (!(strlen(my_kws.bintargfile)>0) && !my_kws.adaptivebins) ) {
            cout << "Fatal error: WE simulation not set up correctly" << endl; exit(EXIT_FAILURE); } }
    if (my_kws.enh_method==2) { // kPS simulation
        if (!(my_kws.tau>0.) || (!(strlen(my_kws.binfile)>0) && !my_kws.adaptivebins) || my_kws.nelim<=0 ) {
            cout << "Fatal error: kPS simulation not set up correctly" << endl; exit(EXIT_FAILURE); } }

    return my_kws;
}
