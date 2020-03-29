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
            } else if (vecstr[5]=="MILES") {
                my_kws.enh_method=6;
            } else if (vecstr[6]=="TPS") {
                my_kws.enh_method=7;
            } else { exit(EXIT_FAILURE); }
        } else if (vecstr[0]=="NNODES") {
            my_kws.n_nodes=stoi(vecstr[1]);
        } else if (vecstr[0]=="NEDGES") {
            my_kws.n_edges=stoi(vecstr[1]);
        } else if (vecstr[0]=="NABPATHS") {
            my_kws.nabpaths=stoi(vecstr[1]);
        } else if (vecstr[0]=="MAXIT") {
            my_kws.maxit=stoi(vecstr[1]);
        } else if (vecstr[0]=="MINAFILE") {
            my_kws.minafile=vecstr[1];
            my_kws.nA=stoi(vecstr[2]);
        } else if (vecstr[0]=="MINBFILE") {
            my_kws.minbfile=vecstr[1];
            my_kws.nB=stoi(vecstr[2]);
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
        } else if (vecstr[0]=="NELIM") {
            my_kws.nelim=stoi(vecstr[1]);
        } else if (vecstr[0]=="PFOLD") {
            my_kws.pfold=true;
        } else if (vecstr[0]=="TRANSNPROBS") {
            my_kws.transnprobs=true;
        } else if (vecstr[0]=="BRANCHPROBS") {
            my_kws.branchprobs=true;
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
    if (my_kws.n_nodes<=0 || my_kws.n_edges<=0 || my_kws.nA<=0 || my_kws.nB<=0) {
        cout << "keywords> error: transition network parameters not set correctly" << endl; exit(EXIT_FAILURE); }
    if (my_kws.nabpaths<=0 || my_kws.maxit<=0) {
        cout << "keywords> error: termination condition not specified correctly" << endl; exit(EXIT_FAILURE); }
    if (my_kws.commsfile!=nullptr && my_kws.ncomms<=1) {
        cout << "keywords> error: there must be at least two communities in the specified partitioning" << endl; exit(EXIT_FAILURE); }
    if (my_kws.kmc_method<=0 && my_kws.enh_method!=2) {
        cout << "keywords> error: must choose a kMC trajectory method except with kPS" << endl; exit(EXIT_FAILURE); }
    if (my_kws.enh_method==-1) {
        cout << "keywords> error: Enhanced kMC method not chosen correctly" << endl; exit(EXIT_FAILURE); }
    if (my_kws.transnprobs && my_kws.enh_method!=2) {
        cout << "keywords> error: edge weights must be read in as transition rates if not using kPS" << endl; exit(EXIT_FAILURE); }
    if (my_kws.enh_method==1) { // WE simulation
        if (!(my_kws.tau>0.) || (my_kws.commsfile!=nullptr && !my_kws.adaptivecomms) || \
            (my_kws.commstargfile!=nullptr && !my_kws.adaptivecomms) ) {
            cout << "keywords> error: WE simulation not set up correctly" << endl; exit(EXIT_FAILURE); } }
    if (my_kws.enh_method==2) { // kPS simulation
        if ((!(my_kws.tau>0.) && !my_kws.branchprobs) || (my_kws.commsfile==nullptr && !my_kws.adaptivecomms) || my_kws.nelim<=0 || \
            (my_kws.kpskmcsteps>0 && !my_kws.branchprobs) || (my_kws.pfold && my_kws.ncomms!=3)) {
            cout << "keywords> error: kPS simulation not set up correctly" << endl; exit(EXIT_FAILURE); } }

    return my_kws;
}
