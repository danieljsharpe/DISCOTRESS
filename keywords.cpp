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
    ifstream kw_f(kw_file);
    while (getline(kw_f,line)) {
        vector<string> vecstr;
        stringstream ss(line);
        while (getline(ss,token,*delim)) { // collect keyword + args in vector
            vecstr.emplace_back(token); }

        // PROCESS KEYWORDS
        if (vecstr[0]=="TRAJ") {
            if (vecstr[1]=="GILLESPIE") {
                my_kws.kmc_method=1;
            } else if (vecstr[1]=="TAULEAPING") {
                my_kws.kmc_method=2;
            } else if (vecstr[2]=="BKL") {
                my_kws.kmc_method=3;
            } else { exit(EXIT_FAILURE); }
        } else if (vecstr[0]=="ENHANCED") {
            if (vecstr[1]=="NONE") {
                my_kws.enh_method=0;
            } else if (vecstr[1]=="WE") {
                my_kws.enh_method=1;
            } else if (vecstr[1]=="FFS") {
                my_kws.enh_method=2;
            } else { exit(EXIT_FAILURE); }
        } else if (vecstr[0]=="TAU") {
            my_kws.tau=stod(vecstr[1]);
        } else if (vecstr[0]=="NWALKERS") {
            my_kws.nwalkers=stoi(vecstr[1]);
        } else if (vecstr[0]=="BINFILE") {
            my_kws.binfile=vecstr[1].c_str();
            if (vecstr.size()>2) { my_kws.startbin=stoi(vecstr[2]); }
        } else {
            cout << "Fatal error: Unrecognised keyword: " << vecstr[0] << endl;
            exit(EXIT_FAILURE);
        }
    }
    kw_f.close();

    // CHECK NECESSARY KEYWORDS
    if (my_kws.kmc_method!=1 && my_kws.kmc_method!=2 && my_kws.kmc_method!=3) {
        cout << "Fatal error: KMC method not chosen correctly" << endl; exit(EXIT_FAILURE); }
    if (my_kws.enh_method!=0 && my_kws.enh_method!=1 && my_kws.enh_method!=2) {
        cout << "Fatal error: Enhanced KMC method not chosen correctly" << endl; exit(EXIT_FAILURE); }
    if (my_kws.enh_method==1) { // WE simulation
        if (!(my_kws.tau>0.) || !(my_kws.nwalkers>0) || !(strlen(my_kws.binfile)>0)) {
            cout << "Fatal error: WE simulation not set up correctly" << endl; exit(EXIT_FAILURE); } }

    return my_kws;
}
