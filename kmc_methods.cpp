/*
File containing kinetic Monte Carlo simulation algorithms to propagate individual trajectories
*/

#include "kmc_methods.h"
#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

Walker::~Walker() {}

/* write trajectory and path quantities to file */
void Walker::dump_walker_info(int path_no, bool endofpath) {
    if (curr_node==nullptr) throw exception();
    ofstream walker_f;
    walker_f.open("walker."+to_string(path_no)+".dat",ios_base::app);
    walker_f.setf(ios::right,ios::adjustfield); walker_f.width(5);
    walker_f << curr_node->node_id << "  " << curr_node->comm_id;
    walker_f.precision(16); walker_f.width(12);
    walker_f << "    " << t << "    " << k  << "    " << p << "    "  << s << endl;
    if (!endofpath) return;
    ofstream tpdistrib_f;
    tpdistrib_f.open("tp_distribns.dat",ios_base::app);
    tpdistrib_f.precision(16); tpdistrib_f.width(12);
    tpdistrib_f << path_no << "    " << t << "    " << k << "    " << p << "    " << s << endl;
}

/* reset path quantities */
void Walker::reset_walker_info() {
    k=0; p=-numeric_limits<double>::infinity(); t=0.; s=0.;
    curr_node=nullptr;
}

KMC_Enhanced_Methods::KMC_Enhanced_Methods() {}

KMC_Enhanced_Methods::~KMC_Enhanced_Methods() {}

/* function to set the KMC_Standard_Method member function to propagate individual trajectories */
void KMC_Enhanced_Methods::set_standard_kmc(void(*kmcfuncptr)(Walker&)) {
    kmc_func = kmcfuncptr;
}

void KMC_Enhanced_Methods::find_bin_onthefly() {

}

/* Increment number of A<-B and B<-B paths simulated. If desired, update the vectors containing counts needed to
   calculate transition path statistics for communities */
void KMC_Enhanced_Methods::update_tp_stats(bool abpath, bool update) {
    n_traj++; if (abpath) n_ab++;
    if (update) {
    for (bool comm_visit: visited) {
        if (comm_visit) continue;
    }
    }
    fill(visited.begin(),visited.end(),false);
}

/* calculate the transition path statistics for communities from the observed counts during the simulation */
void KMC_Enhanced_Methods::calc_tp_stats() {
    cout << "kmc_enhanced_methods> calculating transition path statistics" << endl;
    for (int i=0;i<nbins;i++) {
        committors[i] = static_cast<double>(ab_successes[i])/static_cast<double>(ab_successes[i]+ab_failures[i]+1); // quack
        tp_densities[i] = static_cast<double>(tp_visits[i])/static_cast<double>(n_ab);
    }
    write_tp_stats();
}

/* write transition path statistics to file */
void KMC_Enhanced_Methods::write_tp_stats() {
    cout << "kmc_enhanced_methods> writing transition path statistics to file" << endl;
    ofstream tpstats_f;
    tpstats_f.open("tp_stats.dat");
    for (int i=0;i<nbins;i++) {
        tpstats_f << setw(4) << i << setw(8) << tp_visits[i] << setw(8) << ab_successes[i] << setw(8) << ab_failures[i];
        tpstats_f << fixed << setprecision(6);
        tpstats_f << "    " << setw(10) << tp_densities[i] << setw(10) << committors[i] << endl;
    }
}

STD_KMC::STD_KMC(const Network& ktn, int maxn_abpaths, int maxit) {
    this->maxn_abpaths=maxn_abpaths; this->maxit=maxit;
}

STD_KMC::~STD_KMC() {}

/* main loop to drive a standard kMC simulation (no enhanced sampling) */
void STD_KMC::run_enhanced_kmc(const Network &ktn) {
    cout << "std_kmc> beginning standard kMC simulation" << endl;
    Walker walker={walker_id:0,bin_curr:0,bin_prev:0,k:0,active:true,p:0.,t:0.,s:0.}; // this method uses only a single walker
    walker.curr_node = &ktn.nodes[0]; // quack choose properly
    int n_ab=0, n_kmcit=0;
    while ((n_ab<maxn_abpaths) && (n_kmcit<maxit)) { // algorithm terminates after max. no. of iterations of the standard kMC algorithm
        cout << "iteration " << n_kmcit << endl;
        (*kmc_func)(walker);
        if (walker.curr_node->aorb==-1) n_ab++;
        n_kmcit++;
    }
    cout << "std_kmc> walker time: " << walker.t << " activity: " << walker.k << " entropy flow: " << walker.s << " log path prob: " << walker.p << endl;
    cout << "std_kmc> finished standard kMC simulation" << endl;
}

KMC_Standard_Methods::KMC_Standard_Methods() {}
KMC_Standard_Methods::~KMC_Standard_Methods() {}

void KMC_Standard_Methods::bkl(Walker &walker) {
    // propagate trajectory
    double rand_no = KMC_Standard_Methods::rand_unif_met(19); // quack // random number used to select transition
    Edge *edgeptr = walker.curr_node->top_from;
    const Node *old_node = walker.curr_node;
    double prev_p = 0.; // previous accumulated branching probability
    while (edgeptr!=nullptr) { // loop over FROM edges and check random number against accumulated branching probability
        if (edgeptr->t > rand_no) break;
        prev_p = edgeptr->t;
        edgeptr = edgeptr->next_from;
    }
    cout << "  moving to node: " << edgeptr->to_node->node_id << endl;
    walker.curr_node = edgeptr->to_node; // advance trajectory
    // update path quantities
    walker.k++; // dynamical activity (no. of steps)
    walker.p += log(edgeptr->t-prev_p); // log path probability
    walker.t += -(1./exp(old_node->k_esc))*log(KMC_Standard_Methods::rand_unif_met(19)); // sample transition time
    walker.s += edgeptr->rev_edge->k-edgeptr->k; // entropy flow
}

void KMC_Standard_Methods::rejection_kmc(Walker &walker) {

}

void KMC_Standard_Methods::leapfrog(Walker &walker) {

}

/* draw a uniform random number between 0 and 1, used in Metropolis conditions etc. */
double KMC_Standard_Methods::rand_unif_met(int seed) {
    static default_random_engine generator(seed);
    static uniform_real_distribution<double> unif_real_distrib(0.,1.);
    return unif_real_distrib(generator);
}
