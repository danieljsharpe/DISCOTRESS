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
void Walker::dump_walker_info(int path_no, bool transnpath, bool newpath, bool writetraj) {
    if (curr_node==nullptr) throw exception();
    if (writetraj) {
    ofstream walker_f;
    string walker_fname="walker."+to_string(this->walker_id)+"."+to_string(path_no)+".dat";
    if (!newpath) { walker_f.open(walker_fname,ios_base::app);
    } else { walker_f.open(walker_fname,ios_base::trunc); }
    walker_f.setf(ios::right,ios::adjustfield); walker_f.setf(ios::fixed,ios::floatfield); // walker_f.fill('x');
    walker_f << setw(7) << curr_node->node_id << setw(8) << curr_node->comm_id << setw(15) << k;
    walker_f.precision(12); // walker_f.width(18);
    walker_f << setw(30) << t << setw(30) << p << setw(30) << s << endl;
    }
    if (!transnpath) return;
    ofstream tpdistrib_f;
    tpdistrib_f.open("tp_distribns.dat",ios_base::app);
    tpdistrib_f.setf(ios::right,ios::adjustfield); tpdistrib_f.setf(ios::fixed,ios::floatfield);
    tpdistrib_f.precision(12);
    tpdistrib_f << setw(15) << path_no << setw(15) << k << setw(30) << t << setw(30) << p << setw(30) << s << endl;
}

/* reset path quantities */
void Walker::reset_walker_info() {
    k=0; p=-numeric_limits<double>::infinity(); t=0.; s=0.;
    curr_node=nullptr;
}

KMC_Enhanced_Methods::KMC_Enhanced_Methods() {}

KMC_Enhanced_Methods::~KMC_Enhanced_Methods() {}

/* sample an initial node (from the B set) and set this node as the starting node of the walker */
Node *KMC_Enhanced_Methods::get_initial_node(const Network &ktn, Walker &walker) {

    Node *node_b=nullptr; // sampled starting node
    double pi_B = -numeric_limits<double>::infinity(); // (log) occupation probability of all nodes in initial set B
    vector<pair<Node*,double>> b_probs(ktn.nodesB.size()); // accumulated probs of selecting starting node
    set<Node*>::iterator it_set = ktn.nodesB.begin();
    if (ktn.nodesB.size()==1) { // there is only one node in the starting set
        node_b=*it_set;
        pi_B=(*it_set)->pi;
    } else if (!ktn.initcond) { // no initial condition was set, choose node in set B in proportion to stationary probs
        while (it_set!=ktn.nodesB.end()) {
            pi_B = log(exp(pi_B)+exp((*it_set)->pi));
            it_set++; }
        it_set = ktn.nodesB.begin();
        double cum_prob=0.;
        while (it_set!=ktn.nodesB.end()) {
            cum_prob += exp((*it_set)->pi-pi_B);
            b_probs.push_back(make_pair((*it_set),cum_prob));
            it_set++; }
    } else { // choose node in set B in proportion to specified initial condition probs
        pi_B=0.; // for specified initial condition, sum of probabilities should be unity
        double cum_prob=0.; int i=0;
        while (it_set!=ktn.nodesB.end()) {
            cum_prob += exp(ktn.init_probs[i]);
            b_probs.push_back(make_pair(*it_set,cum_prob)); i++;
        }
    }
    if (node_b==nullptr) { // if there was more than one node in B, sample the initial node
        double rand_no = KMC_Standard_Methods::rand_unif_met(seed);
        vector<pair<Node*,double>>::iterator it_vec = b_probs.begin();
        while (it_vec!=b_probs.end()) {
            if ((*it_vec).second>=rand_no) { node_b=(*it_vec).first; break; }
            it_vec++; }
        if (it_vec==b_probs.end()) node_b=(*it_vec).first;
    }
    walker.curr_node=&(*node_b);
    walker.p=node_b->pi-pi_B; // factor in path probability corresponding to initial occupation of node
    if (ktn.ncomms>0) visited[node_b->comm_id]=true;
    return node_b;
}

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
    if (!update) return;
    int i=0; // community ID
    for (bool comm_visit: visited) {
        if (comm_visit) {
            if (abpath) { ab_successes[i]++;
            } else { ab_failures[i]++; }
        }
        i++;
    }
    fill(visited.begin(),visited.end(),false);
}

/* calculate the transition path statistics for communities from the observed counts during the simulation */
void KMC_Enhanced_Methods::calc_tp_stats(int ncomms) {
    cout << "kmc_enhanced_methods> calculating transition path statistics" << endl;
    for (int i=0;i<ncomms;i++) {
        committors[i] = static_cast<double>(ab_successes[i])/static_cast<double>(ab_successes[i]+ab_failures[i]);
        tp_densities[i] = static_cast<double>(ab_successes[i])/static_cast<double>(n_ab);
    }
    write_tp_stats(ncomms);
}

/* write transition path statistics to file */
void KMC_Enhanced_Methods::write_tp_stats(int ncomms) {
    cout << "kmc_enhanced_methods> writing transition path statistics to file" << endl;
    ofstream tpstats_f;
    tpstats_f.open("tp_stats.dat");
    for (int i=0;i<ncomms;i++) {
        tpstats_f << setw(7) << i << setw(10) << ab_successes[i] << setw(10) << ab_failures[i];
        tpstats_f << fixed << setprecision(6);
        tpstats_f << setw(20) << tp_densities[i] << setw(14) << committors[i] << endl;
    }
}

STD_KMC::STD_KMC(const Network& ktn, int maxn_abpaths, int maxit, double tintvl, int seed) {
    // quack move this somewhere more general
    this->maxn_abpaths=maxn_abpaths; this->maxit=maxit; this->tintvl=tintvl; this->seed=seed;
    if (ktn.ncomms>0) {
        visited.resize(ktn.ncomms); fill(visited.begin(),visited.end(),false);
        tp_densities.resize(ktn.ncomms); committors.resize(ktn.ncomms);
        ab_successes.resize(ktn.ncomms); ab_failures.resize(ktn.ncomms);
    }
}

STD_KMC::~STD_KMC() {}

/* main loop to drive a standard kMC simulation (no enhanced sampling) */
void STD_KMC::run_enhanced_kmc(const Network &ktn) {
    cout << "\nstd_kmc> beginning standard kMC simulation" << endl;
    Walker walker={walker_id:0,bin_curr:0,bin_prev:0,k:0,active:true,p:0.,t:0.,s:0.}; // only a single walker
    double dummy_randno = KMC_Standard_Methods::rand_unif_met(seed); // seed generator
    Node *dummy_node = get_initial_node(ktn,walker);
    walker.dump_walker_info(0,false,true,tintvl>=0.);
    n_ab=0; n_traj=0; int n_kmcit=1;
    double next_tintvl=tintvl; // next time interval for writing trajectory data
    bool leftb=false; // flag indicates if trajectory has left initial set of states yet
    while ((n_ab<maxn_abpaths) && (n_kmcit<=maxit)) { // algo terminates after max no of iterations of the standard kMC algorithm
        (*kmc_func)(walker);
        visited[walker.curr_node->comm_id]=true;
        if (!leftb && walker.curr_node->aorb!=1) leftb=true;
        walker.dump_walker_info(n_ab,walker.curr_node->aorb==-1,false, \
            (walker.curr_node->aorb==-1 || tintvl==0. || (tintvl>0. && walker.t>=next_tintvl)));
        if (tintvl>0. && walker.t>=next_tintvl) { // reached time interval for dumping trajectory data, calc next interval
            while (walker.t>=next_tintvl) next_tintvl+=tintvl; }
        n_kmcit++;
        if (walker.curr_node->aorb==-1 || (walker.curr_node->aorb==1 && leftb)) {
            // traj has reached absorbing macrostate A or has returned to B
            update_tp_stats(walker.curr_node->aorb==-1,ktn.ncomms>0);
            walker.reset_walker_info();
            dummy_node = get_initial_node(ktn,walker);
            if ((n_ab<maxn_abpaths) && (n_kmcit<maxit) && (tintvl>=0.)) walker.dump_walker_info(n_ab,false,true);
            leftb=false; next_tintvl=tintvl;
        }
    }
    cout << "std_kmc> standard kMC simulation terminated after " << n_kmcit-1 << " kMC iterations. Simulated " \
         << n_ab << " transition paths" << endl;
    if (ktn.ncomms>0) calc_tp_stats(ktn.ncomms);
}

KMC_Standard_Methods::KMC_Standard_Methods() {}
KMC_Standard_Methods::~KMC_Standard_Methods() {}

/* function to take a single kMC step using the BKL algorithm */
void KMC_Standard_Methods::bkl(Walker &walker) {
    // propagate trajectory
    double rand_no = KMC_Standard_Methods::rand_unif_met(); // quack // random number used to select transition
    Edge *edgeptr = walker.curr_node->top_from;
    const Node *old_node = walker.curr_node;
    double prev_p = 0.; // previous accumulated branching probability
    while (edgeptr!=nullptr) { // loop over FROM edges and check random number against accumulated branching probability
        if (edgeptr->t > rand_no) break;
        prev_p = edgeptr->t;
        edgeptr = edgeptr->next_from;
    }
    walker.curr_node = edgeptr->to_node; // advance trajectory
    // update path quantities
    walker.k++; // dynamical activity (no. of steps)
    walker.p += log(edgeptr->t-prev_p); // log path probability
    walker.t += -(1./exp(old_node->k_esc))*log(KMC_Standard_Methods::rand_unif_met()); // quack // sample transition time
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
