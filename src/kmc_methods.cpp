/*
File containing kinetic Monte Carlo simulation algorithms to propagate individual trajectories

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

#include "kmc_methods.h"
#include <random>
#include <queue>
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
    walker_f << setw(7) << curr_node->node_id << setw(7) << curr_node->comm_id << setw(30) << k;
    walker_f.precision(10); // walker_f.width(18);
    walker_f << setw(60) << t << setw(35) << p << setw(20) << s << endl;
    }
    if (!transnpath) return;
    ofstream tpdistrib_f;
    tpdistrib_f.open("tp_distribns.dat",ios_base::app);
    tpdistrib_f.setf(ios::right,ios::adjustfield); tpdistrib_f.setf(ios::fixed,ios::floatfield);
    tpdistrib_f.precision(10);
    tpdistrib_f << setw(14) << path_no << setw(30) << k << setw(60) << t << setw(35) << p << setw(20) << s << endl;
}

/* reset path quantities */
void Walker::reset_walker_info() {
    k=0; p=-numeric_limits<long double>::infinity(); t=0.; s=0.;
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
            cum_prob += ktn.init_probs[i];
            b_probs.push_back(make_pair(*it_set,cum_prob));
            i++; it_set++;
        }
    }
    if (node_b==nullptr) { // if there was more than one node in B, sample the initial node
        double rand_no = KMC_Enhanced_Methods::rand_unif_met(seed);
        vector<pair<Node*,double>>::iterator it_vec = b_probs.begin();
        while (it_vec!=b_probs.end()) {
            if ((*it_vec).second>=rand_no) { node_b=(*it_vec).first; break; }
            it_vec++; }
        if (it_vec==b_probs.end()) node_b=(*it_vec).first;
    }
    if (node_b==nullptr) throw exception();
    walker.curr_node=&(*node_b);
    walker.p=node_b->pi-pi_B; // factor in path probability corresponding to initial occupation of node
    if (ktn.nbins>0) walker.visited[node_b->bin_id]=true;
    return node_b;
}

/* function to set the KMC_Standard_Method member function to propagate individual trajectories */
void KMC_Enhanced_Methods::set_standard_kmc(void(*kmcfuncptr)(Walker&)) {
    kmc_func = kmcfuncptr;
}

/* use breadth-first search (BFS) procedure to find a community on-the-fly, based on a maximum size of the community
   and a specified transition rate cutoff */
vector<int> KMC_Enhanced_Methods::find_comm_onthefly(const Network &ktn, const Node *init_node, \
        double adaptminrate, int maxsz) {

    vector<int> nodes_in_comm(ktn.n_nodes); // store flags indicating if node is of community or is part of absorbing boundary
    queue<int> nbr_queue; // queue of node IDs to visit in the BFS procedure
    nbr_queue.push(init_node->node_id);
    int nv=0; // number of nodes in the community being built up
    while (!nbr_queue.empty() && nv<maxsz) {
        int curr_node_id = nbr_queue.front();
        nbr_queue.pop();
        nodes_in_comm[curr_node_id-1]=2; nv++; // indicates that node is part of the current community
        const Node *nodeptr = &ktn.nodes[curr_node_id-1];
        const Edge *edgeptr = nodeptr->top_from;
        while (edgeptr!=nullptr) {
            if (edgeptr->deadts || nodes_in_comm[edgeptr->to_node->node_id-1]==2) { // removed edge or node already in comm
                edgeptr=edgeptr->next_from; continue; }
            if (exp(edgeptr->k)>adaptminrate && edgeptr->to_node->aorb!=-1) { // queue neighbouring node to be added into community
                if (nodes_in_comm[edgeptr->to_node->node_id-1]==0) { // node is not already queued
                    nbr_queue.push(edgeptr->to_node->node_id);
                }
            }
            // mark node as belonging to absorbing boundary (for now)
            nodes_in_comm[edgeptr->to_node->node_id-1]=3;
            edgeptr=edgeptr->next_from;
        }
    }
    return nodes_in_comm;
}

/* Increment number of A<-B and B<-B paths simulated. If desired, update the vectors containing counts needed to
   calculate transition path statistics for bins */
void KMC_Enhanced_Methods::update_tp_stats(Walker &walker, bool abpath, bool update) {
    n_traj++; if (abpath) n_ab++;
    if (!update) return;
    int i=0; // bin ID
    for (bool bin_visit: walker.visited) {
        if (bin_visit) {
            if (abpath) { ab_successes[i]++;
            } else { ab_failures[i]++; }
        }
        i++;
    }
    fill(walker.visited.begin(),walker.visited.end(),false);
}

/* calculate the transition path statistics for bins from the observed counts during the simulation */
void KMC_Enhanced_Methods::calc_tp_stats(int nbins) {
    cout << "kmc_enhanced_methods> calculating transition path statistics" << endl;
    for (int i=0;i<nbins;i++) {
        committors[i] = static_cast<double>(ab_successes[i])/static_cast<double>(ab_successes[i]+ab_failures[i]);
        tp_densities[i] = static_cast<double>(ab_successes[i])/static_cast<double>(n_ab);
    }
    write_tp_stats(nbins);
}

/* write transition path statistics to file */
void KMC_Enhanced_Methods::write_tp_stats(int nbins) {
    cout << "kmc_enhanced_methods> writing transition path statistics to file" << endl;
    ofstream tpstats_f;
    tpstats_f.open("tp_stats.dat");
    for (int i=0;i<nbins;i++) {
        tpstats_f << setw(7) << i << setw(10) << ab_successes[i] << setw(10) << ab_failures[i];
        tpstats_f << fixed << setprecision(12);
        tpstats_f << setw(26) << tp_densities[i] << setw(20) << committors[i] << endl;
    }
}

/* draw a uniform random number between 0 and 1, used in Metropolis conditions etc. */
long double KMC_Enhanced_Methods::rand_unif_met(int seed) {
    static default_random_engine generator(seed);
    static uniform_real_distribution<long double> unif_real_distrib(0.,1.);
    return unif_real_distrib(generator);
}

STD_KMC::STD_KMC(const Network& ktn, int maxn_abpaths, int maxit, double tintvl, int seed) {
    // quack move this somewhere more general
    this->walker.accumprobs=ktn.accumprobs;
    this->maxn_abpaths=maxn_abpaths; this->maxit=maxit; this->tintvl=tintvl; this->seed=seed;
    if (ktn.ncomms>0) {
        walker.visited.resize(ktn.nbins); fill(walker.visited.begin(),walker.visited.end(),false);
        tp_densities.resize(ktn.nbins); committors.resize(ktn.nbins);
        ab_successes.resize(ktn.nbins); ab_failures.resize(ktn.nbins);
    }
}

STD_KMC::~STD_KMC() {}

/* main loop to drive a standard kMC simulation (no enhanced sampling) */
void STD_KMC::run_enhanced_kmc(const Network &ktn) {
    cout << "\nstd_kmc> beginning standard kMC simulation" << endl;
    long double dummy_randno = KMC_Enhanced_Methods::rand_unif_met(seed); // seed generator
    Node *dummy_node = get_initial_node(ktn,walker);
    walker.dump_walker_info(0,false,true,tintvl>=0.);
    n_ab=0; n_traj=0; int n_kmcit=1;
    double next_tintvl=tintvl; // next time interval for writing trajectory data
    bool leftb=false; // flag indicates if trajectory has left initial set of nodes yet
    while ((n_ab<maxn_abpaths) && (n_kmcit<=maxit)) { // algo terminates after max no of iterations of the standard kMC algorithm
        (*kmc_func)(walker);
        if (ktn.nbins>0) walker.visited[walker.curr_node->bin_id]=true;
        if (!leftb && walker.curr_node->aorb!=1) leftb=true;
        walker.dump_walker_info(n_ab,walker.curr_node->aorb==-1,false, \
            (walker.curr_node->aorb==-1 || tintvl==0. || (tintvl>0. && walker.t>=next_tintvl)));
        if (tintvl>0. && walker.t>=next_tintvl) { // reached time interval for dumping trajectory data, calc next interval
            while (walker.t>=next_tintvl) next_tintvl+=tintvl; }
        n_kmcit++;
        if (walker.curr_node->aorb==-1 || (walker.curr_node->aorb==1 && leftb)) {
            // traj has reached absorbing macrostate A or has returned to B
            update_tp_stats(walker,walker.curr_node->aorb==-1,ktn.nbins>0);
            if (walker.curr_node->aorb==1) continue;
            walker.reset_walker_info();
            dummy_node = get_initial_node(ktn,walker);
            if ((n_ab<maxn_abpaths) && (n_kmcit<maxit) && (tintvl>=0.)) walker.dump_walker_info(n_ab,false,true);
            leftb=false; next_tintvl=tintvl;
        }
    }
    cout << "std_kmc> standard kMC simulation terminated after " << n_kmcit-1 << " kMC iterations. Simulated " \
         << n_ab << " transition paths" << endl;
    if (ktn.nbins>0) calc_tp_stats(ktn.nbins);
}

KMC_Standard_Methods::KMC_Standard_Methods() {}
KMC_Standard_Methods::~KMC_Standard_Methods() {}

/* function to take a single kMC step using the BKL algorithm */
void KMC_Standard_Methods::bkl(Walker &walker) {
    // propagate trajectory
    double rand_no = KMC_Enhanced_Methods::rand_unif_met(); // random number used to select transition
    Edge *edgeptr = walker.curr_node->top_from;
    const Node *old_node = walker.curr_node;
    long double prev_cum_p = 0.; // previous accumulated branching probability
    long double p; // branching probability of accepted move
    while (edgeptr!=nullptr) { // loop over FROM edges and check random number against accumulated branching probability
        if (walker.accumprobs) { // branching probability values are cumulative
            if (edgeptr->t > rand_no) { p = edgeptr->t-prev_cum_p; break; }
            prev_cum_p = edgeptr->t;
        } else {
            if (edgeptr->t+prev_cum_p > rand_no) { p = edgeptr->t; break; }
            prev_cum_p += edgeptr->t;
        }
        edgeptr = edgeptr->next_from;
    }
    walker.curr_node = edgeptr->to_node; // advance trajectory
    // update path quantities
    walker.k++; // dynamical activity (no. of steps)
    walker.p += log(p); // log path probability
    walker.t += -(1./exp(old_node->k_esc))*log(KMC_Enhanced_Methods::rand_unif_met()); // sample transition time
    walker.s += edgeptr->rev_edge->k-edgeptr->k; // entropy flow
}
