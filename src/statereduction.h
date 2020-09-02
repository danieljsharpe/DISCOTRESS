/*
File containing functions handling the processing and output of the state reduction methods


This file is a part of DISCOTRESS, a software package to simulate the dynamics on arbitrary continuous- and discrete-time Markov chains (CTMCs).
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

#ifndef __STATEREDUCTION_H_INCLUDED__
#define __STATEREDUCTION_H_INCLUDED__

#include <cmath>
#include <string>

using namespace std;

/* calculate the A<-B and B<-A committor probabilities and write to files */
void KPS::calc_committor(const Network& ktn) {

    cout << "kps> calculating committor probabilities from the graph transformation" << endl;
    vector<long double> q_ab_vals(ktn.n_nodes); // A<-B committor
    vector<long double> q_ba_vals(ktn.n_nodes); // B<-A committor
    /* nodes in A that are not directly connected to the intermediate set I do not appear in the ktn_kps Network object -
       use a nodemask to do bookkeeping. The A<-B committor probability for any node in A is 1, and the A<-B committor
       probability of any *internal* node in B is 0. The committor probability of any initial node at the *boundary*
       of B is computed from the committor probabilities of neighbouring nodes in the set I, and the transition
       probabilities to these nodes. The same logic applies to the A<-B committor probabilities */
    vector<bool> nodemask(ktn.n_nodes,false);
    /* committor probabilities for intermediate nodes */
    for (const Node& node: ktn_kps->nodes) {
        nodemask[node.node_id-1]=true;
        if (node.aorb==-1 || node.aorb==1) continue; // node in A or B
        const Edge *edgeptr = node.top_from;
        long double q_ab=0., q_ba=0.;
        while (edgeptr!=nullptr) {
            if (edgeptr->deadts || edgeptr->to_node->eliminated) {
                edgeptr=edgeptr->next_from; continue; }
            if (edgeptr->to_node->aorb==-1) { // node in A
                q_ab += edgeptr->t;
            } else if (edgeptr->to_node->aorb==1) { // node in B
                q_ba += edgeptr->t;
            } else { // all remaining edges should be to nodes in A or B
                throw exception();
            }
            edgeptr=edgeptr->next_from;
        }
        q_ab_vals[node.node_id-1] = q_ab;
        q_ba_vals[node.node_id-1] = q_ba;
    }
    /* committor probabilities for endpoint nodes */
    for (int i=0;i<ktn.n_nodes;i++) {
        if (nodemask[i] && ktn.nodes[i].aorb==0) { continue; // intermediate nodes have all been accounted for
        } else if (!nodemask[i] && ktn.nodes[i].aorb==-1) { // internal node of A
            q_ab_vals[ktn.nodes[i].node_id-1]=1.; q_ba_vals[ktn.nodes[i].node_id-1]=0.;
        } else if (nodemask[i] && ktn.nodes[i].aorb==-1) { // boundary node of A
            q_ab_vals[ktn.nodes[i].node_id-1]=1.;
            q_ba_vals[ktn.nodes[i].node_id-1] = KPS::committor_boundary_node(ktn,ktn.nodes[i].node_id,q_ba_vals,-1);
        } else if (nodemask[i] && ktn.nodes[i].aorb==1) { // boundary or internal node of B
            q_ba_vals[ktn.nodes[i].node_id-1]=1.;
            q_ab_vals[ktn.nodes[i].node_id-1] = KPS::committor_boundary_node(ktn,ktn.nodes[i].node_id,q_ab_vals,1);
        } else { // the nodemask should have tracked all nodes except the internal nodes of A
            throw exception();
        }
    }
    Wrapper_Method::write_vec<long double>(q_ab_vals,"committor_AB.dat");
    Wrapper_Method::write_vec<long double>(q_ba_vals,"committor_BA.dat");
    cout << "kps> finished writing committor probabilities to files" << endl;
}

/* calculate the committor probability for an initial node at the boundary of the initial state, which is /= 0 */
long double KPS::committor_boundary_node(const Network& ktn, int node_id, const vector<long double> q_vals, int aorb) {
    const Node& node = ktn.nodes[node_id-1];
    long double q_val=0.;
    const Edge *edgeptr = node.top_from;
    while (edgeptr!=nullptr) {
        if (!(edgeptr->deadts || edgeptr->to_node->aorb==aorb)) {
            q_val += edgeptr->t*q_vals[edgeptr->to_node->node_id-1]; }
        edgeptr=edgeptr->next_from;
    }
    return q_val;
}

/* compute and write absorption probabilities. NB At this point, pi values of the ktn_kps object should, for
   initial nodes, have been overwritten to the initial probability distribution values */
void KPS::calc_absprobs() {
    cout << "kps> calculating absorption probabilities from the graph transformation" << endl;
    write_renormalised_probs("absorption.dat"); // compute probabilities that trajectories initialised at node i are absorbed at node j
    if (gth || ktn_kps->nodesA.size()==1) return; // stat probs are not known at this point, or only one possible target state to hit
    /* find total hitting probabilities for absorbing nodes given the initial probability distribution */
    ofstream hitprob_f; hitprob_f.open("hitting_probs.dat");
    hitprob_f.setf(ios::right,ios::adjustfield); hitprob_f.setf(ios::scientific,ios::floatfield);
    hitprob_f.precision(10);
    vector<Node>::iterator it_nodevec = ktn_kps->nodes.begin();
    while (it_nodevec!=ktn_kps->nodes.end()) {
        if (it_nodevec->aorb!=-1) { it_nodevec++; continue; } // not an absorbing node
        long double b=0.; // hitting probability for this absorbing node
        Edge *edgeptr = it_nodevec->top_to;
        while (edgeptr!=nullptr) {
            if (edgeptr->deadts || edgeptr->from_node->aorb!=1) { edgeptr=edgeptr->next_to; continue; } // skip non-initial nodes
            b += exp(edgeptr->from_node->pi)*edgeptr->t;
            edgeptr = edgeptr->next_to;
        }
        hitprob_f << setw(5) << it_nodevec->node_id << setw(18) << b << endl;
        it_nodevec++;
    }
    cout << "kps> finished writing absorption probabilities to files" << endl;
}

/* write the elements of the fundamental matrix for an absorbing (i.e. reducible) Markov chain */
void KPS::calc_fundamentalred(const Network &ktn) {
    cout << "kps> calculating expected numbers of node visits from the graph transformation" << endl;
    write_renormalised_probs("transient_visits.dat"); // compute expected number of visits of node j given that trajectories are initialised at node i
    if (gth) return; // the stationary distribution is not known at this point
    // ...
    cout << "kps> finished writing expected numbers of node visits to files" << endl;
}

/* write the elements of the vector of MFPTs to the absorbing state. NB At this point, pi values of the ktn_kps
   object should, for initial nodes, have been overwritten to the initial probability distribution values */
void KPS::calc_mfpt() {
    cout << "kps> writing MFPTs for transitions from all non-absorbing nodes to file" << endl;
    ofstream mfpt_f; mfpt_f.open("mfpts.dat"); mfpt_f.setf(ios::scientific,ios::floatfield);
    mfpt_f.precision(10);
    vector<Node>::iterator it_nodevec = ktn_kps->nodes.begin();
    long double mfpt_ab = 0.L; // calculate total A<-B MFPT given the initial probability distribution within the set B
    while (it_nodevec!=ktn_kps->nodes.end()) {
        if (it_nodevec->aorb==-1) { it_nodevec++; continue; // skip absorbing nodes, for which the MFPT is not defined
        } else if (it_nodevec->aorb==1) { mfpt_ab += exp(it_nodevec->pi)*mfpt_vals[it_nodevec->node_pos]; }
        mfpt_f << setw(5) << it_nodevec->node_id << setw(18) << mfpt_vals[it_nodevec->node_pos] << endl;
        it_nodevec++;
    }
    if (gth) return; // the stationary distribution is not known at this point
    cout << "kps> the A<-B MFPT is:" << string(10,' ') << setw(18) << scientific << setprecision(10) << mfpt_ab << endl;
    cout << "kps> finished writing MFPTs to file" << endl;
}

/* compute and write the stationary probabilities determined by the GTH algorithm */
void KPS::calc_gth() {
    cout << "kps> writing stationary probabilities determined by the GTH algorithm to file" << endl;
    cout << "mu is: " << mu << endl;
    vector<long double> gth_pi_vals(ktn_kps->n_nodes);
    for (vector<Node>::iterator it_nodevec=ktn_kps->nodes.begin();it_nodevec!=ktn_kps->nodes.end();++it_nodevec) {
        it_nodevec->pi *= 1.L/mu; gth_pi_vals[it_nodevec->node_pos] = it_nodevec->pi; }
    Wrapper_Method::write_vec<long double>(gth_pi_vals,"stat_prob_gth.dat");
    cout << "kps> finished writing stationary distribution to file" << endl;
}


/* rewrite the stationary probabilities of the ktn_kps network to reflect the initial distribution */
void KPS::rewrite_stat_probs(const Network &ktn) {
    set<const Node*>::iterator it_set = ktn.nodesB.begin();
    if (ktn.nodesB.size()==1) { // there is only one node in the initial set
        ktn_kps->nodes[nodemap[(*it_set)->node_id]-1].pi=0.L;
    } else if (ktn.initcond) { // specified initial probability distribution from file
        int i=0;
        while (it_set!=ktn.nodesB.end()) {
            ktn_kps->nodes[nodemap[(*it_set)->node_id]-1].pi = ktn.init_probs[i];
            i++; it_set++;
        }
    } else { // local equilibrium distribution within initial set
        long double pi_B = -numeric_limits<long double>::infinity(); // (log) occupation probability of all nodes in initial set B
        while (it_set!=ktn.nodesB.end()) {
            pi_B = log(exp(pi_B)+exp((*it_set)->pi));
            it_set++;
        }
        it_set = ktn.nodesB.begin();
        while (it_set!=ktn.nodesB.end()) {
            ktn_kps->nodes[nodemap[(*it_set)->node_id]-1].pi -= pi_B;
            it_set++;
        }
    }
}

/* write the elements of the graph-transformed network to a file */
void KPS::write_renormalised_probs(string fname) {
    ofstream elems_f; elems_f.open(fname);
    elems_f.setf(ios::right,ios::adjustfield); elems_f.setf(ios::scientific,ios::floatfield);
    elems_f.precision(10);
    for (vector<Node>::iterator it_nodevec=ktn_kps->nodes.begin();it_nodevec!=ktn_kps->nodes.end();++it_nodevec) {
        if (fundamentalred && !it_nodevec->flag) continue;
        if (!it_nodevec->eliminated && it_nodevec->aorb!=-1) { // print self-loop of non-absorbing node if node is non-eliminated
            elems_f << setw(5) << it_nodevec->node_id << setw(5) << it_nodevec->node_id << setw(18) << it_nodevec->t << endl;
        }
        Edge *edgeptr = it_nodevec->top_from;
        while (edgeptr!=nullptr) {
            if (edgeptr->deadts || edgeptr->to_node->eliminated) { edgeptr=edgeptr->next_from; continue; }
            elems_f << setw(5) << edgeptr->from_node->node_id << setw(5) << edgeptr->to_node->node_id \
                    << setw(18) << edgeptr->t << endl;
            edgeptr = edgeptr->next_from;
        }
    }
}

#endif
