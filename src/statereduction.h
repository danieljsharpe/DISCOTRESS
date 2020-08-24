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
        q_ab_vals[node.node_id-1] = q_ab/(q_ab+q_ba);
        q_ba_vals[node.node_id-1] = q_ba/(q_ab+q_ba);
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

/* calculate the committor probability for an initial node at the boundary of the initial state, which is =/= 0 */
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

/* compute and write absorption probabilities */
void KPS::calc_absprobs(const Network &ktn) {}

/* write the elements of the fundamental matrix for an absorbing (i.e. reducible) Markov chain */
void KPS::write_fundamentalred(const Network &ktn) {}
