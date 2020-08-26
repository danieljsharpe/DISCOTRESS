/*
File containing functions to construct and handle a Network (transition network) object

This file is a part of DISCOTRESS, a software package to simulate the dynamics on arbitrary continuous- and discrete-time Markov chains (CTMCs and DTMCs).
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

#include "ktn.h"
#include <vector>
#include <queue>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

Network::Network(int nnodes, int nedges) {
    nodes.resize(nnodes); n_nodes=nnodes;
    edges.resize(2*nedges); n_edges=nedges;
}

Network::~Network() {}

Network::Network(const Network &ktn) {
    n_nodes=ktn.n_nodes; n_edges=ktn.n_edges;
    nodes.resize(n_nodes); edges.resize(n_edges);
    for (int i=0;i<n_nodes;i++) nodes[i] = ktn.nodes[i];
    for (int i=0;i<n_edges;i++) edges[i] = ktn.edges[i];
    n_dead=ktn.n_dead; ncomms=ktn.ncomms;
    branchprobs=ktn.branchprobs; tau=ktn.tau;
    // now sort out pointers
    for (int i=0;i<n_edges;i++) {
        edges[i].from_node = &nodes[ktn.edges[i].from_node->node_pos];
        edges[i].to_node = &nodes[ktn.edges[i].to_node->node_pos];
        add_to_edge(ktn.edges[i].to_node->node_pos,i);
        add_from_edge(ktn.edges[i].from_node->node_pos,i);
        edges[i].rev_edge = &edges[ktn.edges[i].rev_edge->edge_id];
    }
}

Node::Node() {}

Node::~Node() {}

/* copy constructor for Node copies properties but not pointers to Edge types */
Node::Node(const Node &node) {
    node_id=node.node_id; comm_id=node.comm_id; aorb=node.aorb;
    udeg=node.udeg; eliminated=node.eliminated;
    k_esc=node.k_esc; t=node.t; pi=node.pi;
}

/* print mean waiting times for nodes to file */
void Network::dumpwaittimes() {
    ofstream tau_f; // file containing mean waiting times for nodes
    tau_f.open("meanwaitingtimes.dat");
    tau_f.setf(ios::scientific,ios::floatfield); tau_f.precision(20);
    for (const Node &node: nodes) {
        tau_f << 1./exp(node.k_esc) << endl; }
}

/* update the Nodes and Edges of the Network data structure to contain transition probabibilities
   calculated from the linearised transition probabibility matrix */
void Network::get_tmtx_lin(long double tau) {
    cout << "ktn> calculating linearised transition probability matrix at lag time: " << tau << endl;
    for (auto &node: nodes) {
        if (tau>1./exp(node.k_esc)) throw Ktn_exception(); // value of tau does not give stochastic matrix
        node.t = 1.-(exp(node.k_esc)*tau);
        node.k_esc = log(tau); // override node escape rates
    }
    for (auto &edge: edges) {
        if (edge.deadts) continue;
        edge.t = exp(edge.k)*tau; }
    this->tau=tau;
}

/* the transition probabilities are calculated as the branching probabilities */
void Network::get_tmtx_branch() {
    cout << "ktn> calculating branching probability matrix" << endl;
    branchprobs=true;
    for (auto & edge: edges) {
        if (edge.deadts) continue;
        edge.t = exp(edge.k-edge.from_node->k_esc); }
    for (auto &node: nodes) {
        node.t = 0.; // branching probability matrix contains no self-loops
        long double cum_t=0.; // accumulated branching probability
        Edge *edgeptr = node.top_from;
        while (edgeptr!=nullptr) {
            if (!edgeptr->deadts) { cum_t += edgeptr->t; }
            edgeptr=edgeptr->next_from;
        }
        if (abs(cum_t-1.)>1.E-14) throw Ktn_exception(); // transition probabilities do not give stochastic matrix
    }
}

/* calculate the branching probabilities, resort edges into decreasing order of branching probabilities, and set the transition
   probabilities as the accumulated values. This is for optimisation of the rejection-free kMC algorithm (BKL) */
void Network::get_cum_branchprobs() {
    this->get_tmtx_branch();
    cout << "ktn> calculating accumulative branching probabilities" << endl;
    accumprobs=true;
    for (auto &node: nodes) {
        Edge *edgeptr = node.top_from;
        auto cmp = [](Edge *l, Edge *r) { return l->t > r->t; };
        priority_queue<Edge*,vector<Edge*>,decltype(cmp)> edge_pq(cmp); // priority queue of edges (based on branching probability)
        while (edgeptr!=nullptr) {
            Edge *next_edge = edgeptr->next_from;
            edge_pq.push(&(*edgeptr));
            del_from_edge(node.node_id-1);
            edgeptr->next_from = nullptr;
            edgeptr = next_edge;
        }
        while (!edge_pq.empty()) { // edges will be added so that the final linked list is in order of decreasing t
            this->add_from_edge(node.node_id-1,edge_pq.top()->edge_id);
            edge_pq.pop();
        }
        edgeptr = node.top_from;
        long double cum_t=0.;
        while (edgeptr!=nullptr) { // calculate accumulated branching probabilities
            if (edgeptr->deadts) { edgeptr=edgeptr->next_from; continue; }
            long double prev_cum_t = cum_t;
            cum_t += edgeptr->t;
            edgeptr->t += prev_cum_t;
            edgeptr=edgeptr->next_from;
        }
        if (abs(cum_t-1.)>1.E-16) throw Ktn_exception();
    }
}

// delete node i 
void Network::del_node(int i) {
    if (nodes[i].eliminated) throw Ktn_exception();
    Edge *edgeptr;
    edgeptr = nodes[i].top_to;
    while (edgeptr!=nullptr) {
        del_to_edge(i);
        edgeptr = edgeptr->next_to;
    }
    edgeptr = nodes[i].top_from;
    while (edgeptr!=nullptr) {
        del_from_edge(i);
        edgeptr = edgeptr->next_from;
    }
    nodes[i].eliminated = true;
    tot_nodes--;
}

// edge j goes TO node i
void Network::add_to_edge(int i, int j) {
    if (nodes[i].top_to != nullptr) {
        edges[j].next_to = nodes[i].top_to;
        nodes[i].top_to = &edges[j]; }
    else {
        nodes[i].top_to = &edges[j];
        nodes[i].top_to->next_to = nullptr; }
    tot_edges++;
}

// edge j goes FROM node i
void Network::add_from_edge(int i, int j) {
    if (nodes[i].top_from != nullptr) {
        edges[j].next_from = nodes[i].top_from;
        nodes[i].top_from = &edges[j]; }
    else {
        nodes[i].top_from = &edges[j];
        nodes[i].top_from->next_from = nullptr; }
    tot_edges++;
    nodes[i].udeg++;
}

// delete the top TO edge for node i
void Network::del_to_edge(int i) {
    if (nodes[i].top_to != nullptr) {
        if (nodes[i].top_to->next_to != nullptr) {
            nodes[i].top_to = nodes[i].top_to->next_to;
        } else {
            nodes[i].top_to = nullptr;
        }
        tot_edges--;
    } else {
        throw Ktn_exception();
    }
}

// delete the top FROM edge for node i
void Network::del_from_edge(int i) {
    if (nodes[i].top_from != nullptr) {
        if (nodes[i].top_from->next_from != nullptr) {
            nodes[i].top_from = nodes[i].top_from->next_from;
        } else {
            nodes[i].top_from = nullptr;
        }
        tot_edges--;
    } else {
        throw Ktn_exception();
    }
    nodes[i].udeg--;
}

// delete TO edge with edge_id j for node i
void Network::del_spec_to_edge(int i, int j) {
    Edge *edgeptr; Edge *edgeptr_prev = nullptr;
    bool edge_exists = false;
    if (nodes[i].top_to==nullptr) throw Ktn_exception();
    edgeptr = nodes[i].top_to;
    while (edgeptr!=nullptr) {
        if (edgeptr->edge_id==j) {
            if (edgeptr_prev != nullptr) {
                edgeptr_prev->next_to = edgeptr->next_to;
            } else if (edgeptr->next_to != nullptr) {
                nodes[i].top_to = edgeptr->next_to;
            } else if (edgeptr->next_to == nullptr) {
                nodes[i].top_to = nullptr;
            }
            edge_exists = true;
            break;
        }
        edgeptr_prev = edgeptr;
        edgeptr = edgeptr->next_to;
    }
    if (!edge_exists) throw Ktn_exception();
}

// delete FROM edge with edge_id j for node i
void Network::del_spec_from_edge(int i, int j) {
    Edge *edgeptr; Edge *edgeptr_prev = nullptr;
    bool edge_exists = false;
    if (nodes[i].top_from==nullptr) throw Ktn_exception();
    edgeptr = nodes[i].top_from;
    while (edgeptr!=nullptr) {
        if (edgeptr->edge_id==j) {
            if (edgeptr_prev != nullptr) {
                edgeptr_prev->next_from = edgeptr->next_from;
            } else if (edgeptr->next_from != nullptr) {
                nodes[i].top_from = edgeptr->next_from;
            } else if (edgeptr->next_from ==nullptr) {
                nodes[i].top_from = nullptr;
            }
            edge_exists = true;
            break;
        }
        edgeptr_prev = edgeptr;
        edgeptr = edgeptr->next_from;
    }
    if (!edge_exists) throw Ktn_exception();
}

// update edge with edge_id j so that it now points TO node i
void Network::update_to_edge(int i, int j) {
    Edge *edgeptr;
    edgeptr = &edges[j];
    int old_to = edgeptr->to_node->node_id;
    edgeptr->to_node = &nodes[i];
    del_spec_to_edge(old_to-1,edgeptr->edge_id);
    add_to_edge(i,j);
}

// update edge with edge_id j so that it now points FROM i
void Network::update_from_edge(int i, int j) {
    Edge *edgeptr;
    edgeptr = &edges[j];
    int old_from = edgeptr->from_node->node_id;
    edgeptr->from_node = &nodes[i];
    del_spec_from_edge(old_from-1,edgeptr->edge_id);
    add_from_edge(i,j);
}

/* calculate the escape rate for node */
void Network::calc_k_esc(Node &node) {
    Edge *edgeptr;
    node.k_esc = -numeric_limits<long double>::infinity();
    edgeptr = node.top_from;
    while (edgeptr!=nullptr) {
        if (!edgeptr->deadts) node.k_esc = log(exp(node.k_esc)+exp(edgeptr->k));
        edgeptr = edgeptr->next_from;
    }
}

/* calculate the self-loop transition probability for node */
void Network::calc_t_selfloop(Node &node) {
    Edge *edgeptr;
    node.t=1.;
    edgeptr = node.top_from;
    while (edgeptr!=nullptr) {
        if (!edgeptr->deadts) node.t -= edgeptr->t;
        edgeptr = edgeptr->next_from;
    }
}

/* calculate the net flux along an edge and its reverse edge */
long double Network::calc_net_flux(Edge &edge) {
    if (edge.rev_edge==nullptr) throw Network::Ktn_exception();
    if (!((edge.to_node->node_id==edge.rev_edge->from_node->node_id) || \
          (edge.from_node->node_id==edge.rev_edge->to_node->node_id))) throw Network::Ktn_exception();
    return exp(edge.k+edge.from_node->pi)-exp(edge.rev_edge->k+edge.to_node->pi);
}

/* set the vector specifying the initial probabilities of nodes in set B (if different from relative stationary probabilities) */
void Network::set_initcond(const vector<double> &init_probs) {

    if (init_probs.size()!=nodesB.size()) throw Network::Ktn_exception();
    double p_tot=0.;
    for (const auto &p: init_probs) p_tot+=p;
    if (abs(p_tot-1.)>1.E-10) throw Network::Ktn_exception();
    initcond=true;
    this->init_probs=init_probs;
}

/* update the Network object pointed to by the ktn argument to include an additional edge (with index k in the edges array)
   connecting from_node and to_node */
void Network::add_edge_network(Network *ktn, Node &from_node, Node &to_node, int k) {
    (ktn->edges[k]).from_node = &from_node;
    (ktn->edges[k]).to_node = &to_node;
    ktn->add_from_edge(from_node.node_id-1,k);
    ktn->add_to_edge(to_node.node_id-1,k);
}

/* set up the kinetic transition network */
void Network::setup_network(Network& ktn, const vector<pair<int,int>> &ts_conns, \
        const vector<long double> &ts_wts, const vector<double> &stat_probs, const vector<int> &nodesinA, \
        const vector<int> &nodesinB, bool transnprobs, long double tau, int ncomms, const vector<int> &comms,
        const vector<int> &bins) {

    cout << "ktn> constructing nodes and edges of transition network from vectors" << endl;
    if (!((ts_conns.size()==ktn.n_edges) || (ts_wts.size()==2*ktn.n_edges) || \
         (stat_probs.size()==ktn.n_nodes) || \
         (!comms.empty() && comms.size()==ktn.n_nodes))) throw Network::Ktn_exception();
    if (transnprobs) {
        cout << "ktn> interpreting edge weights as transition probabilities at a lag time: " << tau << endl;
        ktn.tau=tau; }
    ktn.ncomms=ncomms;
    if (!comms.empty()) ktn.comm_sizes.resize(ncomms);
    double tot_pi = -numeric_limits<double>::infinity();
    for (int i=0;i<ktn.n_nodes;i++) {
        ktn.nodes[i].node_id = i+1; ktn.nodes[i].node_pos = i;
        if (!comms.empty()) {
            ktn.nodes[i].comm_id = comms[i];
            ktn.nodes[i].bin_id = bins[i];
            ktn.comm_sizes[comms[i]]++;
            if (comms[i]>=ncomms) throw Ktn_exception();
            if (bins[i]+1>ktn.nbins) ktn.nbins=bins[i]+1;
        }
        ktn.nodes[i].pi = stat_probs[i];
        tot_pi = log(exp(tot_pi) + exp(stat_probs[i]));
    }
    tot_pi = exp(tot_pi);
    if (abs(tot_pi-1.)>1.E-10) {
        cout << "ktn> Error: total equilibrium probabilities of nodes is: " << tot_pi << " =/= 1." << endl;
        throw Network::Ktn_exception(); }
    for (int i=0;i<ktn.n_edges;i++) {
        ktn.edges[2*i].edge_id = 2*i;
        ktn.edges[(2*i)+1].edge_id = (2*i)+1;
        if (ts_conns[i].first == ts_conns[i].second) { // "dead" transition state (dangling node)
            cout << "ktn> warning: transition state " << i << " is dead" << endl;
            ktn.edges[2*i].deadts = true;
            ktn.edges[(2*i)+1].deadts = true;
            ktn.n_dead++;
            continue;
        } else {
            ktn.edges[2*i].deadts = false;
            ktn.edges[(2*i)+1].deadts = false;
        }
        if (!transnprobs) { // ts_wts are transition rates
            ktn.edges[2*i].k = ts_wts[2*i];
            ktn.edges[(2*i)+1].k = ts_wts[(2*i)+1];
        } else { // ts_wts are transition probabilities
            ktn.edges[2*i].k = 0.; ktn.edges[(2*i)+1].k = 0.; // dummy values for transition rates
            ktn.edges[2*i].t = ts_wts[2*i];
            ktn.edges[(2*i)+1].t = ts_wts[(2*i)+1];
        }
        ktn.edges[2*i].from_node = &ktn.nodes[ts_conns[i].first-1];
        ktn.edges[2*i].to_node = &ktn.nodes[ts_conns[i].second-1];
        ktn.edges[(2*i)+1].from_node = &ktn.nodes[ts_conns[i].second-1];
        ktn.edges[(2*i)+1].to_node = &ktn.nodes[ts_conns[i].first-1];

        ktn.add_to_edge(ts_conns[i].second-1,2*i);
        ktn.add_from_edge(ts_conns[i].first-1,2*i);
        ktn.add_to_edge(ts_conns[i].first-1,(2*i)+1);
        ktn.add_from_edge(ts_conns[i].second-1,(2*i)+1);

        ktn.edges[2*i].rev_edge = &ktn.edges[(2*i)+1];
        ktn.edges[(2*i)+1].rev_edge = &ktn.edges[2*i];
    }
    for (int i=0;i<ktn.n_nodes;i++) {
        if (!transnprobs) { // continuous-time, in general, nodes have different waiting times (escape rates)
            calc_k_esc(ktn.nodes[i]);
        } else {
            calc_t_selfloop(ktn.nodes[i]);
            ktn.nodes[i].k_esc = log(1./tau); // continuous- or discrete-time, all nodes have same escape rate
        }
    }
    for (int i=0;i<nodesinA.size();i++) {
        if (nodesinA[i]>ktn.n_nodes) throw Ktn_exception();
        ktn.nodes[nodesinA[i]-1].aorb = -1;
        ktn.nodesA.insert(&ktn.nodes[nodesinA[i]-1]);
    }
    for (int i=0;i<nodesinB.size();i++) {
        if (nodesinB[i]>ktn.n_nodes) throw Ktn_exception();
        ktn.nodes[nodesinB[i]-1].aorb = 1;
        ktn.nodesB.insert(&ktn.nodes[nodesinB[i]-1]);
    }
    /* the community IDs of A (and of B) nodes should be consistent, but note that it is *not* required that the
       number of initial (B) nodes is equal to the number of nodes in that community - i.e. (I \cup B) is allowed to
       be a single community */
    if (!comms.empty() && !nodesinA.empty()) { // community IDs of A (and of B) nodes should be consistent
        if (ktn.nodesA.size()!=ktn.comm_sizes[(*ktn.nodesA.begin())->comm_id]) throw Ktn_exception();
        int commA=(*ktn.nodesA.begin())->comm_id;
        for (set<const Node*>::iterator it_set=ktn.nodesA.begin();it_set!=ktn.nodesA.end();++it_set) {
            if ((*it_set)->comm_id!=commA) throw Ktn_exception(); }
        int commB=(*ktn.nodesB.begin())->comm_id;
        for (set<const Node*>::iterator it_set=ktn.nodesB.begin();it_set!=ktn.nodesB.end();++it_set) {
            if ((*it_set)->comm_id!=commB) throw Ktn_exception(); }
    }
    cout << "ktn> finished setting up transition network data structure" << endl;
}
