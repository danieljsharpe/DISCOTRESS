#include "ktn.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <limits>

using namespace std;

Network::Network(int nnodes, int nedges) {
    nodes.resize(nnodes);
    edges.resize(2*nedges);
    n_nodes=nnodes; n_edges=nedges;
}

Network::~Network() {}

Network::Network(const Network &ktn) {}

/* update the Nodes and Edges of the Network data structure to contain transition probabibilities
   calculated from the linearised transition probabibility matrix */
void Network::get_tmtx_lin(double tau) {
    for (auto &node: nodes) {
        if (tau>1./exp(node.k_esc)) throw Ktn_exception(); // value of tau does not give stochastic matrix
        node.t = 1.-(exp(node.k_esc)*tau); }
    for (auto &edge: edges) {
        edge.t = exp(edge.k)*tau; }
}

// delete node i 
void Network::del_node(int i) {
    if (nodes[i].deleted) throw Ktn_exception();
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
    nodes[i].deleted = true;
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
}

// delete TO edge with ts_id j for node i
void Network::del_spec_to_edge(int i, int j) {
    Edge *edgeptr; Edge *edgeptr_prev = nullptr;
    bool ts_exists = false;
    if (nodes[i].top_to==nullptr) throw Ktn_exception();
    edgeptr = nodes[i].top_to;
    while (edgeptr!=nullptr) {
        if (edgeptr->ts_id==j) {
            if (edgeptr_prev != nullptr) {
                edgeptr_prev->next_to = edgeptr->next_to;
            } else if (edgeptr->next_to != nullptr) {
                nodes[i].top_to = edgeptr->next_to;
            } else if (edgeptr->next_to == nullptr) {
                nodes[i].top_to = nullptr;
            }
            ts_exists = true;
            break;
        }
        edgeptr_prev = edgeptr;
        edgeptr = edgeptr->next_to;
    }
    if (!ts_exists) throw Ktn_exception();
}

// delete FROM edge with ts_id j for node i
void Network::del_spec_from_edge(int i, int j) {
    Edge *edgeptr; Edge *edgeptr_prev = nullptr;
    bool ts_exists = false;
    if (nodes[i].top_from==nullptr) throw Ktn_exception();
    edgeptr = nodes[i].top_from;
    while (edgeptr!=nullptr) {
        if (edgeptr->ts_id==j) {
            if (edgeptr_prev != nullptr) {
                edgeptr_prev->next_from = edgeptr->next_from;
            } else if (edgeptr->next_from != nullptr) {
                nodes[i].top_from = edgeptr->next_from;
            } else if (edgeptr->next_from ==nullptr) {
                nodes[i].top_from = nullptr;
            }
            ts_exists = true;
            break;
        }
        edgeptr_prev = edgeptr;
        edgeptr = edgeptr->next_from;
    }
    if (!ts_exists) throw Ktn_exception();
}

// update edge so that it now points TO i
void Network::update_to_edge(int i, int j) {
    Edge *edgeptr;
    edgeptr = &edges[j];
    int old_to = edgeptr->to_node->node_id;
    edgeptr->to_node = &nodes[i];
    del_spec_to_edge(old_to-1,edgeptr->ts_id);
    add_to_edge(i,j);
}

// update edge so that it now points FROM i
void Network::update_from_edge(int i, int j) {
    Edge *edgeptr;
    edgeptr = &edges[j];
    int old_from = edgeptr->from_node->node_id;
    edgeptr->from_node = &nodes[i];
    del_spec_from_edge(old_from-1,edgeptr->ts_id);
    add_from_edge(i,j);
}

/* calculate the escape rate for node i */
void Network::calc_k_esc(Node &node) {
    Edge *edgeptr;
    node.k_esc = -numeric_limits<double>::infinity();
    edgeptr = node.top_from;
    while (edgeptr!=nullptr) {
        if (!edgeptr->deadts) node.k_esc = log(exp(node.k_esc)+exp(edgeptr->k));
        edgeptr = edgeptr->next_from;
    }
}

/* calculate the net flux along an edge and its reverse edge */
void Network::calc_net_flux(Edge &edge) {
    if (edge.rev_edge==nullptr) throw Network::Ktn_exception();
    if (!((edge.to_node->node_id==edge.rev_edge->from_node->node_id) || \
          (edge.from_node->node_id==edge.rev_edge->to_node->node_id))) throw Network::Ktn_exception();
    edge.j = exp(edge.k+edge.from_node->pi)-exp(edge.rev_edge->k+edge.to_node->pi);
    edge.rev_edge->j = -edge.j;
}

/* set up the kinetic transition network */
void Network::setup_network(Network& ktn, const vector<pair<int,int>> &ts_conns, \
        const vector<double> &ts_wts, const vector<double> &stat_probs, const vector<int> &nodesinA, \
        const vector<int> &nodesinB, const vector<int> &comms) {

    if (!((ts_conns.size()==ktn.n_edges) || (ts_wts.size()==2*ktn.n_edges) || \
         (stat_probs.size()==ktn.n_nodes))) throw Network::Ktn_exception();
    double tot_pi = -numeric_limits<double>::infinity();
    for (int i=0;i<ktn.n_nodes;i++) {
        ktn.nodes[i].node_id = i+1;
        ktn.nodes[i].comm_id = comms[i];
        ktn.nodes[i].pi = stat_probs[i];
        tot_pi = log(exp(tot_pi) + exp(stat_probs[i]));
    }
    tot_pi = exp(tot_pi);
    if (abs(tot_pi-1.)>1.E-10) {
        cout << "Error: total equilibrium probabilities of minima is: " << tot_pi << " =/= 1." << endl;
        throw Network::Ktn_exception(); }
    for (int i=0;i<ktn.n_edges;i++) {
        ktn.edges[2*i].ts_id = i+1;
        ktn.edges[(2*i)+1].ts_id = i+1;
        ktn.edges[2*i].edge_pos = 2*i;
        ktn.edges[(2*i)+1].edge_pos = (2*i)+1;
        if (ts_conns[i].first == ts_conns[i].second) { // "dead" transition state (dangling node)
            ktn.edges[2*i].deadts = true;
            ktn.edges[(2*i)+1].deadts = true;
            ktn.n_dead++;
            continue;
        } else {
            ktn.edges[2*i].deadts = false;
            ktn.edges[(2*i)+1].deadts = false;
        }
        ktn.edges[2*i].k = ts_wts[2*i];
        ktn.edges[(2*i)+1].k = ts_wts[(2*i)+1];
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

        calc_net_flux(ktn.edges[2*i]);
    }
    for (int i=0;i<ktn.n_nodes;i++) {
        calc_k_esc(ktn.nodes[i]);
    }
    for (int i=0;i<nodesinA.size();i++) {
        ktn.nodes[nodesinA[i]-1].aorb = -1;
        ktn.nodesA.insert(ktn.nodes[nodesinA[i]-1]);
    }
    for (int i=0;i<nodesinB.size();i++) {
        ktn.nodes[nodesinB[i]-1].aorb = 1;
        ktn.nodesB.insert(ktn.nodes[nodesinB[i]-1]);
    }

    cout << "ktn> finished reading in kinetic transition network" << endl;
}
