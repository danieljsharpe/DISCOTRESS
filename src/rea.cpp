/*
File containing functions relating to the recursive enumeration algorithm (REA) for the k shortest paths problem.

Here the REA is used to determine the k first passage A<-B paths of the Markov chain with the highest probability. See:
V. M. Jimenez and A. Marzal: "Computing the k shortest paths: a new algorithm and experimental comparison," in Algorithm Engineering:
3rd International Workshop, WAE '99, London, UK, ed. J. S. Vitter and C. D. Zaroliagis (Springer Berlin, Heidelberg, 1999) pp. 15-29.

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

#include "kmc_methods.h"
#include <iostream>

using namespace std;

REA::REA(const Network &ktn, const Wrapper_args &wrapper_args) : Wrapper_Method(wrapper_args) {
    cout << "rea> using the REA to determine the " << wrapper_args.nabpaths << " highest-probability paths" << endl;
    shortest_paths.resize(ktn.n_nodes); candidate_paths.resize(ktn.n_nodes);
    for (int i=0;i<ktn.n_nodes;i++) { // allocate arrays for candidate and assigned shortest paths
        shortest_paths[i].resize(wrapper_args.nabpaths);
        for (int k=1;k<wrapper_args.nabpaths+1;k++) {
            shortest_paths[i][k-1] = { walker_id:0,path_no:k,k:0,t:0.L, \
                                       p:-numeric_limits<long double>::infinity(),s:0.L,prev_node:nullptr,curr_node:&ktn.nodes[i] };
            shortest_paths[i][k-1].visited.resize(wrapper_args.nbins);
            fill(shortest_paths[i][k-1].visited.begin(),shortest_paths[i][k-1].visited.end(),false);
        }
        candidate_paths[i].resize(ktn.nodes[i].udeg); // edges are bidirectional, so in- and out-degrees of nodes are the same
        for (int j=0;j<ktn.nodes[i].udeg;j++) candidate_paths[i][j]=nullptr;
    }
}

REA::~REA() {}

void REA::run_enhanced_kmc(const Network &ktn, Traj_Method *traj_method_obj) {

    source = *ktn.nodesB.begin(); // NB there is only a single source node
    sink = *ktn.nodesA.begin(); // NB there is only a single sink node
    dijkstra(ktn);
    for (int k=2;k<nabpaths+1;k++) next_path(ktn,ktn.nodes[0],k); // main loop of REA
    if (ktn.ncomms>0) {
//        for (int i=0;i<nabpaths;i++) print_visits();
    }
}

/* compute the first shortest path from the source node to all other nodes using Dijkstra's algorithm */
void REA::dijkstra(const Network& ktn) {
    cout << "Dijkstra's algorithm" << endl;
}

/* for shortest paths k>=2, and given that the 1,...,(k-1)-th shortest paths to node v have been computed, find
   the k-th shortest path to node v */
void REA::next_path(const Network &ktn, const Node &vnode, int k) {
    cout << "next_path()" << endl;
}

/* add the union of: { k-th shortest path to node u } \cup node v
   as a possible candidate for the next shortest path to node v  */
void REA::add_candidate(int u, int v, int k) {
    cout << "add_candidate()" << endl;
}

/* write the nodes that are visited along a shortest path to the sink (target) node to an output file */
void REA::print_visits(const Walker&) {
    cout << "print_visits()" << endl;
}
