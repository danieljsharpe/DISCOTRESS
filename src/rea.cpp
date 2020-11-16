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
#include <cmath>
#include <iostream>

using namespace std;

REA::REA(const Network &ktn, bool discretetime, const Wrapper_args &wrapper_args) : Wrapper_Method(wrapper_args) {
    cout << "rea> using the REA to determine the " << wrapper_args.nabpaths << " highest-probability paths" << endl;
    this->discretetime=discretetime;
    source_node = *ktn.nodesB.begin(); // NB there is only a single source node
    sink_node = *ktn.nodesA.begin(); // NB there is only a single sink node
    shortest_paths.resize(ktn.n_nodes); candidate_paths.resize(ktn.n_nodes);
    for (int i=0;i<ktn.n_nodes;i++) { // allocate arrays for candidate and assigned shortest paths
        shortest_paths[i].resize(wrapper_args.nabpaths);
        for (int k=1;k<wrapper_args.nabpaths+1;k++) { // path cost is initially infinite and predecessor node not set for all paths
            shortest_paths[i][k-1] = { walker_id:0,path_no:k,k:0,t:0.L, \
                                       p:numeric_limits<long double>::infinity(),s:0.L,prev_node:nullptr,curr_node:&ktn.nodes[i] };
            shortest_paths[i][k-1].visited.resize(wrapper_args.nbins);
            fill(shortest_paths[i][k-1].visited.begin(),shortest_paths[i][k-1].visited.end(),false);
        }
        candidate_paths[i].resize(ktn.nodes[i].udeg); // edges are bidirectional, so in- and out-degrees of nodes are the same
        for (int j=0;j<ktn.nodes[i].udeg;j++) {
            candidate_paths[i][j].first=nullptr; candidate_paths[i][j].second=nullptr;
        }
    }
    if (!ktn.nbins>0) return;
    for (int k=1;k<wrapper_args.nabpaths+1;k++) { // sink (target) node is visited along all shortest paths to sink node
        shortest_paths[sink_node->node_id-1][k-1].visited[sink_node->bin_id]=true; } // this is needed for writing visits.x.dat files
}

REA::~REA() {}

void REA::run_enhanced_kmc(const Network &ktn, Traj_Method *traj_method_obj) {

    dijkstra(ktn); // find shortest (i.e. highest-probability) path
    print_dijkstra(); // write shortest path to output file
    for (int k=2;k<nabpaths+1;k++) next_path(ktn,ktn.nodes[0],k); // main loop of REA
    if (ktn.nbins>0) {
        for (int k=1;k<nabpaths+1;k++) print_visits(&shortest_paths[sink_node->node_id-1][k-1]); }
}

/* compute the first shortest path from the source node to all other nodes using Dijkstra's algorithm */
void REA::dijkstra(const Network& ktn) {
    cout << "Dijkstra's algorithm" << endl;
    int m, n;
    vector<bool> insptree(ktn.n_nodes,false); // vector for bookkeeping which nodes have been incorporated into the shortest path tree
    // initialisation
    const Node *curr_node=source_node;
    shortest_paths[curr_node->node_id-1][0].p=0.L;
    // main loop for Dijkstra's algorithm
    for (int i=0;i<ktn.n_nodes;i++) {
        n=curr_node->node_id-1;
        insptree[n]=true;
        const Edge *edgeptr=curr_node->top_from;
        while (edgeptr!=nullptr) { // loop over outgoing edges
            m=edgeptr->to_node->node_id-1;
            if (shortest_paths[n][0].p - 1.L*log(edgeptr->t) < shortest_paths[m][0].p) {
                // update path values
                shortest_paths[m][0].p = shortest_paths[n][0].p - 1.L*log(edgeptr->t);
                shortest_paths[m][0].k = shortest_paths[n][0].k + 1;
                shortest_paths[m][0].t = shortest_paths[n][0].t + edgeptr->from_node->t_esc;
                if (!discretetime) { shortest_paths[m][0].s = shortest_paths[n][0].s + (edgeptr->rev_edge->k-edgeptr->k);
                } else { shortest_paths[m][0].s = shortest_paths[n][0].s + log (edgeptr->rev_edge->t/edgeptr->t); }
                shortest_paths[m][0].prev_node = curr_node; // set previous node in shortest path tree
                for (int j=0;j<ktn.nbins;j++) { // update boolean values representing bin visits
                    if (j==curr_node->node_id-1) { shortest_paths[m][0].visited[j]=true;
                    } else if (j!=sink_node->node_id-1) { shortest_paths[m][0].visited[j] = shortest_paths[n][0].visited[j]; }
                }
            }
            edgeptr=edgeptr->next_from;
        }
        // find node with current lowest shortest path cost
        long double mincost=numeric_limits<long double>::infinity();
        for (int j=0;j<ktn.n_nodes;j++) {
            if (!insptree[j] && (shortest_paths[j][0].p < mincost)) {
                mincost = shortest_paths[j][0].p;
                curr_node = &ktn.nodes[j];
            }
        }
    }
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

/* print the first shortest path, determined by Dijkstra's algorithm, by tracing the shortest path tree */
void REA::print_dijkstra() {
    cout << "print_dijkstra()" << endl;
    ofstream spath_f;
    string spath_fname="shortest_path.dat";
    spath_f.open(spath_fname,ios_base::trunc);
    spath_f.setf(ios::right,ios::adjustfield); spath_f.setf(ios::scientific,ios::floatfield);
    spath_f.precision(10);
    const Node *curr_node = sink_node;
    const Walker* curr_path;
    while (curr_node!=nullptr) {
        curr_path = &shortest_paths[curr_node->node_id-1][0];
        spath_f << setw(7) << curr_node->node_id << setw(7) << curr_node->bin_id;
        spath_f << setw(25) << curr_path->t << setw(30) << curr_path->k << setw(25) << curr_path->p << setw(25) << curr_path->s << endl;
        curr_node=shortest_paths[curr_node->node_id-1][0].prev_node;
    }
}

/* write the nodes that are visited along a shortest path to the sink (target) node to an output file */
void REA::print_visits(const Walker *walker) {
    cout << "print_visits()" << endl;
    ofstream visits_f;
    string visits_fname="visits."+to_string(walker->path_no)+".dat";
    visits_f.open(visits_fname,ios_base::trunc);
    int i=0;
    for (auto inpath: walker->visited) {
        if (inpath) visits_f << i << endl; i++; }
}
