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

REA::REA(const Network &ktn, bool discretetime, bool writerea, const Wrapper_args &wrapper_args) : Wrapper_Method(wrapper_args) {
    cout << "rea> using the REA to determine the " << wrapper_args.nabpaths << " highest-probability paths" << endl;
    this->discretetime=discretetime; this->writerea=writerea;
    source_node = *ktn.nodesB.begin(); // NB there is only a single source node
    sink_node = *ktn.nodesA.begin(); // NB there is only a single sink node
    shortest_paths.resize(ktn.n_nodes); candidate_paths.resize(ktn.n_nodes);
    for (int i=0;i<ktn.n_nodes;i++) { // allocate arrays for candidate and assigned shortest paths
        shortest_paths[i].resize(wrapper_args.nabpaths);
        for (int k=1;k<wrapper_args.nabpaths+1;k++) { // path cost is initially infinite and predecessor node not set for all paths
            /* for the REA Wrapper_Method, the members of the Walker objects are interpreted as follows:
               the [path_no]-th shortest path to curr_node, for which the path length/time/action/entropy are stored in members k/t/p/s,
               is the union of:    { [walker_id]-th shortest path to node prev_node } \cup curr_node
               and is stored as the element:    shortest_paths[curr_node->node_id-1][path_no-1] */
            shortest_paths[i][k-1] = { walker_id:0,path_no:k,k:0,t:0.L, \
                                       p:numeric_limits<long double>::infinity(),s:0.L,prev_node:nullptr,curr_node:&ktn.nodes[i] };
        }
        candidate_paths[i].resize(ktn.nodes[i].udeg); // edges are bidirectional, so in- and out-degrees of nodes are the same
        for (int j=0;j<ktn.nodes[i].udeg;j++) {
            candidate_paths[i][j].first=nullptr; candidate_paths[i][j].second=nullptr;
        }
    }
}

REA::~REA() {}

void REA::run_enhanced_kmc(const Network &ktn, Traj_Method *traj_method_obj) {

    dijkstra(ktn); // find shortest (i.e. highest-probability) path
    for (int k=2;k<nabpaths+1;k++) next_path(*sink_node,k); // main loop of REA
    for (int k=1;k<nabpaths+1;k++) shortest_paths[sink_node->node_id-1][k-1].dump_fpp_properties();
    if (writerea) print_shortest_paths(); // print the k shortest paths to the sink node
}

/* compute the first shortest path from the source node to all other nodes using Dijkstra's algorithm */
void REA::dijkstra(const Network& ktn) {
    if (debug) cout << "Dijkstra's algorithm" << endl;
    int m, n;
    vector<bool> insptree(ktn.n_nodes,false); // vector for bookkeeping which nodes have been incorporated into the shortest path tree
    // initialisation
    const Node *curr_node=source_node;
    shortest_paths[curr_node->node_id-1][0].p=0.L;
    // main loop for Dijkstra's algorithm
    for (int i=0;i<ktn.n_nodes;i++) {
        const Edge *edgeptr=curr_node->top_from;
        if (*curr_node==*sink_node) goto find_next_node; // sink_node cannot be a predecessor of any other node in the shortest path tree, skip
        n=curr_node->node_id-1;
        insptree[n]=true;
        while (edgeptr!=nullptr) { // loop over outgoing edges
            m=edgeptr->to_node->node_id-1;
            if (shortest_paths[n][0].p - 1.L*log(edgeptr->t) < shortest_paths[m][0].p) {
                // update path values
                shortest_paths[m][0].p = shortest_paths[n][0].p - 1.L*log(edgeptr->t);
                shortest_paths[m][0].k = shortest_paths[n][0].k + 1;
                shortest_paths[m][0].t = shortest_paths[n][0].t + edgeptr->from_node->t_esc;
                if (!discretetime) shortest_paths[m][0].s = shortest_paths[n][0].s + (edgeptr->rev_edge->k-edgeptr->k);
                shortest_paths[m][0].prev_node = curr_node; // set previous node in shortest path tree
                shortest_paths[m][0].walker_id = 1;
            }
            edgeptr=edgeptr->next_from;
        }
        find_next_node: {}
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
void REA::next_path(const Node &vnode, int k) {
    if (debug) cout << "next_path() for node: " << vnode.node_id << " path no: " << k << endl;
    const Edge *edgeptr;
    const Node *unode; int k1;
    /* second shortest path to node v, initialise a set of candidate paths based on the first shortest path tree */
    if (k==2) {
        // loop over nodes with edges to node v
        edgeptr = vnode.top_to;
        while (edgeptr!=nullptr) {
            unode = edgeptr->from_node;
            // we are interested in first passage paths; the sink node cannot be a predecessor in any shortest path to any node
            if (*unode==*sink_node) { edgeptr=edgeptr->next_to; continue; }
            if (vnode==*source_node || !(*(shortest_paths[vnode.node_id-1][0].prev_node)==*unode)) {
                /* from_node is not the predecessor of v in the shortest path tree, add the union of:
                   { 1st shortest path to node u } \cup node v    as a possible candidate for the next shortest path to node v */
                add_candidate(&shortest_paths[unode->node_id-1][0],edgeptr);
            }
            edgeptr=edgeptr->next_to;
        }
    }
    /* second shortest path to source node must have the form of the union of:
       { 1st shortest path to predecessor of source node } \cup source node,   all of which have already been added to the list of
       candidate paths, so skip to selection of best candidate path */
    if (vnode==*source_node && k==2) goto select_candidate; // skip to selection of candidate path
    /* find node u and path no. k1 that satisfy:  (k-1)-th shortest path to node v is the union of:
       { k1-th shortest path to node u } \cup node v */
    unode = shortest_paths[vnode.node_id-1][k-2].prev_node;
    k1 = shortest_paths[vnode.node_id-1][k-2].walker_id;
    if (debug) cout << "node u: " << unode->node_id << " k1: " << k1 << endl;
    // if the (k1+1)-th shortest path to node u has not yet been computed, then compute it with a recursive call to next_path()
    if (shortest_paths[unode->node_id-1][k1].prev_node==nullptr) next_path(*unode,k1+1);
    // at this point, the (k1+1)-th shortest path to node u should now have been determined */
    if (shortest_paths[unode->node_id-1][k1].prev_node==nullptr) {
        cout << "rea> error: failed to determine the " << k1+1 << "-th shortest path to node " << unode->node_id << endl; exit(EXIT_FAILURE); }
    // loop over edges from node u to find the u->v edge
    edgeptr = unode->top_from;
    while (edgeptr!=nullptr) {
        if (*(edgeptr->to_node)==vnode) break; edgeptr=edgeptr->next_from; }
    if (edgeptr==nullptr) {
        cout << "rea> error: there is not a direct transition from node " << unode->node_id << " to node " << vnode.node_id << endl; exit(EXIT_FAILURE); }
    /* add the union of:
       { (k1+1)-th shortest path to node u } \cup node v    as a possible candidate for the next shortest path to node v */
    add_candidate(&shortest_paths[unode->node_id-1][k1],edgeptr);
    // program skips directly to here after initialising list of candidate paths when vnode is the source node and k==2
    select_candidate: select_candidate(vnode,k); // find the candidate for the next (i.e. k-th) shortest path to node v with the lowest cost
}

/* add the union of: { k-th shortest path to node u } \cup node v,   where u/v are the to/from nodes associated with uvedge, respectively,
   and the corresponding path to node u is pointed to by cand_path,   as a possible candidate for the next shortest path to node v  */
void REA::add_candidate(const Walker *cand_path, const Edge *uvedge) {
    if (debug) {
        cout << "add_candidate() to node: " << uvedge->to_node->node_id \
             << "    path no. " << cand_path->path_no << " from node: " << uvedge->from_node->node_id << endl;
    }
    int v = uvedge->to_node->node_id;
    bool foundempty=false; // boolean value to indicate if an available space in the candidate_paths array has been found
    for (int i=0;i<uvedge->to_node->udeg;i++) {
        if (candidate_paths[v-1][i].first==nullptr) {
            candidate_paths[v-1][i].first=cand_path;
            candidate_paths[v-1][i].second=uvedge;
            foundempty=true;
            break;
        }
    }
    if (!foundempty) { cout << "rea> error: number of candidate paths exceeds max. possible number" << endl; exit(EXIT_FAILURE); }
}

/* select the best candidate (i.e. that with lowest cost) and assign as the k-th shortest path to node v, and remove the chosen
   candidate path from the list */
void REA::select_candidate(const Node &vnode, int k) {
    if (debug) cout << "in select_candidate for node: " << vnode.node_id << " path no.: " << k << endl;
    int v = vnode.node_id;
    int m=-1; long double mincost=numeric_limits<long double>::infinity();
    for (int i=0;i<vnode.udeg;i++) { // loop over candidate paths
        if (candidate_paths[v-1][i].first==nullptr) continue; // empty space in list of candidate paths
        if (debug) {
            cout << "  idx: " << i << endl; cout << "    path no. : " << candidate_paths[v-1][i].first->path_no \
                 << "    to node u: " << candidate_paths[v-1][i].first->curr_node->node_id << endl;
            cout << "    cost of  { path to u } cup v: " << candidate_paths[v-1][i].first->p - log(candidate_paths[v-1][i].second->t) << endl;
        }
        if (candidate_paths[v-1][i].first->p - log(candidate_paths[v-1][i].second->t) < mincost) {
            m=i; // m represents the index of the best candidate path in the list of all candidate paths for node v
            mincost = candidate_paths[v-1][i].first->p - log(candidate_paths[v-1][i].second->t);
        }
    }
    if (m<0) { cout << "rea> error: no candidates for next shortest path to node " << vnode.node_id << endl; exit(EXIT_FAILURE); }
    if (debug) cout << "  selected candidate idx m: " << m << "    mincost: " << mincost << endl;
    // assign the properties of the best candidate path to the k-th shortest path to node v
    shortest_paths[v-1][k-1].p = candidate_paths[v-1][m].first->p - log(candidate_paths[v-1][m].second->t);
    shortest_paths[v-1][k-1].k = candidate_paths[v-1][m].first->k + 1;
    shortest_paths[v-1][k-1].t = candidate_paths[v-1][m].first->t + candidate_paths[v-1][m].second->from_node->t_esc;
    if (!discretetime) {
        shortest_paths[v-1][k-1].s = candidate_paths[v-1][m].first->s + \
                (candidate_paths[v-1][m].second->rev_edge->k - candidate_paths[v-1][m].second->k);
    }
    shortest_paths[v-1][k-1].prev_node = candidate_paths[v-1][m].second->from_node;
    shortest_paths[v-1][k-1].walker_id = candidate_paths[v-1][m].first->path_no;
    // delete the selected candidate path from the list
    candidate_paths[v-1][m].first=nullptr; candidate_paths[v-1][m].second=nullptr;
}

/* print the k shortest paths to the sink node from the source node by tracing the elements in the array of the k shortest paths
   to all nodes (note that the paths are therefore printed backwards) */
void REA::print_shortest_paths() {
    if (debug) cout << "print_shortest_paths()" << endl;
    const Walker *walker;
    for (int k=1;k<nabpaths+1;k++) {
        ofstream spath_f;
        string spath_fname="shortest_path."+to_string(k)+".dat";
        spath_f.open(spath_fname,ios_base::trunc);
        spath_f.setf(ios::right,ios::adjustfield); spath_f.setf(ios::scientific,ios::floatfield);
        spath_f.precision(10);
        // start from k-th shortest path to sink node and loop to trace back through the k shortest paths array
        walker = &shortest_paths[sink_node->node_id-1][k-1];
        while (true) {
            // print path information
            spath_f << setw(7) << walker->curr_node->node_id << setw(7) << walker->curr_node->comm_id;
            spath_f << setw(25) << walker->t << setw(30) << walker->k << setw(25) << walker->p << setw(25) << walker->s << endl;
            if (walker->prev_node==nullptr) break;
            // find parent path of current path in k shortest paths array
            walker = &shortest_paths[walker->prev_node->node_id-1][walker->walker_id-1];
        }
    }
}
