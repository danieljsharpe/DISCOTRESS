/* tests and sanity checks to aid debugging */

#ifndef __DEBUG_TESTS_H_INCLUDED__
#define __DEBUG_TESTS_H_INCLUDED__

#include "ktn.h"
#include <iostream>

using namespace std;

void run_debug_tests(const Network& ktn) {

    // check implementation of transition network data structure
    cout << "\ndebug> ktn info: no. of nodes: " << ktn.n_nodes << " no. of edges: " << ktn.n_edges << endl;
    if (ktn.n_nodes<20) {
    cout << "debug> network is small, printing all transition probabilities" << endl;
    for (int i=0;i<ktn.n_nodes;i++) {
        cout << "node: " << i+1 << endl;
        cout << "  to: " << i+1 << "  t: " << ktn.nodes[i].t << endl;
        Edge *edgeptr = ktn.nodes[i].top_from;
        while (edgeptr!=nullptr) {
            if (!edgeptr->deadts) cout << "  to: " << edgeptr->to_node->node_id << "  t: " << edgeptr->t << endl;
            edgeptr = edgeptr->next_from;
        }
    }
    }
    int test_node=21; // index of node to print info for
    if (test_node>ktn.n_nodes) {
        cout << "debug> skipping tests, node with index " << test_node << " does not exist" << endl; return; }
    cout << "debug> printing info for node " << test_node << endl;
    cout << "  node id: " << ktn.nodes[test_node-1].node_id << endl;
    cout << "  community id: " << ktn.nodes[test_node-1].comm_id << endl;
    cout << "  aorb endpoint set flag: " << ktn.nodes[test_node-1].aorb << endl;
    cout << "  log stationary probability: " << ktn.nodes[test_node-1].pi << endl;
    cout << "  log escape rate: " << ktn.nodes[test_node-1].k_esc << " escape rate: " << exp(ktn.nodes[test_node-1].k_esc) << endl;
    cout << "  self-transition probability: " << ktn.nodes[test_node-1].t << endl;

    cout << "\ndebug> printing outgoing edges for node " << test_node << endl;
    Edge *edgeptr = ktn.nodes[test_node-1].top_from;
    Edge *edgeptr_prev;
    while (edgeptr!=nullptr) {
        cout << "\n  ts id: " << edgeptr->ts_id << " edge id: " << edgeptr->edge_pos << endl;
        cout << "    from node: " << edgeptr->from_node->node_id << " to node: " << edgeptr->to_node->node_id << endl;
        cout << "    log transition rate: " << edgeptr->k << " transition rate: " << exp(edgeptr->k) << endl;
        cout << "    transition probability: " << edgeptr->t << endl;
        cout << "    flux: " << edgeptr->j << endl;
        edgeptr_prev = edgeptr;
        edgeptr = edgeptr->rev_edge;
        cout << "  reverse edge..." << endl;
        cout << "  ts id: " << edgeptr->ts_id << " edge id: " << edgeptr->edge_pos << endl;
        cout << "    from node: " << edgeptr->from_node->node_id << " to node: " << edgeptr->to_node->node_id << endl;
        cout << "    log transition rate: " << edgeptr->k << " transition rate: " << exp(edgeptr->k) << endl;
        cout << "    transition probability: " << edgeptr->t << endl;
        cout << "    flux: " << edgeptr->j << endl;
        edgeptr=edgeptr_prev->next_from;
    }
}

#endif
