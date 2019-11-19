/* tests and sanity checks to aid debugging */

#include "ktn.h"
#include <iostream>

using namespace std;

void run_debug_tests(const Network& ktn) {

    // check implementation of transition network data structure
    int test_node=22, test_edge=30; // indices of node and edge to print info for
    cout << "\ndebug> printing info for node " << test_node << endl;
    cout << "  node id: " << ktn.nodes[test_node-1].node_id << endl;
    cout << "  community id: " << ktn.nodes[test_node-1].comm_id << endl;
    cout << "  aorb endpoint set flag: " << ktn.nodes[test_node-1].aorb << endl;
    cout << "  log stationary probability: " << ktn.nodes[test_node-1].pi << endl;
    cout << "  log escape rate: " << ktn.nodes[test_node-1].k_esc << endl;
    cout << "  self-transition probability: " << ktn.nodes[test_node-1].t << endl;

    cout << "\ndebug> printing outgoing edges for node " << test_node << endl;
    Edge *edgeptr = ktn.nodes[test_node-1].top_to;
    while (edgeptr!=nullptr) {
        cout << "  ts id: " << edgeptr->ts_id << " edge id: " << edgeptr->edge_pos << endl;
        cout << "    from node: " << edgeptr->from_node->node_id << " to node: " << edgeptr->to_node->node_id << endl;
        cout << "    log transition rate: " << edgeptr->k << endl;
        cout << "    transition probability: " << edgeptr->k << endl;
        cout << "    flux: " << edgeptr->j << endl;
        edgeptr=edgeptr->next_to;
    }
}
