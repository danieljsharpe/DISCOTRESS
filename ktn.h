/*
Custom data structure for representing and reading in a ktn
*/

#ifndef __KTN_H_INCLUDED__
#define __KTN_H_INCLUDED__

#include <set>
#include <exception>
#include <vector>

#include <iostream>

using namespace std;

class KMC_Suite;

struct Node;

struct Edge {
    int ts_id;
    int edge_pos; // position of the TS in the edges array
    double k; // (log) transition rate
    double t; // transition probability
    double j; // net flux
    bool deadts; // indicates TS only linked to one minimum or is otherwise deleted
    Node *to_node=nullptr;
    Node *from_node=nullptr;
    Edge *next_to=nullptr;
    Edge *next_from=nullptr;
    Edge *rev_edge=nullptr; // reverse edge (all edges are bidirectional)

    inline Edge operator+(const Edge& other_edge) const {
        Edge new_edge{.ts_id=ts_id, .edge_pos=edge_pos, .k=k+other_edge.k, \
            .t=t+other_edge.t, .j=j+other_edge.j, .deadts=deadts, .to_node=to_node, \
            .from_node=from_node, .next_to=next_to, .next_from=next_from, .rev_edge=rev_edge};
        return new_edge;
    }

    inline Edge& operator=(const Edge& other_edge) {
        ts_id=other_edge.ts_id;
        k=other_edge.k; t=other_edge.t; j=other_edge.j; deadts=other_edge.deadts;
    }
};

struct Node {
    int node_id;
    int comm_id = -1; // community ID (-1 indicates null value)
    int aorb = 0; // indicates set to which node belongs: -1 for A, +1 for B, 0 for I
    int udeg = 0; // (unweighted) node (out-) degree
    bool eliminated = false; // node has been eliminated from the network (in graph transformation) (or otherwise deleted)
    double k_esc; // (log) escape rate from node (sum of outgoing transition rates)
    double t; // self-transition probability
    double pi; // (log) occupation probability (usually the stationary/equilibrium probability)
    Edge *top_to=nullptr;
    Edge *top_from=nullptr;

    Node();
    ~Node();
    Node(const Node&);

    const inline bool operator<(const Node& other_node) const {
        return (node_id<other_node.node_id);
    }

    /* in assignment operator for Node, do not copy the pointers to Edge objects */
    inline Node& operator=(const Node& other_node) {
        node_id=other_node.node_id; comm_id=other_node.comm_id; aorb=other_node.aorb;
        udeg=other_node.udeg; eliminated=other_node.eliminated;
        k_esc=other_node.k_esc; t=other_node.t; pi=other_node.pi;
    }
};

// structure containing the kinetic transition network
struct Network {

    friend class KMC_Suite;

    public:

    Network(int,int);
    ~Network();
    Network(const Network&);

    void del_node(int);
    void add_to_edge(int,int);
    void add_from_edge(int,int);
    void del_to_edge(int);
    void del_from_edge(int);
    void del_spec_to_edge(int,int);
    void del_spec_from_edge(int,int);
    void update_to_edge(int,int);
    void update_from_edge(int,int);
    static void calc_k_esc(Node&);
    static void calc_net_flux(Edge&);
    void get_tmtx_lin(double); // calculate the linearised transition probability matrix
    static void setup_network(Network&,const vector<pair<int,int>>&,const vector<double>&, \
        const vector<double>&,const vector<int>&, const vector<int>&,const vector<int>& = {});

    vector<Node> nodes;
    vector<Edge> edges; // note that this vector contains two entries for forward and reverse transitions for
                        // each pair of nodes

    struct Ktn_exception {
        const char * what () const throw () { return "KTN Exception"; }
    };

    int n_nodes, n_edges; // number of nodes and bidirectional edges (not including self-loops)
    int tot_nodes=0, tot_edges=0;
    int n_dead=0; // number of dead/deleted edges
    int ncomms; // total number of communities
    set<Node*> nodesA, nodesB; // A and B endpoint nodes (A<-B)

    inline Network& operator=(const Network& other_network) {
        nodes=other_network.nodes; edges=other_network.edges;
        n_nodes=other_network.n_nodes; n_edges=other_network.n_edges;
        n_dead=other_network.n_dead; ncomms=other_network.ncomms;
        nodesA=other_network.nodesA; nodesB=other_network.nodesB;
        return *this;
    }

};

#endif
