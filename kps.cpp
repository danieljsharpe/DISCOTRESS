/*
File containing functions relating to kinetic path sampling
*/

#include "kmc_methods.h"
#include <queue>
#include <limits>
#include <cmath>
#include <iostream>

using namespace std;

KPS::KPS(const Network& ktn, int n_abpaths, int n_kpsmaxit, int nelim, double tau, int nbins, int kpskmcsteps, \
         bool adaptivebins, bool initcond, bool debug) {

    cout << "kps> running kPS with parameters:\n  lag time: " << tau << " \tmax. no. of eliminated nodes: " \
         << nelim << "\n  no. of bins: " << nbins << " \tno. of kMC steps after kPS iteration: " << kpskmcsteps \
         << "\n  adaptive binning (y/n): " << adaptivebins << endl;
    this->nelim=nelim; this->nbins=nbins; this->tau=tau; this->kpskmcsteps=kpskmcsteps;
    this->adaptivebins=adaptivebins; this->initcond=initcond;
    this->n_abpaths=n_abpaths; this->n_kpsmaxit=n_kpsmaxit; this->debug=debug;
    basin_ids.resize(ktn.n_nodes);
}

KPS::~KPS() {}

/* main loop of the kinetic path sampling algorithm */
void KPS::run_enhanced_kmc(const Network &ktn) {
    cout << "kps> beginning kPS simulation" << endl;
    int n_ab=0, n_kpsit=0;
    while ((n_ab<n_abpaths) and (n_kpsit<n_kpsmaxit)) {
        setup_basin_sets(ktn);
        graph_transformation(ktn);
        alpha = sample_absorbing_node();
        iterative_reverse_randomisation();
//        if (alpha->aorb==-1) n_ab++; // trajectory has reached endpoint absorbing macrostate A
        n_kpsit++;
        delete ktn_kps; delete ktn_kps_orig;
        delete ktn_l; delete ktn_u;
    }
    cout << "kps> finished kPS simulation" << endl;
}

/* Reset data of previous kPS iteration and find the microstates of the current trapping basin */
void KPS::setup_basin_sets(const Network &ktn) {

    cout << "kps> setting up basin sets" << endl;
    N_c=0; N=0; N_B=0; N_e=0;
    if (!alpha) { // first iteration of A-B path, need to set starting node
        if (!initcond) { // no initial condition was set, choose microstate in set B in proportion to stationary probabilities
            if (ktn.nodesB.size()==1000) { // quack
            auto it_set = ktn.nodesB.begin();
            const Node *tmpnodeptr=*it_set;
            epsilon = const_cast<Node*>(tmpnodeptr);
            } else {
            double pi_B = -numeric_limits<double>::infinity();
            set<Node*>::iterator it_set = ktn.nodesB.begin();
            while (it_set!=ktn.nodesB.end()) {
                pi_B = log(exp(pi_B)+exp((*it_set)->pi));
                it_set++; }
            vector<pair<Node*,double>> eps_probs(ktn.nodesB.size()); // accumulated probabilities of selecting starting node
            it_set = ktn.nodesB.begin();
            double cum_prob=0.;
            while (it_set!=ktn.nodesB.end()) {
                cum_prob += exp((*it_set)->pi-pi_B);
                eps_probs.push_back(make_pair((*it_set),cum_prob));
                it_set++; }
            double rand_no = KMC_Standard_Methods::rand_unif();
            vector<pair<Node*,double>>::iterator it_vec = eps_probs.begin();
            while (it_vec!=eps_probs.end()) {
                if ((*it_vec).second>=rand_no) epsilon=(*it_vec).first;
                it_vec++; }
            }
        } else {
            // ...
        }
    } else {
        epsilon=alpha;
    }
    alpha=nullptr;
    fill(basin_ids.begin(),basin_ids.end(),0); // reset basin IDs (zero flag indicates absorbing nonboundary node)
    if (!adaptivebins) { // basin IDs are based on community IDs
        // find all nodes of the current occupied pre-set community, mark these nodes as transient noneliminated nonboundary
        if (debug) cout << "basin nodes:" << endl;
        for (int i=0;i<ktn.n_nodes;i++) {
            if (ktn.nodes[i].comm_id==epsilon->comm_id) {
                if (debug) cout << "  " << i+1;
                basin_ids[i]=3; N_B++; N_e+=ktn.nodes[i].udeg; }
        }
        if (debug) cout << endl << "absorbing nodes:" << endl;
        // find all absorbing boundary nodes
        for (int i=0;i<ktn.n_nodes;i++) {
            if (basin_ids[i]!=3) continue;
            Edge *edgeptr = ktn.nodes[i].top_from;
            while (edgeptr!=nullptr) {
                if (edgeptr->deadts) {edgeptr=edgeptr->next_from; continue; }
                if (edgeptr->to_node->comm_id!=epsilon->comm_id && !basin_ids[edgeptr->to_node->node_id-1]) {
                    basin_ids[edgeptr->to_node->node_id-1]=1;
                    N_e+=ktn.nodes[edgeptr->to_node->node_id-1].udeg; N_c++;
                    if (debug) cout << "  " << edgeptr->to_node->node_id;
                }
                edgeptr=edgeptr->next_from;
            }
        }
        cout << endl;
    } else {
        // ...
    }
    eliminated_nodes.clear(); nodemap.clear();
    eliminated_nodes.reserve(!(N_B>nelim)?N_B:nelim);
    if (debug) {
        cout << "number of eliminated nodes: " << (!(N_B>nelim)?N_B:nelim) << endl;
        cout << "number of nodes in basin: " << N_B << " number of absorbing boundary nodes: " << N_c << endl;
        cout << "number of edges of subnetwork: " << N_e << endl;
        cout << "epsilon: " << epsilon->node_id << endl;
        cout << "currently occupied community id: " << epsilon->comm_id << endl;
    }
}

/* Iterative reverse randomisation procedure to stochastically sample the hopping matrix
   H^(0) corresponding to T^(0), given H^(N) and the {T^(n)} for 0 <= n <= N.
   Return a sampled time for the stochastic escape trajectory. */
double KPS::iterative_reverse_randomisation() {

    cout << "kps> iterative reverse randomisation" << endl;
    cout << "N is: " << N << endl;
    for (int i=N;i>0;i--) {
//        cout << "Node: " << eliminated_nodes[i-1]->node_id << " k: " << eliminated_nodes[i-1]->top_from->k << endl;
    }
    double t_esc = KPS::gamma_distribn(0,0.); // time for escape trajectory
    return t_esc;
}

/* Sample a node at the absorbing boundary of the current trapping basin, by the
   categorical sampling procedure based on T^(0) and T^(N) */
Node *KPS::sample_absorbing_node() {

    cout << "kps> sample absorbing node" << endl;
    cout << "epsilon: " << epsilon->node_id << endl;
    int curr_comm_id = epsilon->comm_id;
    Node *next_node, *curr_node;
    curr_node = epsilon;
    do {
    if (basin_ids[epsilon->node_id=1]==1) { // eliminated node
    Edge *edgeptr = epsilon->top_from;
    Node *next_node;
    while (edgeptr!=nullptr) {

        edgeptr=edgeptr->next_from;
    }
    }
    else if (basin_ids[epsilon->node_id-1]==3) { // transient noneliminated node

    }
    else {
        cout << "kps> something went wrong in sample_absorbing_node()" << endl; exit(EXIT_FAILURE); }
    // set curr_node to next_node
    break;
    } while (epsilon->comm_id==curr_comm_id);
    // set alpha to curr_node
    return nullptr;
}

/* Graph transformation to eliminate up to N nodes of the current trapping basin.
   Calculates the set of N-1 transition probability matrices {T^(n)} for 0 < n <= N.
   The transition network input to this function is the full network, and get_subnetwork() returns T^(0).
   The graph transformation is performed by performing a LU-decomposition of T^(0) */
void KPS::graph_transformation(const Network &ktn) {

    if (debug) cout << "kps> graph transformation" << endl;
    ktn_kps=get_subnetwork(ktn);
    ktn_kps_orig=get_subnetwork(ktn);
    ktn_l = new Network(N_B+N_c,0);
    ktn_u = new Network(N_B+N_c,0);
    for (int i=0;i<ktn_kps->n_nodes;i++) {
        ktn_l->nodes[i] = ktn_kps->nodes[i];
        ktn_u->nodes[i] = ktn_kps->nodes[i];
        ktn_l->nodes[i].t=0.; ktn_u->nodes[i].t=0.; // the "transn probs" in the L and U TNs are the values to "undo" GT
    }
    auto cmp = [](Node *l,Node *r) { return l->udeg >= r->udeg; };
    priority_queue<Node*,vector<Node*>,decltype(cmp)> gt_pq(cmp); // priority queue of nodes (based on out-degree)
    for (vector<Node>::iterator it_nodevec=ktn_kps->nodes.begin();it_nodevec!=ktn_kps->nodes.end();++it_nodevec) {
        if (it_nodevec->comm_id!=epsilon->comm_id) continue;
        gt_pq.push(&(*it_nodevec));
    }
    h = vector<int>(gt_pq.size(),0); // reset flicker vector
    while (!gt_pq.empty() && N<nelim) {
        Node *node_elim=gt_pq.top();
//        cout << "N: " << N << " Node: " << node_elim->node_id << " priority: " << node_elim->udeg \
               << " k: " << node_elim->top_from->k << endl;
        gt_pq.pop();
        node_elim = &ktn_kps->nodes[N]; // eliminate nodes in order of IDs
        gt_iteration(node_elim);
        eliminated_nodes.push_back(node_elim);
        N++;
    }
    if (N!=(!(N_B>nelim)?N_B:nelim)) {
        cout << "kps> fatal error: lost track of number of eliminated nodes" << endl; exit(EXIT_FAILURE); }
    if (debug) cout << "kps> finished graph transformation" << endl;
}

/* return the subnetwork corresponding to the active trapping basin and absorbing boundary nodes, to be transformed
   in the graph transformation phase of the kPS algorithm */
Network *KPS::get_subnetwork(const Network& ktn) {

    if (debug) cout << "kps> get_subnetwork: create TN of " << N_B+N_c << " nodes and " << N_e << " edges" << endl;
    Network *ktnptr = new Network(N_B+N_c,N_e);
    ktnptr->edges.resize(N_e); // quack why is this necessary?
    int j=0;
    for (int i=0;i<ktn.n_nodes;i++) {
        if (!basin_ids[i]) continue;
        nodemap[i]=j; j++;
        ktnptr->nodes[j-1] = ktn.nodes[i];
    }
    int m=0, n=0;
    vector<bool> edgemask(2*ktn.n_edges,false);
    // note that the indices of the edge array in the subnetwork are not in a meaningful order
    for (map<int,int>::iterator it_map=nodemap.begin();it_map!=nodemap.end();++it_map) {
        n++;
        const Node *node = &ktn.nodes[it_map->first];
        // for absorbing node, do not incl any FROM edges, or any TO edges for non-basin nbr nodes, in the subnetwork
        if (node->comm_id!=epsilon->comm_id) continue;
        const Edge *edgeptr = node->top_from;
        while (edgeptr!=nullptr) {
            if (edgeptr->deadts || edgemask[edgeptr->edge_pos]) {edgeptr=edgeptr->next_from; continue; }
            ktnptr->edges[m] = *edgeptr; // edge of subnetwork inherits properties (transn rate etc) of node in full network
            ktnptr->edges[m].edge_pos = m;
            ktnptr->edges[m].from_node = &ktnptr->nodes[nodemap[edgeptr->from_node->node_id-1]];
            ktnptr->edges[m].to_node = &ktnptr->nodes[nodemap[edgeptr->to_node->node_id-1]];
            ktnptr->add_from_edge(nodemap[edgeptr->from_node->node_id-1],m);
            ktnptr->add_to_edge(nodemap[edgeptr->to_node->node_id-1],m);
//            Network::add_edge_network(ktnptr,ktnptr->nodes[nodemap[edgeptr->from_node->node_id-1]], \
                ktnptr->nodes[nodemap[edgeptr->to_node->node_id-1]],m);
            m++; edgemask[edgeptr->edge_pos]=true;
            const Edge *edgeptr_rev = edgeptr->rev_edge;
            if (edgeptr_rev->deadts || edgemask[edgeptr_rev->edge_pos]) {
                edgeptr=edgeptr->next_from; continue; }
            // reverse edge
            ktnptr->edges[m] = *edgeptr_rev;
            ktnptr->edges[m].edge_pos = m;
            ktnptr->edges[m].from_node = &ktnptr->nodes[nodemap[edgeptr_rev->from_node->node_id-1]];
            ktnptr->edges[m].to_node = &ktnptr->nodes[nodemap[edgeptr_rev->to_node->node_id-1]];
            ktnptr->add_from_edge(nodemap[edgeptr_rev->from_node->node_id-1],m);
            ktnptr->add_to_edge(nodemap[edgeptr_rev->to_node->node_id-1],m);
//            Network::add_edge_network(ktnptr,ktnptr->nodes[nodemap[edgeptr_rev->from_node->node_id-1]], \
                ktnptr->nodes[nodemap[edgeptr_rev->to_node->node_id-1]],m);
            ktnptr->edges[m-1].rev_edge = &ktnptr->edges[m];
            ktnptr->edges[m].rev_edge = &ktnptr->edges[m-1];
            m++; edgemask[edgeptr_rev->edge_pos]=true;
        }
    }
    cout << "added " << n << " nodes and " << m << " edges to subnetwork" << endl;
    if (n!=N_B+N_c || m!=N_e) { cout << "kps> something went wrong in get_subnetwork()" << endl; exit(EXIT_FAILURE); }
//    Network *ktnptr = &const_cast<Network&>(ktn); // quack
    return ktnptr;
}

/* a single iteration of the graph transformation method. Argument is a pointer to the node to be
   eliminated from the network to which the ktn_kps pointer refers.
   The networks "L" and "U" required to undo the graph transformation iterations are updated */
void KPS::gt_iteration(Node *node_elim) {

    double factor=0.; // equal to (1-T_{nn})
    if (node_elim->t>0.999) { // loop over neighbouring edges to maintain numerical precision
        Edge *edgeptr = node_elim->top_from;
        while (edgeptr!=nullptr) {
            if (!(edgeptr->deadts || edgeptr->to_node->eliminated)) factor += edgeptr->t;
            edgeptr=edgeptr->next_from;
        }
    } else { factor=1.-node_elim->t; }
    cout << "eliminating node: " << node_elim->node_id << endl;
    // objects to queue all non-eliminated neighbours of the eliminated node, and update relevant edges
    vector<Node*> nodes_nbrs;
    typedef struct {
        bool dirconn; // flag indicates if node is directly connected to current node being considered
        double t_fromn; // transition probability from eliminated node to this node
        double t_ton; // transition probability to eliminated node from this node
    } nbrnode;
    map<int,nbrnode> nbrnode_map; // map contains all non-elimd nodes directly connected to the current elimd node;
    // update the weights for all edges from the elimd node to non-elimd nbr nodes, and self-loops of nbr nodes
    Edge *edgeptr = node_elim->top_from;
    while (edgeptr!=nullptr) {
        if (edgeptr->deadts || edgeptr->to_node->eliminated) {
            edgeptr=edgeptr->next_from; continue; }
        edgeptr->to_node->flag=true;
        nodes_nbrs.push_back(edgeptr->to_node); // queue non-eliminated neighbour node
        nbrnode_map[edgeptr->to_node->node_id]=(nbrnode){false,edgeptr->t,edgeptr->rev_edge->t};
        cout << "  to node: " << edgeptr->to_node->node_id << endl;
        // update L and U networks here...
        cout << "    old node t: " << edgeptr->to_node->t << "  incr in node t: " \
             << (edgeptr->t)*(edgeptr->rev_edge->t)/factor << endl;
        edgeptr->to_node->t += (edgeptr->t)*(edgeptr->rev_edge->t)/factor; // update self-loop of non-elimd nbr node
        // update L and U networks here...
        cout << "    old edge t: " << edgeptr->t << "  incr in t: " << (edgeptr->t)*(node_elim->t)/factor << endl;
        edgeptr->t += (edgeptr->t)*(node_elim->t)/factor; // update edge from elimd node to non-elimd nbr node
        edgeptr=edgeptr->next_from;
    }
    cout << "  updating edges between pairs of nodes both directly connected to the eliminated node..." << endl;
    // update the weights for all pairs of nodes directly connected to the eliminated node
    int old_n_edges = ktn_kps->n_edges; // number of edges in the network before we start adding edges in the GT algorithm
    for (vector<Node*>::iterator it_nodevec=nodes_nbrs.begin();it_nodevec!=nodes_nbrs.end();++it_nodevec) {
        edgeptr = (*it_nodevec)->top_from; // loop over neighbouring nodes of this non-eliminated, non-absorbing node
        while (edgeptr!=nullptr) {
            // skip nodes not directly connected to elimd node, and also skip elimd nodes
            if (edgeptr->deadts || edgeptr->to_node->eliminated || !edgeptr->to_node->flag) {
                edgeptr=edgeptr->next_from; continue; }
            cout << "  node " << (*it_nodevec)->node_id << " is directly connected to node " << edgeptr->to_node->node_id << endl;
            nbrnode_map[edgeptr->to_node->node_id].dirconn=true; // this pair of elimd nodes are directly connected
            if (edgeptr->edge_pos>old_n_edges) { // skip nodes for which a new edge has already been added
                cout << "    edge already added" << endl;
                edgeptr=edgeptr->next_from; continue; }
            // update L and U networks here...
            cout << "    old edge t: " << edgeptr->t << "  incr in t: " \
                 << (nbrnode_map[edgeptr->from_node->node_id].t_ton)*\
                    (nbrnode_map[edgeptr->to_node->node_id].t_fromn)/factor << endl;
            edgeptr->t += (nbrnode_map[edgeptr->from_node->node_id].t_ton)*\
                (nbrnode_map[edgeptr->to_node->node_id].t_fromn)/factor;
            edgeptr=edgeptr->next_from;
        }
        for (map<int,nbrnode>::iterator it_map=nbrnode_map.begin();it_map!=nbrnode_map.end();++it_map) {
            if (it_map->first==(*it_nodevec)->node_id) continue;
            if ((it_map->second).dirconn) { (it_map->second).dirconn=false; continue; } // reset flag
            cout << "  node " << (*it_nodevec)->node_id << " is not directly connected to node " << it_map->first << endl;
            cout << "  t of new edge: " << (it_map->second).t_fromn*nbrnode_map[(*it_nodevec)->node_id].t_ton/factor << endl;
            cout << "  t of new reverse edge: " << (it_map->second).t_ton*nbrnode_map[(*it_nodevec)->node_id].t_fromn/factor << endl;
            // update L and U networks here...
            // nodes are directly connected to the elimd node but not to one another, add an edge in the transformed network
            (ktn_kps->edges).push_back(Edge());
            (ktn_kps->edges).back().t = (it_map->second).t_fromn*nbrnode_map[(*it_nodevec)->node_id].t_ton/factor;
            size_t pos = (ktn_kps->edges).size()-1;
            (ktn_kps->edges).back().edge_pos = pos;
//            Network::add_edge_network(ktn_kps,&ktn_kps->nodes[nodemap[edgeptr->from_node->node_id-1]], \
                &ktn_kps->nodes[nodemap[edgeptr->to_node->node_id-1]],pos);
            (ktn_kps->edges).back().from_node = &ktn_kps->nodes[nodemap[(*it_nodevec)->node_id-1]];
            (ktn_kps->edges).back().to_node = &ktn_kps->nodes[nodemap[it_map->first-1]];
            ktn_kps->add_from_edge(nodemap[(*it_nodevec)->node_id-1],pos);
            ktn_kps->add_to_edge(nodemap[it_map->first-1],pos);
            ktn_kps->n_edges++; // increment number of edges in the network
            // reverse edge
            (ktn_kps->edges).push_back(Edge());
            (ktn_kps->edges).back().t = (it_map->second).t_ton*nbrnode_map[(*it_nodevec)->node_id].t_fromn/factor;
            (ktn_kps->edges).back().edge_pos = pos+1;
            (ktn_kps->edges).back().from_node = &ktn_kps->nodes[nodemap[it_map->first-1]];
            (ktn_kps->edges).back().to_node = &ktn_kps->nodes[nodemap[(*it_nodevec)->node_id-1]];
            ktn_kps->add_from_edge(nodemap[it_map->first-1],pos+1);
            ktn_kps->add_to_edge(nodemap[(*it_nodevec)->node_id-1],pos+1);
//            Network::add_edge_network(ktn_kps,&ktn_kps->nodes[nodemap[edgeptr_rev->from_node->node_id-1]], \
                &ktn_kps->nodes[nodemap[edgeptr_rev->to_node->node_id-1]],pos+1);
            ktn_kps->n_edges++;
            ktn_kps->edges[pos+1].rev_edge = &ktn_kps->edges[pos];
            ktn_kps->edges[pos].rev_edge = &ktn_kps->edges[pos+1];
        }
    }
    // reset the flags
    edgeptr = node_elim->top_from;
    while (edgeptr!=nullptr) {
        edgeptr->to_node->flag=false;
        edgeptr = edgeptr->next_from;
    }
    node_elim->eliminated=true; // this flag negates the need to zero the weights to the eliminated node
}

/* undo a single iteration of the graph transformation. Argument is a pointer to the node to be un-eliminated from the network */
void KPS::undo_gt_iteration(Node *node_elim) {

}

/* Gamma distribution with shape parameter a and rate parameter 1./b */
double KPS::gamma_distribn(int a, double b) {

    return 1.;
}

/* Binomial distribution with trial number h and success probability p.
   Returns the number of successes after h Bernoulli trials. */
int KPS::binomial_distribn(int h, double p) {

    if (!((h>=0) || ((p>=0.) && (p<=1.)))) throw exception();
    if (h==0)  { return 0;
    } else if (p==1.) { return h; }
    // ...
    return 1;
}

/* Negative binomial distribution with success number r and success probability p.
   Returns the number of failures before the r-th success. */
int KPS::negbinomial_distribn(int r, double p) {

    if (!((r>=0) || ((p>0.) && (p<=1.)))) throw exception();
    if ((r==0) || (p==1.)) return 0;
    // ...
    return 1;
}

/* Exponential distribution with rate parameter 1./tau */
double KPS::exp_distribn(double tau) {

    return 1.;
}
