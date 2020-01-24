/*
File containing functions relating to kinetic path sampling
*/

#include "kmc_methods.h"
#include <queue>
#include <cmath>
#include <random>
#include <iostream>

using namespace std;

KPS::KPS(const Network &ktn, int maxn_abpaths, int maxit, int nelim, double tau, int nbins, int kpskmcsteps, \
         bool adaptivebins, bool initcond, int seed, bool debug) {

    cout << "kps> running kPS with parameters:\n  lag time: " << tau << " \tmax. no. of eliminated nodes: " \
         << nelim << "\n  no. of bins: " << nbins << " \tno. of kMC steps after kPS iteration: " << kpskmcsteps \
         << "\n  adaptive binning (y/n): " << adaptivebins \
         << "\n  random seed: " << seed << " \tdebug printing: " << debug << endl;
    this->nelim=nelim; this->nbins=nbins; this->tau=tau; this->kpskmcsteps=kpskmcsteps;
    this->adaptivebins=adaptivebins; this->initcond=initcond;
    this->maxn_abpaths=maxn_abpaths; this->maxit=maxit; this->seed=seed; this->debug=debug;
    basin_ids.resize(ktn.n_nodes);
    // quack need to move this somewhere more general
    if (nbins>0) {
        visited.resize(nbins); fill(visited.begin(),visited.end(),false);
        tp_visits.resize(nbins); tp_densities.resize(nbins);
        ab_successes.resize(nbins); ab_failures.resize(nbins); committors.resize(nbins);
    }
}

KPS::~KPS() {}

void KPS::test_ktn(const Network &ktn) {
    cout << "debug> ktn info: no. of nodes: " << ktn.n_nodes << " no. of edges: " << ktn.n_edges << endl;
    for (int i=0;i<ktn.n_nodes;i++) {
        cout << "node: " << i+1 << endl;
        if (!ktn.nodes[i].eliminated) cout << "  to: " << i+1 << "  t: " << ktn.nodes[i].t << endl;
        Edge *edgeptr = ktn.nodes[i].top_from;
        while (edgeptr!=nullptr) {
            if (!edgeptr->deadts && !edgeptr->to_node->eliminated) cout << "  to: " << edgeptr->to_node->node_id << "  t: " << edgeptr->t \
                << "  pos: " << edgeptr->edge_pos << endl;
            edgeptr = edgeptr->next_from;
        }
    }
    cout << "\n" << endl;
}

/* main loop of the kinetic path sampling algorithm */
void KPS::run_enhanced_kmc(const Network &ktn) {
    cout << "\nkps> beginning kPS simulation" << endl;
    n_ab=0; n_traj=0; int n_kpsit=0;
    while ((n_ab<maxn_abpaths) and (n_kpsit<maxit)) { // algorithm terminates when max no of kPS basin escapes have been simulated
        setup_basin_sets(ktn);
        graph_transformation(ktn);
        Node *dummy_alpha = sample_absorbing_node();
        alpha = &ktn.nodes[dummy_alpha->node_id-1];
        visited[alpha->comm_id]=true;
        double t_esc = iterative_reverse_randomisation();
        update_path_quantities(t_esc,alpha);
        walker.dump_walker_info(n_traj,alpha->aorb==-1);
        if (alpha->aorb==-1 || alpha->aorb==1) { // traj has reached absorbing macrostate A or has returned to B
            update_tp_stats(alpha->aorb==-1,!adaptivebins);
            walker.reset_walker_info();
        }
        n_kpsit++;
        epsilon=alpha; alpha=nullptr;
        delete ktn_kps; delete ktn_kps_orig;
        delete ktn_l; delete ktn_u;
    }
    cout << "\nkps> kPS simulation terminated after " << n_kpsit << " iterations. Simulated " << n_ab << " transition paths" << endl;
    if (!adaptivebins) calc_tp_stats(); // calculate committors and transn path densities for communities and write to file
    cout << "kps> walker time: " << walker.t << "  activity: " << walker.k \
         << "  log path prob: " << walker.p << "  entropy flow: " << walker.s << endl;
    cout << "kps> finished kPS simulation" << endl;
}

/* Reset data of previous kPS iteration and find the microstates of the current trapping basin */
void KPS::setup_basin_sets(const Network &ktn) {

    if (debug) cout << "\nkps> setting up basin sets" << endl;
    N_c=0; N=0; N_B=0; N_e=0;
    double pi_B = -numeric_limits<double>::infinity(); // (log) occupation probability of all nodes in initial set B
    if (!epsilon) { // first iteration of A<-B path, need to set starting node
        if (!initcond) { // no initial condition was set, choose node in set B in proportion to stationary probs
            if (ktn.nodesB.size()==1000) { // quack should be ==1
            auto it_set = ktn.nodesB.begin();
            epsilon=*it_set;
            } else {
            set<Node*>::iterator it_set = ktn.nodesB.begin();
            while (it_set!=ktn.nodesB.end()) {
                pi_B = log(exp(pi_B)+exp((*it_set)->pi));
                it_set++; }
            vector<pair<Node*,double>> eps_probs(ktn.nodesB.size()); // accumulated probs of selecting starting node
            it_set = ktn.nodesB.begin();
            double cum_prob=0.;
            while (it_set!=ktn.nodesB.end()) {
                cum_prob += exp((*it_set)->pi-pi_B);
                eps_probs.push_back(make_pair((*it_set),cum_prob));
                it_set++; }
            double rand_no = KMC_Standard_Methods::rand_unif_met(seed);
            vector<pair<Node*,double>>::iterator it_vec = eps_probs.begin();
            while (it_vec!=eps_probs.end()) {
                if ((*it_vec).second>=rand_no) { epsilon=(*it_vec).first; break; }
                it_vec++; }
            }
        } else { // choose node in set B in proportion to specified initial condition probs
            // ...
        }
    }
    fill(basin_ids.begin(),basin_ids.end(),0); // reset basin IDs (zero flag indicates absorbing nonboundary node)
    if (!adaptivebins) { // basin IDs are based on community IDs
        // find all nodes of the current occupied pre-set community, mark these nodes as transient noneliminated
        if (debug) cout << "basin nodes:" << endl;
        for (int i=0;i<ktn.n_nodes;i++) {
            if (ktn.nodes[i].comm_id==epsilon->comm_id) {
                if (debug) cout << "  " << i+1;
                basin_ids[i]=2; N_B++; N_e+=ktn.nodes[i].udeg; }
        }
        if (debug) cout << endl << "absorbing nodes:" << endl;
        // find all absorbing boundary nodes
        for (int i=0;i<ktn.n_nodes;i++) {
            if (basin_ids[i]!=2) continue;
            Edge *edgeptr = ktn.nodes[i].top_from;
            while (edgeptr!=nullptr) {
                if (edgeptr->deadts) {edgeptr=edgeptr->next_from; continue; }
                if (edgeptr->to_node->comm_id!=epsilon->comm_id && !basin_ids[edgeptr->to_node->node_id-1]) {
                    basin_ids[edgeptr->to_node->node_id-1]=3; // flag absorbing boundary node
                    N_e+=ktn.nodes[edgeptr->to_node->node_id-1].udeg; N_c++;
                    if (debug) cout << "  " << edgeptr->to_node->node_id;
                }
                edgeptr=edgeptr->next_from;
            }
        }
        if (debug) cout << endl;
    } else {
        // ...
    }
    eliminated_nodes.clear(); nodemap.clear();
    eliminated_nodes.reserve(!(N_B>nelim)?N_B:nelim);
    visited[epsilon->comm_id]=true;
    walker.curr_node=&(*epsilon);
    walker.p=epsilon->pi-pi_B; // factor in path probability corresponding to probability of initial node
    walker.dump_walker_info(n_traj,false);
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

    if (debug) {
        cout << "\nkps> iterative reverse randomisation" << endl;
        cout << "N is: " << N << endl; cout << "node alpha: " << alpha->node_id << endl; }
    // update kMC hop counts for noneliminated nodes
    for (auto &node: ktn_kps->nodes) {
        if (node.eliminated) continue;
        node.h += KPS::negbinomial_distribn(node.h,1.-ktn_kps_orig->nodes[node.node_id-1].t,seed); }
    // main loop of the iterative reverse randomisation procedure
    for (int i=N;i>0;i--) {
        Node *curr_node = eliminated_nodes[i-1];
        vector<pair<Node*,Edge*>> nodes_nbrs = undo_gt_iteration(curr_node);
        map<int,int> fromn_hops; // map stores number of kMC hops from i-th node to noneliminated nodes
        for (vector<pair<Node*,Edge*>>::iterator it_nodevec=nodes_nbrs.begin();it_nodevec!=nodes_nbrs.end();++it_nodevec) {
            if (((*it_nodevec).first)->eliminated || (*it_nodevec).first==curr_node) continue;
            fromn_hops[((*it_nodevec).first)->node_id]=0;
        }
        /* sample transitions from eliminated to noneliminated nodes, not incl the i-th eliminated node, and also
           update transitions from eliminated nodes to the i-th node, except the self-loop of the i-th node.
           Note that only nodes directly connected to the i-th node are affected. */
        for (vector<pair<Node*,Edge*>>::iterator it_nodevec=nodes_nbrs.begin();it_nodevec!=nodes_nbrs.end();++it_nodevec) {
            if (!((*it_nodevec).first)->eliminated || (*it_nodevec).first==curr_node) continue;
            Edge *edgeptr=((*it_nodevec).first)->top_from;
            int hx=0; // number of transitions from eliminated node to the i-th node
            while (edgeptr!=nullptr) {
                if (!(edgeptr->deadts || edgeptr->to_node->eliminated)) {
                    double ratio=edgeptr->t/(edgeptr->t+edgeptr->dt);
                    int h_prev = edgeptr->h;
                    edgeptr->h = KPS::binomial_distribn(edgeptr->h,ratio,seed);
                    hx += h_prev-edgeptr->h;
                    fromn_hops[edgeptr->to_node->node_id] += hx;
                    if (debug) cout << "from: " << edgeptr->from_node->node_id << "  to: " << edgeptr->to_node->node_id \
                                    << "  R: " << ratio << "  old h: " << h_prev << "  new h: " << edgeptr->h << endl;
                }
                edgeptr=edgeptr->next_from;
            }
            ((*it_nodevec).second)->rev_edge->h = hx;
        }
        // update transitions from the i-th node to noneliminated nodes
        for (vector<pair<Node*,Edge*>>::iterator it_nodevec=nodes_nbrs.begin();it_nodevec!=nodes_nbrs.end();++it_nodevec) {
            if (((*it_nodevec).first)->eliminated || (*it_nodevec).first==curr_node) continue;
            ((*it_nodevec).second)->h += fromn_hops[((*it_nodevec).first)->node_id];
            if (debug) cout << "from: " << curr_node->node_id << "  to: " << ((*it_nodevec).first)->node_id \
                            << "  new h: " << ((*it_nodevec).second)->h << endl;
        }
        // reset flags for neighbouring nodes
        for (vector<pair<Node*,Edge*>>::iterator it_nodevec=nodes_nbrs.begin();it_nodevec!=nodes_nbrs.end();++it_nodevec) {
            ((*it_nodevec).first)->flag=false; }
        // sample the number of self-hops for the i-th node
        Edge *edgeptr = curr_node->top_from;
        while (edgeptr!=nullptr) {
            if (!(edgeptr->deadts || edgeptr->to_node->eliminated)) {
                curr_node->h += edgeptr->h; }
            edgeptr=edgeptr->next_from;
        }
        curr_node->h += KPS::negbinomial_distribn(curr_node->h,1.-curr_node->t,seed);
        if (debug) cout << "number of self-hops for node " << curr_node->node_id << ":  " << curr_node->h << endl;
    }
    // count the number of hops and stochastically draw the time for the escape trajectory
    double t_esc=0.; // sampled time for basin escape trajectory
    if (!ktn_kps->branchprobs) { // transiitons from all nodes are associated with the same waiting time
        int nhops=0; // total number of kMC hops
        for (const auto &node: ktn_kps->nodes) nhops += node.h;
        for (const auto &edge: ktn_kps->edges) {
            if (edge.deadts) continue; nhops += edge.h; }
        t_esc = KPS::gamma_distribn(nhops,tau,seed);
    } else { // the waiting times are different between nodes, and there are no self-loops
        for (const auto &node: ktn_kps->nodes) {
            int nhops=0;
            Edge *edgeptr = node.top_from;
            while (edgeptr!=nullptr) {
                if (!edgeptr->deadts) nhops += edgeptr->h;
                edgeptr = edgeptr->next_from;
            }
            t_esc += KPS::gamma_distribn(nhops,1./exp(node.k_esc),seed);
        }
    }
    if (debug) cout << "kps> finished iterative reverse randomisation" << endl;
    return t_esc;
}

/* Sample a node at the absorbing boundary of the current trapping basin, by the
   categorical sampling procedure based on T^(0) and T^(N) */
Node *KPS::sample_absorbing_node() {

    if (debug) cout << "\nkps> sample absorbing node, epsilon: " << epsilon->node_id << endl;
    int curr_comm_id = epsilon->comm_id;
    Node *next_node, *curr_node, *dummy_node;
    /* NB epsilon points to a node in the original network. At the start of each iteration of the following loop,
       curr_node points to a node in the transformed network. It is swapped for a node in the original subnetwork if
       it is a noneliminated node */
    curr_node = &ktn_kps->nodes[nodemap[epsilon->node_id-1]];
    do {
        if (debug) cout << "curr_node is: " << curr_node->node_id << endl;
        double rand_no = KMC_Standard_Methods::rand_unif_met(seed);
        double cum_t = 0.; // accumulated transition probability
        bool nonelimd = false; // flag indicates if the current node is transient noneliminated
        double factor = 0.;
        if (!curr_node->eliminated) {
            if (debug) cout << "  node has not been eliminated" << endl;
            dummy_node = &(*curr_node);
            curr_node = &ktn_kps_orig->nodes[curr_node->node_id-1]; // now points to a node in the untransformed subnetwork
            nonelimd = true;
        }
        // sample the next node using the appropriate probability distribution vector
        Edge *edgeptr = curr_node->top_from;
        if (nonelimd) factor=calc_gt_factor(curr_node);
        while (edgeptr!=nullptr) {
            if (edgeptr->deadts || edgeptr->to_node->eliminated) {
                edgeptr=edgeptr->next_from; continue; }
            cum_t += edgeptr->t;
            if (nonelimd) cum_t += (edgeptr->t)*(curr_node->t)/factor;
            if (debug) cout << "    to node: " << edgeptr->to_node->node_id << "  edgeptr->t: " << edgeptr->t \
                 << "  extra contribn: " << (edgeptr->t)*(curr_node->t)/factor << "  cum_t: " << cum_t << endl;
            if (cum_t>rand_no) { next_node = edgeptr->to_node; break; }
            edgeptr=edgeptr->next_from;
        }
        if (next_node==nullptr || cum_t-1.>1.E-08) {
            cout << "kps> GT error detected in sample_absorbing_node()" << endl; exit(EXIT_FAILURE); }
        // increment the number of kMC hops and set the new node
        if (nonelimd) {
            dummy_node->h++;
            curr_node = &ktn_kps->nodes[nodemap[next_node->node_id-1]];
        } else {
            edgeptr->h++;
            curr_node=next_node;
        }
        next_node=nullptr;
    } while (curr_node->comm_id==curr_comm_id);
    if (debug) cout << "after categorical sampling procedure the current node is: " << curr_node->node_id << endl;
    return curr_node;
}

/* Graph transformation to eliminate up to N nodes of the current trapping basin.
   Calculates the set of N-1 transition probability matrices {T^(n)} for 0 < n <= N.
   The transition network input to this function is the full network, and get_subnetwork() returns T^(0).
   The graph transformation is performed by performing a LU-decomposition of T^(0) */
void KPS::graph_transformation(const Network &ktn) {

    if (debug) cout << "\nkps> graph transformation" << endl;
    ktn_kps=get_subnetwork(ktn);
    ktn_kps_orig=get_subnetwork(ktn);
    ktn_l = new Network(N_B+N_c,0);
    ktn_u = new Network(N_B+N_c,0);
    ktn_l->edges.resize(8); ktn_u->edges.resize(5); // quack
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
    while (!gt_pq.empty() && N<nelim) {
        Node *node_elim=gt_pq.top();
        if (debug) cout << "N: " << N << " Node: " << node_elim->node_id << " priority: " << node_elim->udeg \
                        << " k: " << node_elim->top_from->k << endl;
        gt_pq.pop();
        node_elim = &ktn_kps->nodes[N]; // eliminate nodes in order of IDs
        gt_iteration(node_elim);
        basin_ids[node_elim->node_id-1]=1; // flag eliminated node
        eliminated_nodes.push_back(node_elim);
        N++;
        if (debug) { cout << "\nrunning debug tests on transformed network:" << endl; test_ktn(*ktn_kps); }
    }
    if (debug) {
        cout << "\nrunning debug tests on L:" << endl; test_ktn(*ktn_l);
        cout << "\nrunning debug tests on U:" << endl; test_ktn(*ktn_u);
    }
    if (N!=(!(N_B>nelim)?N_B:nelim)) {
        cout << "kps> fatal error: lost track of number of eliminated nodes" << endl; exit(EXIT_FAILURE); }
    if (debug) cout << "kps> finished graph transformation" << endl;
}

/* return the subnetwork corresponding to the active trapping basin and absorbing boundary nodes, to be transformed
   in the graph transformation phase of the kPS algorithm */
Network *KPS::get_subnetwork(const Network& ktn) {

    if (debug) cout << "\nkps> get_subnetwork: create TN of " << N_B+N_c << " nodes and " << N_e << " edges" << endl;
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
    if (debug) cout << "added " << n << " nodes and " << m << " edges to subnetwork" << endl;
    reset_kmc_hop_counts(*ktnptr);
    if (n!=N_B+N_c || m!=N_e) { cout << "kps> something went wrong in get_subnetwork()" << endl; exit(EXIT_FAILURE); }
//    Network *ktnptr = &const_cast<Network&>(ktn); // quack
    return ktnptr;
}

/* a single iteration of the graph transformation method. Argument is a pointer to the node to be
   eliminated from the network to which the ktn_kps pointer refers.
   The networks "L" and "U" required to undo the graph transformation iterations are updated */
void KPS::gt_iteration(Node *node_elim) {

    double factor=calc_gt_factor(node_elim); // equal to (1-T_{nn})
    if (debug) cout << "kps> eliminating node: " << node_elim->node_id << endl;
    // objects to queue all nbrs of the current elimd node, incl all elimd nbrs, and update relevant edges
    vector<Node*> nodes_nbrs;
    typedef struct {
        bool dirconn; // flag indicates if node is directly connected to current node being considered
        double t_fromn; // transition probability from eliminated node to this node
        double t_ton; // transition probability to eliminated node from this node
    } nbrnode;
    map<int,nbrnode> nbrnode_map; // map contains all nodes directly connected to the current elimd node, incl elimd nodes
    // update the self-loops of the L and U networks
    ktn_u->nodes[nodemap[node_elim->node_id-1]].t = -factor;
    ktn_l->nodes[nodemap[node_elim->node_id-1]].t = node_elim->t/factor;
    // update the weights for all edges from the elimd node to non-elimd nbr nodes, and self-loops of non-elimd nbr nodes
    Edge *edgeptr = node_elim->top_from;
    if (debug) cout << "updating edges from the eliminated node..." << endl;
    while (edgeptr!=nullptr) {
        if (edgeptr->deadts) { edgeptr=edgeptr->next_from; continue; }
        edgeptr->to_node->flag=true;
        nodes_nbrs.push_back(edgeptr->to_node); // queue nbr node
        nbrnode_map[edgeptr->to_node->node_id]=(nbrnode){false,edgeptr->t,edgeptr->rev_edge->t};
        if (debug) cout << "  to node: " << edgeptr->to_node->node_id << endl;
        // update L and U networks
        size_t pos;
/*
        (ktn_l->edges).emplace_back(Edge());
        (ktn_l->edges).back().t = edgeptr->rev_edge->t/factor;
        pos = (ktn_l->edges).size()-1;
        (ktn_l->edges).back().edge_pos = pos;
        (ktn_l->edges).back().from_node = &ktn_l->nodes[nodemap[edgeptr->to_node->node_id-1]];
        (ktn_l->edges).back().to_node = &ktn_l->nodes[nodemap[node_elim->node_id-1]];
        ktn_l->add_to_edge(nodemap[node_elim->node_id-1],pos);
        ktn_l->n_edges++;
*/
        ktn_l->edges[ktn_l->n_edges].t = edgeptr->rev_edge->t/factor;
        ktn_l->edges[ktn_l->n_edges].edge_pos = ktn_l->n_edges;
        ktn_l->edges[ktn_l->n_edges].from_node = &ktn_l->nodes[nodemap[edgeptr->to_node->node_id-1]];
        ktn_l->edges[ktn_l->n_edges].to_node = &ktn_l->nodes[nodemap[node_elim->node_id-1]];
        ktn_l->add_to_edge(nodemap[node_elim->node_id-1],ktn_l->n_edges);
        ktn_l->n_edges++;

        if (edgeptr->to_node->eliminated) { // do not update edges to elimd nodes and self-loops for elimd nodes
            edgeptr=edgeptr->next_from; continue; }
/*
        (ktn_u->edges).emplace_back(Edge());
        (ktn_u->edges).back().t = edgeptr->t;
        pos = (ktn_u->edges).size()-1;
        (ktn_u->edges).back().edge_pos = pos;
        (ktn_u->edges).back().from_node = &ktn_u->nodes[nodemap[node_elim->node_id-1]];
        (ktn_u->edges).back().to_node = &ktn_u->nodes[nodemap[edgeptr->to_node->node_id-1]];
        ktn_u->add_from_edge(nodemap[node_elim->node_id-1],pos);
        ktn_u->n_edges++;
*/
        ktn_u->edges[ktn_u->n_edges].t = edgeptr->t;
        ktn_u->edges[ktn_u->n_edges].edge_pos = ktn_u->n_edges;
        ktn_u->edges[ktn_u->n_edges].from_node = &ktn_u->nodes[nodemap[node_elim->node_id-1]];
        ktn_u->edges[ktn_u->n_edges].to_node = &ktn_u->nodes[nodemap[edgeptr->to_node->node_id-1]];
        ktn_u->add_from_edge(nodemap[node_elim->node_id-1],ktn_u->n_edges);
        ktn_u->n_edges++;

        // update subnetwork
        if (debug) cout << "    old node t: " << edgeptr->to_node->t << "  incr in node t: " \
                        << (edgeptr->t)*(edgeptr->rev_edge->t)/factor << endl;
        edgeptr->to_node->t += (edgeptr->t)*(edgeptr->rev_edge->t)/factor; // update self-loop of non-elimd nbr node
        if (debug) cout << "    old edge t: " << edgeptr->t << "  incr in t: " << (edgeptr->t)*(node_elim->t)/factor << endl;
        edgeptr->t += (edgeptr->t)*(node_elim->t)/factor; // update edge from elimd node to non-elimd nbr node
        edgeptr=edgeptr->next_from;
    }
    if (debug) cout << "updating edges between pairs of nodes both directly connected to the eliminated node..." << endl;
    // update the weights for all pairs of nodes directly connected to the eliminated node
    int old_n_edges = ktn_kps->n_edges; // number of edges in the network before we start adding edges in the GT algo
    for (vector<Node*>::iterator it_nodevec=nodes_nbrs.begin();it_nodevec!=nodes_nbrs.end();++it_nodevec) {
        if (debug) cout << "checking node: " << (*it_nodevec)->node_id << endl;
        edgeptr = (*it_nodevec)->top_from; // loop over edges to neighbouring nodes
        while (edgeptr!=nullptr) { // find pairs of nodes that are already directly connected to one another
            // skip nodes not directly connected to elimd node and edges to elimd nodes
            if (edgeptr->deadts || edgeptr->to_node->eliminated || !edgeptr->to_node->flag) {
                edgeptr=edgeptr->next_from; continue; }
            if (debug) cout << "  node " << (*it_nodevec)->node_id << " is directly connected to node " << edgeptr->to_node->node_id << endl;
            nbrnode_map[edgeptr->to_node->node_id].dirconn=true; // this pair of nodes are directly connected
            if (edgeptr->edge_pos>old_n_edges) { // skip nodes for which a new edge has already been added
                if (debug) cout << "    edge already added" << endl;
                edgeptr=edgeptr->next_from; continue; }
            if (debug) cout << "    old edge t: " << edgeptr->t << "  incr in t: " \
                            << (nbrnode_map[edgeptr->from_node->node_id].t_ton)*\
                               (nbrnode_map[edgeptr->to_node->node_id].t_fromn)/factor << endl;
            edgeptr->t += (nbrnode_map[edgeptr->from_node->node_id].t_ton)*\
                (nbrnode_map[edgeptr->to_node->node_id].t_fromn)/factor;
            edgeptr=edgeptr->next_from;
        }
        if (debug) cout << "  checking for neighbours of eliminated node that are not already connected to this node" << endl;
        for (map<int,nbrnode>::iterator it_map=nbrnode_map.begin();it_map!=nbrnode_map.end();++it_map) {
            if (debug) cout << "    checking nbr node: "<< it_map->first << endl;
            if (it_map->first==(*it_nodevec)->node_id) continue; // self-loops of neighbouring nodes already accounted for
            if (ktn_kps->nodes[nodemap[it_map->first-1]].eliminated) continue;
            if ((it_map->second).dirconn) { (it_map->second).dirconn=false; continue; } // reset flag
            if (debug) { cout << "    node " << (*it_nodevec)->node_id << " is not directly connected to node " << it_map->first << endl;
                cout << "    t of new edge: " << (it_map->second).t_fromn*nbrnode_map[(*it_nodevec)->node_id].t_ton/factor << endl; }
            // nodes are directly connected to the elimd node but not to one another, add an edge in the transformed network
            (ktn_kps->edges).push_back(Edge());
            (ktn_kps->edges).back().t = (it_map->second).t_fromn*nbrnode_map[(*it_nodevec)->node_id].t_ton/factor;
            size_t pos = (ktn_kps->edges).size()-1;
            (ktn_kps->edges).back().edge_pos = pos;
            (ktn_kps->edges).back().label = node_elim->node_id;
//            Network::add_edge_network(ktn_kps,&ktn_kps->nodes[nodemap[edgeptr->from_node->node_id-1]], \
                &ktn_kps->nodes[nodemap[edgeptr->to_node->node_id-1]],pos);
            (ktn_kps->edges).back().from_node = &ktn_kps->nodes[nodemap[(*it_nodevec)->node_id-1]];
            (ktn_kps->edges).back().to_node = &ktn_kps->nodes[nodemap[it_map->first-1]];
            ktn_kps->add_from_edge(nodemap[(*it_nodevec)->node_id-1],pos);
            ktn_kps->add_to_edge(nodemap[it_map->first-1],pos);
            ktn_kps->n_edges++; // increment number of edges in the network
            // reverse edge
            (ktn_kps->edges).push_back(Edge());
            if ((*it_nodevec)->eliminated) {
                (ktn_kps->edges).back().t = 0.; // dummy value
            } else {
                if (debug) cout << "    t of new reverse edge: " << (it_map->second).t_ton*nbrnode_map[(*it_nodevec)->node_id].t_fromn/factor << endl;
                (ktn_kps->edges).back().t = (it_map->second).t_ton*nbrnode_map[(*it_nodevec)->node_id].t_fromn/factor;
            }
            (ktn_kps->edges).back().edge_pos = pos+1;
            (ktn_kps->edges).back().label = node_elim->node_id;
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
        if (!edgeptr->deadts) edgeptr->to_node->flag=false;
        edgeptr = edgeptr->next_from;
    }
    node_elim->eliminated=true; // this flag negates the need to zero the weights to the eliminated node
}

/* undo a single iteration of the graph transformation. Argument is a pointer to the node to be un-eliminated from the network */
vector<pair<Node*,Edge*>> KPS::undo_gt_iteration(Node *node_elim) {

    if (debug) cout << "kps> undoing elimination of node " << node_elim->node_id << endl;
    if (!node_elim->eliminated) throw exception(); // node is already noneliminated
    node_elim->eliminated=false;
    // set the self-loop for the restored node
    node_elim->t = -(ktn_l->nodes[nodemap[node_elim->node_id-1]].t)*(ktn_u->nodes[nodemap[node_elim->node_id-1]].t);
    // construct list of elimd+nonelimd nodes neighbouring the restored node, along with corresponding edges from the restored node
    vector<pair<Node*,Edge*>> nodes_nbrs;
    Edge *edgeptr = node_elim->top_from;
    while (edgeptr!=nullptr) {
        if (!edgeptr->deadts) {
            nodes_nbrs.push_back(make_pair(edgeptr->to_node,edgeptr));
            edgeptr->to_node->flag=true;
        }
        edgeptr=edgeptr->next_from;
    }
    if (debug) {
        cout << "list of neighbouring nodes:" << endl;
        for (auto &neptr: nodes_nbrs) cout << "  " << (neptr.first)->node_id;
        cout << endl; }
    // update the remaining edges for pairs of nodes connected to the restored node 
    edgeptr = ktn_l->nodes[nodemap[node_elim->node_id-1]].top_to;
    if (debug) cout << "doing edges FROM neighbouring nodes" << endl;
    while (edgeptr!=nullptr) {
        Edge *edgeptr2 = ktn_kps->nodes[nodemap[edgeptr->from_node->node_id-1]].top_from;
        if (!edgeptr2->from_node->eliminated) { // quack but what if edge is dead?
            if (debug) cout << " neighbour node " << edgeptr2->from_node->node_id << " is noneliminated, relevant L elem: " << edgeptr->t << endl;
            edgeptr2->from_node->dt = edgeptr->t;
            if (debug) cout << " new t of node is: " << edgeptr2->from_node->dt << endl;
        }
        while (edgeptr2!=nullptr) {
            if (debug) cout << "  edge from: " << edgeptr2->from_node->node_id << "  to: " << edgeptr2->to_node->node_id << endl;
            if (edgeptr2->label==node_elim->node_id) edgeptr2->deadts=true;
            if (edgeptr2->deadts) { edgeptr2=edgeptr2->next_from; continue; }
            if (edgeptr2->to_node->flag) {
                if (debug) cout << "    to node is flagged, relevant L elem: " << edgeptr->t << endl;
                edgeptr2->dt = edgeptr->t;
//            } else if (edgeptr2->to_node==node_elim) {
//                cout << "    to node is eliminated node, relevant U elem: " << ktn_u->nodes[nodemap[node_elim->node_id-1]].t << endl;
//                edgeptr2->dt = ktn_u->nodes[nodemap[node_elim->node_id-1]].t;
            }
            edgeptr2 = edgeptr2->next_from;            
        }
        edgeptr=edgeptr->next_to;
    }
    edgeptr = ktn_u->nodes[nodemap[node_elim->node_id-1]].top_from;
    if (debug) cout << "doing edges TO neighbouring nodes" << endl;
    while (edgeptr!=nullptr) {
        Edge *edgeptr2 = ktn_kps->nodes[nodemap[edgeptr->to_node->node_id-1]].top_to;
        if (!edgeptr2->to_node->eliminated) { // quack but what if edge is dead?
            if (debug) cout << " neighbour node: " << edgeptr2->to_node->node_id << " is noneliminated, relevant U elem: " << edgeptr->t << endl;
            edgeptr2->to_node->dt *= edgeptr->t;
            edgeptr2->to_node->t -= edgeptr2->to_node->dt;
            if (debug) cout << " new t of node is: " << edgeptr2->to_node->t << endl;
        }
        while (edgeptr2!=nullptr) {
            if (debug) cout << "  edge from: " << edgeptr2->from_node->node_id << "  to: " << edgeptr2->to_node->node_id << endl;
            if (edgeptr2->label==node_elim->node_id) edgeptr2->deadts=true;
            if (edgeptr2->deadts) {edgeptr2=edgeptr2->next_to; continue; }
            if (edgeptr2->from_node->flag) {
                if (debug) cout << "    from node is flagged, relevant U elem: " << edgeptr->t << endl;
                edgeptr2->dt *= edgeptr->t;
                edgeptr2->t -= edgeptr2->dt;
                if (debug) cout << "      new t of edge is: " << edgeptr2->t << endl;
            } else if (edgeptr2->from_node==node_elim) {
                if (debug) cout << "    from node is eliminated node, relevant L elem: " << ktn_l->nodes[nodemap[node_elim->node_id-1]].t \
                                << "  relevant U elem: " << edgeptr->t << endl;
//                edgeptr2->dt *= ktn_l->nodes[nodemap[node_elim->node_id-1]].t;
//                edgeptr2->t -= edgeptr2->dt;
                edgeptr2->t -= (ktn_l->nodes[nodemap[node_elim->node_id-1]].t)*(edgeptr->t);
                if (debug) cout << "      new t of edge is: " << edgeptr2->t << endl;
            }
            edgeptr2 = edgeptr2->next_to;
        }
        edgeptr=edgeptr->next_from;
    }
    if (debug) test_ktn(*ktn_kps);
    return nodes_nbrs;
}

/* calculate the factor (1-T_{nn}) needed in the elimination of the n-th node in graph transformation */
double KPS::calc_gt_factor(Node *node_elim) {

    double factor=0.; // equal to (1-T_{nn})
    if (node_elim->t>0.999) { // loop over neighbouring edges to maintain numerical precision
        Edge *edgeptr = node_elim->top_from;
        while (edgeptr!=nullptr) {
            if (!(edgeptr->deadts || edgeptr->to_node->eliminated)) factor += edgeptr->t;
            edgeptr=edgeptr->next_from;
        }
    } else { factor=1.-node_elim->t; }
    return factor;
}

/* reset the counts for the number of kMC hops, which appear in the Node and Edge data structures, for
   the Network object passed as a pointer */
void KPS::reset_kmc_hop_counts(Network &ktn) {
    for (auto &node: ktn.nodes) node.h=0;
    for (auto &edge: ktn.edges) edge.h=0;
}

/* Update path quantities along a trajectory, where the (unordered) path is specified by the kMC hop counts
   ("h") in the Node and Edge members of the subnetwork pointed to by ktn_kps */
void KPS::update_path_quantities(double t_esc, const Node *curr_node) {

    if (debug) cout << "kps> updating path quantities" << endl;
    if (ktn_kps==nullptr) throw exception();
    walker.curr_node = &(*curr_node);
    walker.t += t_esc;
    for (const auto &node: ktn_kps->nodes) {
        walker.k += node.h;
        if (!ktn_kps->branchprobs) walker.p += node.h*log(node.t);
    }
    for (const auto &edge: ktn_kps->edges) {
        if (edge.deadts) continue;
        walker.k += edge.h;
        walker.p += edge.h*log(edge.t);
        if (ktn_kps->branchprobs) { // can also calculate the contribution to the entropy flow
            walker.s += (edge.h)*(edge.rev_edge->k-edge.k); }
    }
}

/* Gamma distribution with shape parameter a and rate parameter 1./b */
double KPS::gamma_distribn(int a, double b, int seed) {

    static default_random_engine generator(seed);
    gamma_distribution<double> gamma_distrib(a,b);
    return gamma_distrib(generator);
}

/* Binomial distribution with trial number h and success probability p.
   Returns the number of successes after h Bernoulli trials. */
int KPS::binomial_distribn(int h, double p, int seed) {

    static default_random_engine generator(seed);
    if (!((h>=0) || ((p>=0.) && (p<=1.)))) throw exception();
    if (h==0)  { return 0;
    } else if (p==1.) { return h; }
    binomial_distribution<int> binom_distrib(h,p);
    return binom_distrib(generator);
}

/* Negative binomial distribution with success number r and success probability p.
   Returns the number of failures before the r-th success. */
int KPS::negbinomial_distribn(int r, double p, int seed) {

    static default_random_engine generator(seed);
    if (!((r>=0) || ((p>0.) && (p<=1.)))) throw exception();
    if ((r==0) || (p==1.)) return 0;
    negative_binomial_distribution<int> neg_binom_distrib(r,p);
    return neg_binom_distrib(generator);
}

/* Exponential distribution with rate parameter 1./tau */
double KPS::exp_distribn(double tau, int seed) {

    static default_random_engine generator(seed);
    exponential_distribution<double> exp_distrib(1./tau);
    return exp_distrib(generator);
}
