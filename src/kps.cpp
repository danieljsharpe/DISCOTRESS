/*
File containing functions relating to kinetic path sampling (kPS). See:
M. Athenes and V. V. Bulatov, Phys. Rev. Lett. 113, 230601 (2014).
M. Athenes, S. Kaur, G. Adjanor, T. Vanacker and T. Jourdan, Phys. Rev. Materials 3, 103802 (2019).
D. J. Wales, J. Chem. Phys. 130, 204111 (2009).

The functions in this routine are also used as the basis for the state reduction algorithms to compute MFPTs, mean path lengths, committor and
absorption probabilities, stationary probabilities, the average numbers of node visits and the node visitation probabilities. See:
W. K. Grassmann, M. I. Taksar, and D. P. Heyman, Oper. Res. 33, 1107-1116 (1985).
T. J. Sheskin, Oper. Res. 33, 228-235 (1985).
J. Kohlas, Zeit. Oper. Res. 30, 197-207 (1986).
D. P. Heyman, SIAM J. Matrix Anal. Appl. 16, 954-963 (1995).
D. P. Heyman and D. P. O'Leary, SIAM J. Matrix Anal. Appl. 19, 534-540 (1998).
E. Seneta, SIAM. J. Matrix Anal. Appl. 19, 556-563 (1998).
I. Sonin, Adv. Math. 145, 159-188 (1999).
I. Sonin and J. Thornton, SIAM J. Matrix Anal. Appl. 23, 209-224 (2001).
T. Dayar and N. Akar, SIAM J. Matrix. Anal. Appl. 27, 396-412 (2005).
J. D. Stevenson and D. J. Wales, J. Chem. Phys. 141, 041104 (2014).
J. J. Hunter, Spec. Matrices 4, 151-175 (2016).
J. J. Hunter, Linear Algebra Appl. 549, 100-122 (2018).
R. S. MacKay and J. D. Robinson, Phil. Trans. Roy. Soc. A 376, 20170232 (2018).

This file is a part of DISCOTRESS, a software package to simulate the dynamics on arbitrary continuous- and discrete-time Markov chains (CTMCs).
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
#include "statereduction.h"
#include <queue>
#include <cmath>
#include <random>
#include <numeric>
#include <iostream>

using namespace std;

/* constructor for KPS class */
KPS::KPS(const Network &ktn, int nelim, int kpskmcsteps, bool adaptivecomms, double adaptminrate, \
         const Traj_args &traj_args) : Traj_Method(traj_args) {

    cout << "kps> kPS parameters:\n  max. no. of eliminated nodes: " \
         << nelim << "\n  no. of basins: " << ktn.ncomms << " \tno. of kMC steps after kPS iteration: " << kpskmcsteps \
         << "\n  adaptive definition of communities (y/n): " << adaptivecomms \
         << "\tmin. allowed rate in adaptive communities: " << adaptminrate << endl;
    this->nelim=nelim; this->kpskmcsteps=kpskmcsteps;
    this->adaptivecomms=adaptivecomms; this->adaptminrate=adaptminrate;
    basin_ids.resize(ktn.n_nodes);
}

/* destructor for KPS class */
KPS::~KPS() {
    if (ktn_kps!=nullptr) delete ktn_kps; if (ktn_kps_orig!=nullptr) delete ktn_kps_orig;
    if (ktn_kps_gt!=nullptr) delete ktn_kps_gt;
    if (ktn_l!=nullptr) delete ktn_l; if (ktn_u!=nullptr) delete ktn_u;
    if (sr_args.mfpt) mfpt_vals.clear();
}

/* copy constructor for KPS class */
KPS::KPS(const KPS &kps_obj) : Traj_Method(kps_obj) {
    this->nelim=kps_obj.nelim; this->kpskmcsteps=kps_obj.kpskmcsteps;
    this->adaptivecomms=false; this->adaptminrate=-1.;
    if (kps_obj.statereduction) this->set_statereduction_procs(kps_obj.sr_args);
    this->basin_ids.resize(kps_obj.basin_ids.size());
}

/* call to this function indicates that the purpose fo the computation is state reduction to calculate exact dynamical quantities, and not
   a dynamical simulation; set the state reduction procedures that are to be performed */
void KPS::set_statereduction_procs(const SR_args &sr_args) {
    cout << "kps> one or more state reduction procedures are specified" << endl;
    this->sr_args.absorption=sr_args.absorption; this->sr_args.committor=sr_args.committor;
    this->sr_args.fundamentalirred=sr_args.fundamentalirred; this->sr_args.fundamentalred=sr_args.fundamentalred;
    this->sr_args.gth=sr_args.gth; this->sr_args.mfpt=sr_args.mfpt;
}

void KPS::test_ktn(const Network &ktn) {
    cout << "debug> ktn info: no. of nodes: " << ktn.n_nodes << " no. of edges: " << ktn.n_edges << endl;
    for (int i=0;i<ktn.n_nodes;i++) {
        cout << "node: " << ktn.nodes[i].node_id << endl;
        if (!ktn.nodes[i].eliminated) cout << "  to: " << ktn.nodes[i].node_id << "  t: " << ktn.nodes[i].t << "  h: " << ktn.nodes[i].h << endl;
        Edge *edgeptr = ktn.nodes[i].top_from;
        while (edgeptr!=nullptr) {
            if (!edgeptr->deadts && !edgeptr->to_node->eliminated) {
                cout << "  to: " << edgeptr->to_node->node_id << "  t: " << edgeptr->t \
                     << "  h: " << edgeptr->h << endl; }
            edgeptr = edgeptr->next_from;
        }
    }
    cout << "\n" << endl;
}

/* perform a single kPS basin escape iteration */
void KPS::kmc_iteration(const Network &ktn, Walker &walker) {

    if (!(!adaptivecomms && ktn.ncomms==2 && ktn_kps_orig!=nullptr)) { // for a two-state problem, only need to setup basin and do GT once
        setup_basin_sets(ktn,walker,true);
        graph_transformation(ktn);
    } else {
        setup_basin_sets(ktn,walker,false); // get the new initial node without updating the definition of the basin
    }
    if (statereduction && !sr_args.fundamentalirred && !sr_args.mfpt && !sr_args.gth) {
        return;
    } else if (!statereduction) {
        Node *dummy_alpha = sample_absorbing_node();
        alpha = &ktn.nodes[dummy_alpha->node_id-1];
    }
    long double t_traj = iterative_reverse_randomisation();
    if (statereduction) {
        if (sr_args.mfpt) calc_mfpt();
        if (sr_args.gth) calc_gth();
        return;
    }
    update_path_quantities(walker,t_traj,alpha);
    delete ktn_kps; ktn_kps=nullptr;
    if (!(!adaptivecomms && ktn.ncomms==2)) {
        delete ktn_kps_orig; ktn_kps_orig=nullptr;
        delete ktn_l; delete ktn_u; ktn_l=nullptr; ktn_u=nullptr;
    } else { // restore the graph transformed subnetwork
        ktn_kps = new Network(*ktn_kps_gt);
    }
    epsilon=alpha; alpha=nullptr;
}

/* perform the specified number of kMC iterations, to be executed after a basin escape. The idea is to attempt
   to move away from the transition boundary region of a communtiy before simulating another basin escape.
   Optional argument dt specifies a maximum time for the walker before the loop is forced to break (default value
   infinity), used when the dimensionality reduction wrapper method simulates trajectories of fixed length */
void KPS::do_bkl_steps(const Network &ktn, Walker &walker, long double maxtime) {

    if (adaptivecomms) return;
    int n_kmcit=0;
    while ((n_kmcit<kpskmcsteps || ktn.comm_sizes[epsilon->comm_id]>nelim) && walker.t<maxtime) { // quack force BKL simulation to continue if active community is large
        BKL::bkl(walker,discretetime,ktn.accumprobs,seed);
        alpha=walker.curr_node;
        if (ktn.nbins>0 && !ktn.nodesB.empty()) walker.visited[alpha->bin_id]=true;
        if (alpha->comm_id!=epsilon->comm_id || walker.t>maxtime) { // traj data is not dumped unless comm changes, regardless of tintvl, except if (DIMREDN) max time is exceeded
            this->dump_traj(walker,walker.curr_node->aorb==-1,false,maxtime); }
        epsilon=alpha;
        if (alpha->aorb==-1 || alpha->aorb==1) return; // note that the BKL iterations are terminated if the simulation returns to B
        n_kmcit++;
    }
}

void KPS::reset_nodeptrs() {
    epsilon=nullptr; alpha=nullptr;
}

/* Reset data of previous kPS iteration and find the microstates of the current trapping basin */
void KPS::setup_basin_sets(const Network &ktn, Walker &walker, bool get_new_basin) {

    if (debug) cout << "\nkps> setting up basin sets" << endl;
    bool newpath=false;
    if (statereduction) { // not simulating a trajectory, set epsilon to any node not in A
        for (const Node &node: ktn.nodes) {
            if (node.aorb==-1) continue;
            epsilon=&node; break;
        }
    } else if (!epsilon) { // first iteration of A<-B path, need to set starting node
        epsilon = Wrapper_Method::get_initial_node(ktn,walker,seed);
        if (tintvl>=0.) walker.dump_walker_info(true,0.,walker.curr_node,dumpintvls);
        next_tintvl=tintvl;
    }
    if (!get_new_basin) return; // the basin is not to be updated
    N_c=0; N=0; N_B=0; N_e=0;
    fill(basin_ids.begin(),basin_ids.end(),0); // reset basin IDs (zero flag indicates absorbing nonboundary node)
    if (!adaptivecomms) { // basin IDs are based on community IDs
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
                if (edgeptr->deadts) { edgeptr=edgeptr->next_from; continue; }
                if (edgeptr->to_node->comm_id!=epsilon->comm_id && !basin_ids[edgeptr->to_node->node_id-1]) {
                    basin_ids[edgeptr->to_node->node_id-1]=3; // flag absorbing boundary node
                    N_c++;
                    if (debug) cout << "  " << edgeptr->to_node->node_id;
                }
                if (basin_ids[edgeptr->to_node->node_id-1]==3) N_e++;
                edgeptr=edgeptr->next_from;
            }
        }
        if (debug) cout << endl;
    } else {
        vector<int> nodes_in_comm = Wrapper_Method::find_comm_onthefly(ktn,epsilon,adaptminrate,nelim);
        basin_ids=nodes_in_comm;
        for (int i=0;i<ktn.n_nodes;i++) {
            if (basin_ids[i]==2) {
                N_B++; N_e+=ktn.nodes[i].udeg;
                
            } else if (basin_ids[i]==3) {
                N_c++; int N_e_add=0;
                Edge *edgeptr = ktn.nodes[i].top_from;
                while (edgeptr!=nullptr) {
                    if (!edgeptr->deadts && basin_ids[edgeptr->to_node->node_id-1]==2) N_e_add++;
                    edgeptr=edgeptr->next_from;
                }
                N_e+=N_e_add;
            }
        }
    }
    eliminated_nodes.clear(); nodemap.clear();
    eliminated_nodes.reserve(!(N_B>nelim)?N_B:nelim);
    if (debug) {
        cout << "\nthread no.: " << omp_get_thread_num() << endl;
        cout << "number of eliminated nodes: " << (!(N_B>nelim)?N_B:nelim) << endl;
        cout << "number of nodes in basin: " << N_B << " number of absorbing boundary nodes: " << N_c << endl;
        cout << "number of edges of subnetwork: " << N_e << endl;
        cout << "epsilon: " << epsilon->node_id << endl;
    }
}

/* Iterative reverse randomisation procedure to stochastically sample the hopping matrix
   H^(0) corresponding to T^(0), given H^(N) and the {T^(n)} for 0 <= n <= N.
   Return a sampled time for the stochastic escape trajectory. */
long double KPS::iterative_reverse_randomisation() {

    if (debug) {
        cout << "\nkps> iterative reverse randomisation" << endl;
        cout << "N is: " << N << endl; if (!statereduction) cout << "node alpha: " << alpha->node_id << endl; }
    // main loop of the iterative reverse randomisation procedure
    for (int i=N;i>0;i--) {
        Node *curr_node = &(ktn_kps->nodes[nodemap[eliminated_nodes[i-1]]-1]);
        vector<pair<Node*,Edge*>> nodes_nbrs = undo_gt_iteration(curr_node);
        // reset flags for neighbouring nodes
        for (vector<pair<Node*,Edge*>>::iterator it_nodevec=nodes_nbrs.begin();it_nodevec!=nodes_nbrs.end();++it_nodevec) {
            ((*it_nodevec).first)->flag=false; }
        if (statereduction) continue;
//        cout << "  i: " << i << "    undone GT elimination of node: " << curr_node->node_id << endl;
        // vector stores number of kMC hops from i-th node to noneliminated nodes, other elems are irrelevant
        vector<unsigned long long int> fromn_hops(N_B+N_c);
        for (vector<pair<Node*,Edge*>>::iterator it_nodevec=nodes_nbrs.begin();it_nodevec!=nodes_nbrs.end();++it_nodevec) {
            if (((*it_nodevec).first)->eliminated || (*it_nodevec).first==curr_node) continue;
            fromn_hops[((*it_nodevec).first)->node_pos]=0;
        }
        /* sample transitions from eliminated to noneliminated nodes, not incl the i-th eliminated node, and also
           update transitions from eliminated nodes to the i-th node, except the self-loop of the i-th node.
           Note that only nodes directly connected to the i-th node are affected. */
        for (vector<pair<Node*,Edge*>>::iterator it_nodevec=nodes_nbrs.begin();it_nodevec!=nodes_nbrs.end();++it_nodevec) {
            if (basin_ids[((*it_nodevec).first)->node_id-1]!=1 || (*it_nodevec).first==curr_node) continue;
            unsigned long long int hx=0; // number of transitions from eliminated node to the i-th node
            Edge *edgeptr=((*it_nodevec).first)->top_from;
            if (debug) cout << "from node: " << edgeptr->from_node->node_id << endl;
            // update the self-loop for this node
//            cout << "    stage 1" << endl;
            if (!edgeptr->from_node->eliminated) {
                long double ratio=edgeptr->from_node->t/(edgeptr->from_node->t+edgeptr->from_node->dt);
                unsigned long long int h_prev = edgeptr->from_node->h;
//                cout << "      about to draw from B distribn. h: " << edgeptr->from_node->h << "  ratio: " << ratio << endl;
                edgeptr->from_node->h = KPS::binomial_distribn(edgeptr->from_node->h,ratio,seed);
                hx += h_prev-edgeptr->from_node->h;
                fromn_hops[edgeptr->from_node->node_pos] += h_prev-edgeptr->from_node->h;
                if (debug) cout << " old node h: " << h_prev << "  new node h: " << edgeptr->from_node->h \
                                << "  R: " << ratio << endl;
            }
            edgeptr->from_node->dt=0.L;
            // update edges
//            cout << "    stage 2" << endl;
            while (edgeptr!=nullptr) {
                if (edgeptr->to_node->eliminated || (edgeptr->deadts && edgeptr->label!=curr_node->node_id) \
                    || edgeptr->to_node==curr_node) {
                    edgeptr->dt=0.L; edgeptr=edgeptr->next_from; continue;
                }
                long double ratio;
                if (!edgeptr->deadts) { ratio=edgeptr->t/(edgeptr->t+edgeptr->dt);
                } else { ratio=0.L; }
                unsigned long long int h_prev = edgeptr->h;
//                cout << "      about to draw from B distribn. h: " << edgeptr->from_node->h << "  ratio: " << ratio << endl;                
                edgeptr->h = KPS::binomial_distribn(edgeptr->h,ratio,seed);
                hx += h_prev-edgeptr->h;
                fromn_hops[edgeptr->to_node->node_pos] += h_prev-edgeptr->h;
                if (debug) cout << "  to node : " << edgeptr->to_node->node_id \
                                << "  R: " << ratio << "  old h: " << h_prev << "  new h: " << edgeptr->h << endl;
                edgeptr->dt=0.L; edgeptr=edgeptr->next_from;
            }
            ((*it_nodevec).second)->rev_edge->h = hx; // transitions from eliminated nodes to the i-th node
            if (debug) cout << "  new h to elimd node: " << hx << endl;
        }
//        cout << "    stage 3" << endl;
        // update transitions from the i-th node to noneliminated nodes
        for (vector<pair<Node*,Edge*>>::iterator it_nodevec=nodes_nbrs.begin();it_nodevec!=nodes_nbrs.end();++it_nodevec) {
            if (((*it_nodevec).first)->eliminated || (*it_nodevec).first==curr_node) continue;
            ((*it_nodevec).second)->h += fromn_hops[((*it_nodevec).first)->node_pos];
            if (debug) cout << "from elimd node: " << curr_node->node_id << "  to: " << ((*it_nodevec).first)->node_id \
                            << "  new h: " << ((*it_nodevec).second)->h << endl;
        }
        // sample the number of self-hops for the i-th node
        unsigned long long int nhops=0; // number of kMC hops from the i-th node to alternative nonelimd nodes (ie no self-loops)
        Edge *edgeptr = curr_node->top_from;
        while (edgeptr!=nullptr) {
            if (!(edgeptr->deadts || edgeptr->to_node->eliminated)) {
                nhops += edgeptr->h; }
            edgeptr=edgeptr->next_from;
        }
        long double nb_prob = Network::calc_gt_factor(*curr_node);
//        cout << "    about to draw from NB distribn. nhops: " << nhops << " nb_prob: " << nb_prob << endl;
        curr_node->h = KPS::negbinomial_distribn(nhops,nb_prob,seed);
//        cout << "    curr_node->h is now: " << curr_node->h << endl;
        if (debug) {
            cout << "tot no of hops from node " << curr_node->node_id << " to alt nonelimd nodes: " \
                 << nhops << "  1-t: " << nb_prob << endl;
            cout << "number of self-hops for node " << curr_node->node_id << ":  " << curr_node->h << endl;
            cout << "network after restoring node " << curr_node->node_id << endl; test_ktn(*ktn_kps);
        }
    }
    // count the number of hops and sample the time associated with the escape trajectory
    long double t_traj=0.L; // sampled time for basin escape trajectory
    for (const auto &node: ktn_kps->nodes) {
        unsigned long long int nhops=0;
        nhops += node.h;
        const Edge *edgeptr = node.top_from;
        while (edgeptr!=nullptr) {
            if (!edgeptr->deadts) nhops += edgeptr->h;
            edgeptr = edgeptr->next_from;
        }
        if (discretetime) { t_traj += static_cast<long double>(nhops)*node.t_esc;
        } else { t_traj += KPS::gamma_distribn(nhops,node.t_esc,seed); }
    }
    if (debug) {
        cout << "network after iterative reverse randomisation:" << endl; test_ktn(*ktn_kps);
        cout << "kps> finished iterative reverse randomisation" << endl; }
    return t_traj;
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
    curr_node = &ktn_kps->nodes[nodemap[epsilon->node_id]-1];
    do {
        if (debug) cout << "curr_node is: " << curr_node->node_id << endl;
        double rand_no = Wrapper_Method::rand_unif_met(seed);
        long double cum_t = 0.L; // accumulated transition probability
        bool nonelimd = false; // flag indicates if the current node is transient noneliminated
        long double factor = 0.L;
        if (!curr_node->eliminated) {
            if (debug) cout << "  node has not been eliminated" << endl;
            dummy_node = &(*curr_node);
            curr_node = &ktn_kps_orig->nodes[curr_node->node_id-1]; // now points to a node in the untransformed subnetwork
            nonelimd = true;
        }
        // sample the next node using the appropriate probability distribution vector
        Edge *edgeptr = curr_node->top_from;
        if (nonelimd) factor = Network::calc_gt_factor(*curr_node);
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
            curr_node = &ktn_kps->nodes[nodemap[next_node->node_id]-1];
        } else {
            edgeptr->h++;
            curr_node=next_node;
        }
        next_node=nullptr;
        if (adaptivecomms && basin_ids[curr_node->node_id-1]==3) break; // reached absorbing boundary of on-the-fly community
    } while (curr_node->comm_id==curr_comm_id);
    if (debug) cout << "after categorical sampling procedure the current node is: " << curr_node->node_id << endl;
    return curr_node;
}

/* Graph transformation to eliminate up to N nodes of the current trapping basin.
   Calculates the set of N-1 transition probability matrices {T^(n)} for 0 < n <= N.
   The Markovian network input to this function is the full network, and get_subnetwork() returns T^(0).
   The graph transformation is achieved by performing a LU-decomposition of T^(0) */
void KPS::graph_transformation(const Network &ktn) {

    if (debug) cout << "\nkps> graph transformation" << endl;
    ktn_kps=get_subnetwork(ktn,true);
    ktn_kps->ncomms=ktn.ncomms;
    /* the original, L and U network are not needed for certain state reduction computations, which only require a forward pass phase of GT */
    if (!statereduction || sr_args.fundamentalirred || sr_args.mfpt || sr_args.gth) {
    ktn_kps_orig=get_subnetwork(ktn,false);
    ktn_l = new Network(N_B+N_c,0);
    ktn_u = new Network(N_B+N_c,0);
    vector<int> noderange(N_B); iota(noderange.begin(),noderange.end(),N_c);
    ktn_l->edges.resize((N_B*(N_B+N_c))-N_B);
    ktn_u->edges.resize(accumulate(noderange.begin(),noderange.end(),0));
    for (int i=0;i<ktn_kps->n_nodes;i++) {
        ktn_l->nodes[i] = ktn_kps->nodes[i];
        ktn_u->nodes[i] = ktn_kps->nodes[i];
        ktn_l->nodes[i].node_pos=i; ktn_u->nodes[i].node_pos=i;
        ktn_l->nodes[i].t=0.L; ktn_u->nodes[i].t=0.L; // the "transn probs" in the L and U TNs are the values to "undo" GT
    }
    if (sr_args.mfpt) { mfpt_vals.resize(ktn_kps->n_nodes); fill(mfpt_vals.begin(),mfpt_vals.end(),0.L); }
    }
    /* comparison function for the priority queue. Note that computation of the committor probabilities within the state reduction
       procedure takes place when only nodes of the set A and B remain, so elimination of nodes not in B should be prioritised */
    auto cmp = [committor_val=sr_args.committor](Node *l,Node *r) {
        if (committor_val && l->aorb==1 && r->aorb!=1) { return true;
        } else if (committor_val && l->aorb!=1 && r->aorb==1) { return false; }
        return l->udeg >= r->udeg;
    };
    priority_queue<Node*,vector<Node*>,decltype(cmp)> gt_pq(cmp); // priority queue of nodes (based on out-degree)
    for (vector<Node>::iterator it_nodevec=ktn_kps->nodes.begin();it_nodevec!=ktn_kps->nodes.end();++it_nodevec) {
        if (sr_args.fundamentalred && !it_nodevec->flag) continue; // only eliminate dummy nodes when computing the absorbing fundamental matrix
        if ((!adaptivecomms && it_nodevec->comm_id!=epsilon->comm_id) || \
            (adaptivecomms && basin_ids[it_nodevec->node_id-1]!=2)) continue;
        gt_pq.push(&(*it_nodevec));
    }
    bool done_committor=false;
    while (!gt_pq.empty() && N<nelim) {
        Node *node_elim=gt_pq.top();
//        node_elim = &ktn_kps->nodes[N]; // quack eliminate nodes in order of IDs
        gt_pq.pop();
        if (sr_args.committor && !done_committor && node_elim->aorb==1) { // only nodes not in A and B remain at this point; compute committor probabilities
            calc_committor(ktn); done_committor=true;
        }
        gt_iteration(node_elim);
        basin_ids[node_elim->node_id-1]=1; // flag eliminated node
        eliminated_nodes.push_back(node_elim->node_id);
        N++;
        if (debug) { cout << "\nrunning debug tests on transformed network:" << endl; test_ktn(*ktn_kps); }
        if (sr_args.gth && gt_pq.empty()) { // if GTH, only [the single node in] A remains at this point;
            for (vector<Node>::iterator it_nodevec=ktn_kps->nodes.begin();it_nodevec!=ktn_kps->nodes.end();++it_nodevec) {
                if (!it_nodevec->eliminated) {
                    cout << "unnormalised pi of node: " << it_nodevec->node_id << " is set to 1." << endl;
                    cout << "self-loop for this node: " << it_nodevec->t << endl;
                    it_nodevec->pi=1.L; // note that the GTH stationary probabilities are initially not stored as logs
                    mu=1.L; break;
                }
            }
        }
    }
    if (!adaptivecomms && ktn.ncomms==2 && ktn_kps_gt==nullptr) ktn_kps_gt = new Network(*ktn_kps);
    if (N!=(!(N_B>nelim)?N_B:nelim)) {
        cout << "kps> fatal error: lost track of number of eliminated nodes" << endl; exit(EXIT_FAILURE); }
    if (debug) cout << "kps> finished graph transformation" << endl;
    if (statereduction) rewrite_stat_probs(ktn); // the stationary probs of initial nodes in the ktn_kps object are rewritten to be the initial probs
    if (sr_args.absorption) calc_absprobs(); // only nodes not in A remain at this point; compute absorption probabilities
    if (sr_args.fundamentalred) calc_fundamentalred(ktn); // the remaining edges are the elements of the fundamental matrix for a reducible Markov chain
}

/* return the subnetwork corresponding to the active trapping basin and absorbing boundary nodes, to be transformed
   in the graph transformation phase of the kPS algorithm */
Network *KPS::get_subnetwork(const Network& ktn, bool resize_edgevec) {

    if (debug) cout << "\nkps> get_subnetwork: create TN of " << N_B+N_c << " nodes and " << N_e << " edges" << endl;
    Network *ktnptr = new Network(N_B+N_c,N_e);
    if (resize_edgevec) ktnptr->edges.resize((N_B*(N_B-1))+(2*N_B*N_c));
    ktnptr->branchprobs=ktn.branchprobs;
    int j=0;
    for (int i=0;i<ktn.n_nodes;i++) {
        if (!basin_ids[i]) continue;
        nodemap[i+1]=j+1;
        ktnptr->nodes[j] = ktn.nodes[i];
        ktnptr->nodes[j].node_pos=j; j++;
    }
    int m=0, n=0;
    vector<bool> edgemask(2*ktn.n_edges,false);
    // note that the indices of the edge vector in the subnetwork are not in a meaningful order
    for (auto &node: ktnptr->nodes) {
        n++;
        const Node *node_orig = &ktn.nodes[node.node_id-1];
        // for absorbing node, do not incl any FROM edges, or any TO edges for non-basin nbr nodes, in the subnetwork
        if ((!adaptivecomms && node_orig->comm_id!=epsilon->comm_id) || \
            (adaptivecomms && basin_ids[node_orig->node_id-1]!=2)) continue;
        const Edge *edgeptr = node_orig->top_from;
        while (edgeptr!=nullptr) {
            if (edgeptr->deadts || edgemask[edgeptr->edge_id]) { edgeptr=edgeptr->next_from; continue; }
            ktnptr->edges[m] = *edgeptr; // edge of subnetwork inherits properties (transn rate etc) of node in full network
            ktnptr->edges[m].edge_id = m;
            ktnptr->edges[m].from_node = &ktnptr->nodes[nodemap[edgeptr->from_node->node_id]-1];
            ktnptr->edges[m].to_node = &ktnptr->nodes[nodemap[edgeptr->to_node->node_id]-1];
            ktnptr->add_from_edge(nodemap[edgeptr->from_node->node_id]-1,m);
            ktnptr->add_to_edge(nodemap[edgeptr->to_node->node_id]-1,m);
//            Network::add_edge_network(ktnptr,ktnptr->nodes[nodemap[edgeptr->from_node->node_id]]-1, \
                ktnptr->nodes[nodemap[edgeptr->to_node->node_id]-1],m);
            m++; edgemask[edgeptr->edge_id]=true;
            const Edge *edgeptr_rev = edgeptr->rev_edge;
            if (edgeptr_rev->deadts || edgemask[edgeptr_rev->edge_id]) {
                edgeptr=edgeptr->next_from; continue; }
            // reverse edge
            ktnptr->edges[m] = *edgeptr_rev;
            ktnptr->edges[m].edge_id = m;
            ktnptr->edges[m].from_node = &ktnptr->nodes[nodemap[edgeptr_rev->from_node->node_id]-1];
            ktnptr->edges[m].to_node = &ktnptr->nodes[nodemap[edgeptr_rev->to_node->node_id]-1];
            ktnptr->add_from_edge(nodemap[edgeptr_rev->from_node->node_id]-1,m);
            ktnptr->add_to_edge(nodemap[edgeptr_rev->to_node->node_id]-1,m);
//            Network::add_edge_network(ktnptr,ktnptr->nodes[nodemap[edgeptr_rev->from_node->node_id]-1], \
                ktnptr->nodes[nodemap[edgeptr_rev->to_node->node_id]-1],m);
            ktnptr->edges[m-1].rev_edge = &ktnptr->edges[m];
            ktnptr->edges[m].rev_edge = &ktnptr->edges[m-1];
            m++; edgemask[edgeptr_rev->edge_id]=true;
        }
    }
    if (debug) cout << "added " << n << " nodes and " << m << " edges to subnetwork" << endl;
    if (n!=N_B+N_c || m!=N_e) {
        cout << "kps> something went wrong in get_subnetwork(). Nodes: " << n << " edges: " << m << endl; exit(EXIT_FAILURE); }
    return ktnptr;
}

/* a single iteration of the graph transformation method. Argument is a pointer to the node to be
   eliminated from the network to which the ktn_kps pointer refers.
   The networks "L" and "U" required to undo the graph transformation iterations are updated */
void KPS::gt_iteration(Node *node_elim) {

    long double factor = Network::calc_gt_factor(*node_elim); // equal to (1-T_{nn})
    if (debug) cout << "kps> eliminating node: " << node_elim->node_id << endl;
    // objects to queue all nbrs of the current elimd node, incl all elimd nbrs, and update relevant edges
    vector<Node*> nodes_nbrs;
    typedef struct {
        bool dirconn; // flag indicates if node is directly connected to current node being considered
        long double t_fromn; // transition probability from eliminated node to this node
        long double t_ton; // transition probability to eliminated node from this node
    } nbrnode;
    // vector of which relevant entries are for all nodes directly connected to the current elimd node, incl elimd nodes
    vector<nbrnode> nbrnode_vec(N_B+N_c,(nbrnode){false,0.L,0.L});
    // update the self-loops of the L and U networks
    if (!statereduction || sr_args.fundamentalirred || sr_args.mfpt || sr_args.gth) {
    ktn_u->nodes[node_elim->node_pos].t = -factor;
    ktn_l->nodes[node_elim->node_pos].t = node_elim->t/factor;
    }
    // update the weights for all edges from the elimd node to non-elimd nbr nodes, and self-loops of non-elimd nbr nodes
    Edge *edgeptr = node_elim->top_from;
    if (debug) cout << "updating edges from the eliminated node..." << endl;
    while (edgeptr!=nullptr) {
        if (edgeptr->deadts) { edgeptr=edgeptr->next_from; continue; }
        if (debug) cout << "  to node: " << edgeptr->to_node->node_id << endl;
        edgeptr->to_node->flag=true;
        nodes_nbrs.push_back(edgeptr->to_node); // queue nbr node
        nbrnode_vec[edgeptr->to_node->node_pos].t_fromn=edgeptr->t;
        nbrnode_vec[edgeptr->to_node->node_pos].t_ton=edgeptr->rev_edge->t;
        if (!statereduction || sr_args.fundamentalirred || sr_args.mfpt || sr_args.gth) {
        // update L and U networks
        ktn_l->edges[ktn_l->n_edges].t = edgeptr->rev_edge->t/factor;
        ktn_l->edges[ktn_l->n_edges].edge_id = ktn_l->n_edges;
        ktn_l->edges[ktn_l->n_edges].from_node = &ktn_l->nodes[edgeptr->to_node->node_pos];
        ktn_l->edges[ktn_l->n_edges].to_node = &ktn_l->nodes[node_elim->node_pos];
        ktn_l->add_to_edge(node_elim->node_pos,ktn_l->n_edges);
        ktn_l->n_edges++;
        if (edgeptr->to_node->eliminated) { // do not update edges to elimd nodes and self-loops for elimd nodes
            edgeptr=edgeptr->next_from; continue; }
        ktn_u->edges[ktn_u->n_edges].t = edgeptr->t;
        ktn_u->edges[ktn_u->n_edges].edge_id = ktn_u->n_edges;
        ktn_u->edges[ktn_u->n_edges].from_node = &ktn_u->nodes[node_elim->node_pos];
        ktn_u->edges[ktn_u->n_edges].to_node = &ktn_u->nodes[edgeptr->to_node->node_pos];
        ktn_u->add_from_edge(node_elim->node_pos,ktn_u->n_edges);
        ktn_u->n_edges++;
        }
        // renormalise mean waiting time for the neighbouring node (when noneliminated) if the computation is to compute exact MFPTs
        if (sr_args.mfpt && !edgeptr->to_node->eliminated && edgeptr->to_node->aorb!=-1) {
            edgeptr->to_node->t_esc += (edgeptr->rev_edge->t)*(node_elim->t_esc)/factor; }
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
        bool node1_abs = (basin_ids[(*it_nodevec)->node_id-1]==3);
        edgeptr = (*it_nodevec)->top_from; // loop over edges to neighbouring nodes
        while (edgeptr!=nullptr) { // find pairs of nodes that are already directly connected to one another
            // skip nodes not directly connected to elimd node and edges to elimd nodes
            if (edgeptr->deadts || edgeptr->to_node->eliminated || !edgeptr->to_node->flag || \
                (node1_abs && basin_ids[edgeptr->to_node->node_id-1]==3)) {
                edgeptr=edgeptr->next_from; continue; }
            if (debug) cout << "  node " << (*it_nodevec)->node_id << " is directly connected to node " \
                            << edgeptr->to_node->node_id << endl;
            nbrnode_vec[edgeptr->to_node->node_pos].dirconn=true; // this pair of nodes are directly connected
            if (edgeptr->edge_id>old_n_edges) { // skip nodes for which a new edge has already been added
                if (debug) cout << "    edge already added" << endl;
                edgeptr=edgeptr->next_from; continue; }
            if (debug) cout << "    old edge t: " << edgeptr->t << "  incr in t: " \
                            << (nbrnode_vec[edgeptr->from_node->node_pos].t_ton)*\
                               (nbrnode_vec[edgeptr->to_node->node_pos].t_fromn)/factor << endl;
            edgeptr->t += (nbrnode_vec[edgeptr->from_node->node_pos].t_ton)*\
                (nbrnode_vec[edgeptr->to_node->node_pos].t_fromn)/factor;
            edgeptr=edgeptr->next_from;
        }
        if (debug) cout << "  checking for nbrs of elimd node that are not already connected to this node" << endl;
        for (vector<Node*>::iterator it_nodevec2=nodes_nbrs.begin();it_nodevec2!=nodes_nbrs.end();++it_nodevec2) {
            /* skip self-loops of neighbour nodes (already accounted for), proposed edges TO eliminated nodes (accounted
               for when the reverse direction is found), and proposed edges connecting pairs of absorbing nodes (irrelevant) */
            if (debug) cout << "    checking nbr node: "<< (*it_nodevec2)->node_id << endl;
            if ((*it_nodevec2)==(*it_nodevec) || (*it_nodevec2)->eliminated || \
                (node1_abs && basin_ids[(*it_nodevec2)->node_id-1]==3)) continue;
            int node1_pos=(*it_nodevec)->node_pos, node2_pos=(*it_nodevec2)->node_pos;
            if (nbrnode_vec[node2_pos].dirconn) { nbrnode_vec[node2_pos].dirconn=false; continue; } // reset flag
            if (debug) {
                cout << "    node " << (*it_nodevec)->node_id << " is not directly connected to node " \
                     << (*it_nodevec2)->node_id << "\n    t of new edge: " \
                     << nbrnode_vec[node2_pos].t_fromn*nbrnode_vec[node1_pos].t_ton/factor << endl; }
            // nodes are directly connected to the elimd node but not to one another, add an edge in the transformed network
            ktn_kps->edges[ktn_kps->n_edges].t = nbrnode_vec[node2_pos].t_fromn*nbrnode_vec[node1_pos].t_ton/factor;
            ktn_kps->edges[ktn_kps->n_edges].edge_id = ktn_kps->n_edges;
            ktn_kps->edges[ktn_kps->n_edges].label = node_elim->node_id;
            ktn_kps->edges[ktn_kps->n_edges].from_node = &ktn_kps->nodes[node1_pos];
            ktn_kps->edges[ktn_kps->n_edges].to_node = &ktn_kps->nodes[node2_pos];
            ktn_kps->add_from_edge(node1_pos,ktn_kps->n_edges);
            ktn_kps->add_to_edge(node2_pos,ktn_kps->n_edges);
            ktn_kps->n_edges++;
            // reverse edge
            if ((*it_nodevec)->eliminated) {
                ktn_kps->edges[ktn_kps->n_edges].t = 0.L; // dummy value
            } else {
                if (debug) cout << "    t of new reverse edge: " \
                                << nbrnode_vec[node2_pos].t_ton*nbrnode_vec[node1_pos].t_fromn/factor << endl;
                ktn_kps->edges[ktn_kps->n_edges].t = nbrnode_vec[node2_pos].t_ton*nbrnode_vec[node1_pos].t_fromn/factor;
            }
            ktn_kps->edges[ktn_kps->n_edges].edge_id = ktn_kps->n_edges;
            ktn_kps->edges[ktn_kps->n_edges].label = node_elim->node_id;
            ktn_kps->edges[ktn_kps->n_edges].from_node = &ktn_kps->nodes[node2_pos];
            ktn_kps->edges[ktn_kps->n_edges].to_node = &ktn_kps->nodes[node1_pos];
            ktn_kps->add_from_edge(node2_pos,ktn_kps->n_edges);
            ktn_kps->add_to_edge(node1_pos,ktn_kps->n_edges);

            ktn_kps->edges[ktn_kps->n_edges-1].rev_edge = &ktn_kps->edges[ktn_kps->n_edges];
            ktn_kps->edges[ktn_kps->n_edges].rev_edge = &ktn_kps->edges[ktn_kps->n_edges-1];
            ktn_kps->n_edges++;
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

/* undo a single iteration of the graph transformation.
   Argument is a pointer to the node to be un-eliminated from the network, and which exists in the Network object
   pointed to by ktn_kps */
vector<pair<Node*,Edge*>> KPS::undo_gt_iteration(Node *node_elim) {

    if (debug) cout << "\nkps> undoing elimination of node " << node_elim->node_id << endl;
    if (!node_elim->eliminated) throw exception(); // node is already noneliminated
    node_elim->eliminated=false;
    // set the self-loop for the restored node
    node_elim->t = -(ktn_l->nodes[node_elim->node_pos].t)*(ktn_u->nodes[node_elim->node_pos].t);
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
    edgeptr = ktn_l->nodes[node_elim->node_pos].top_to;
    if (debug) cout << "doing edges FROM neighbouring nodes" << endl;
    while (edgeptr!=nullptr) {
        Edge *edgeptr2 = ktn_kps->nodes[edgeptr->from_node->node_pos].top_from;
        if (!edgeptr2->from_node->eliminated) { // quack but what if edge is dead?
            if (debug) cout << " neighbour node " << edgeptr2->from_node->node_id \
                            << " is noneliminated, relevant L elem: " << edgeptr->t << endl;
            edgeptr2->from_node->dt = edgeptr->t;
            if (debug) cout << " new t of node is: " << edgeptr2->from_node->dt << endl;
        }
        while (edgeptr2!=nullptr) {
            if (debug) cout << "  edge from: " << edgeptr2->from_node->node_id \
                            << "  to: " << edgeptr2->to_node->node_id << endl;
            if (edgeptr2->label==node_elim->node_id) edgeptr2->deadts=true;
            if (edgeptr2->deadts) { edgeptr2=edgeptr2->next_from; continue; }
            if (edgeptr2->to_node->flag) {
                if (debug) cout << "    to node is flagged, relevant L elem: " << edgeptr->t << endl;
                edgeptr2->dt = edgeptr->t;
//            } else if (edgeptr2->to_node==node_elim) {
//                cout << "    to node is eliminated node, relevant U elem: " \
                       << ktn_u->nodes[node_elim->node_pos].t << endl;
//                edgeptr2->dt = ktn_u->nodes[node_elim->node_pos].t;
            }
            edgeptr2 = edgeptr2->next_from;            
        }
        edgeptr=edgeptr->next_to;
    }
    edgeptr = ktn_u->nodes[node_elim->node_pos].top_from;
    if (debug) cout << "doing edges TO neighbouring nodes" << endl;
    while (edgeptr!=nullptr) {
        Edge *edgeptr2 = ktn_kps->nodes[edgeptr->to_node->node_pos].top_to;
        if (!edgeptr2->to_node->eliminated) { // quack but what if edge is dead?
            if (debug) cout << " neighbour node: " << edgeptr2->to_node->node_id \
                            << " is noneliminated, relevant U elem: " << edgeptr->t << endl;
            edgeptr2->to_node->dt *= edgeptr->t;
            edgeptr2->to_node->t -= edgeptr2->to_node->dt;
            if (debug) cout << " new t of node is: " << edgeptr2->to_node->t << endl;
        }
        while (edgeptr2!=nullptr) {
            if (debug) cout << "  edge from: " << edgeptr2->from_node->node_id \
                            << "  to: " << edgeptr2->to_node->node_id << endl;
            if (edgeptr2->label==node_elim->node_id) edgeptr2->deadts=true;
            if (edgeptr2->deadts) {edgeptr2=edgeptr2->next_to; continue; }
            if (edgeptr2->from_node->flag) {
                if (debug) cout << "    from node is flagged, relevant U elem: " << edgeptr->t << endl;
                edgeptr2->dt *= edgeptr->t;
                edgeptr2->t -= edgeptr2->dt;
                if (debug) cout << "      new t of edge is: " << edgeptr2->t << endl;
            } else if (edgeptr2->from_node==node_elim) {
                if (debug) cout << "    from node is eliminated node, relevant L elem: " \
                                << ktn_l->nodes[node_elim->node_pos].t \
                                << "  relevant U elem: " << edgeptr->t << endl;
//                edgeptr2->dt *= ktn_l->nodes[node_elim->node_pos]].t;
//                edgeptr2->t -= edgeptr2->dt;
                edgeptr2->t -= (ktn_l->nodes[node_elim->node_pos].t)*edgeptr->t;
                if (debug) cout << "      new t of edge is: " << edgeptr2->t << endl;
            }
            edgeptr2 = edgeptr2->next_to;
        }
        edgeptr=edgeptr->next_from;
    }
    if (sr_args.mfpt) {
        mfpt_vals[node_elim->node_pos] = node_elim->t_esc;
        Edge *edgeptr = node_elim->top_from;
        while (edgeptr!=nullptr) {
            if (!edgeptr->deadts && !edgeptr->to_node->eliminated && edgeptr->to_node->aorb!=-1) {
                mfpt_vals[node_elim->node_pos] += edgeptr->t*mfpt_vals[edgeptr->to_node->node_pos]; }
            edgeptr=edgeptr->next_from;
        }
        long double factor = Network::calc_gt_factor(*node_elim);
        mfpt_vals[node_elim->node_pos] *= 1.L/factor;
    }
    if (sr_args.gth) {
        cout << "\nrestored node: " << node_elim->node_id << endl;
        cout << "  self-loop: " << node_elim->t << endl;
        long double new_pi=0.L;
        Edge *edgeptr = node_elim->top_to;
        while (edgeptr!=nullptr) {
            if (!edgeptr->deadts && !edgeptr->from_node->eliminated) {
                cout << "  edge from: " << edgeptr->from_node->node_id << "    pi: " << edgeptr->from_node->pi << "   t: " << edgeptr->t << endl;
                new_pi += edgeptr->from_node->pi*edgeptr->t; }
            edgeptr=edgeptr->next_to;
        }
        cout << "    new_pi is: " << new_pi << endl;
        node_elim->pi = new_pi; mu += new_pi;
    }
    while (edgeptr!=nullptr) {
        if (!edgeptr->deadts && !edgeptr->to_node->eliminated) cout << "    to: " << edgeptr->to_node->node_id << "    t: "<< edgeptr->t << endl;
        edgeptr=edgeptr->next_from;
    }
    return nodes_nbrs;
}

/* Update path quantities along a trajectory, where the (unordered) path is specified by the kMC hop counts
   ("h") in the Node and Edge members of the subnetwork pointed to by ktn_kps.
   Transition probabilities associated with nodes and edges should not be accumulated values (this feature
   should only be set for use with pure BKL simulations) */
void KPS::update_path_quantities(Walker &walker, long double t_traj, const Node *curr_node) {

    if (debug) cout << "kps> updating path quantities" << endl;
    if (ktn_kps==nullptr) throw exception();
    walker.prev_node = walker.curr_node;
    walker.curr_node = &(*curr_node);
    walker.t += t_traj;
    for (const auto &node: ktn_kps->nodes) {
        if (!ktn_kps->branchprobs && node.h>0) {
            walker.k += node.h;
            walker.p += -1.L*static_cast<long double>(node.h)*log(node.t);
            // no need to update entropy flow along paths because contribution from self-loop transitions is zero
            if (ktn_kps->ncomms>0 && !walker.visited.empty()) walker.visited[node.bin_id]=true;
        }
        Edge *edgeptr = node.top_from;
        while (edgeptr!=nullptr) {
            if (edgeptr->deadts || edgeptr->h==0) { edgeptr=edgeptr->next_from; continue; }
            walker.k += edgeptr->h;
            walker.p += -1.L*static_cast<long double>(edgeptr->h)*log(edgeptr->t);
            if (ktn_kps->ncomms>0 && !walker.visited.empty()) walker.visited[edgeptr->to_node->bin_id]=true;
            if (!discretetime) {
                walker.s += static_cast<long double>(edgeptr->h)*(edgeptr->rev_edge->k-edgeptr->k);
            } else {
                walker.s += static_cast<long double>(edgeptr->h)*log(edgeptr->rev_edge->t/edgeptr->t);
            }
            edgeptr=edgeptr->next_from;
        }
    }
}

/* Gamma distribution with shape parameter a and rate parameter 1./b */
long double KPS::gamma_distribn(unsigned long long int a, long double b, int seed) {

    static default_random_engine generator(seed);
    gamma_distribution<long double> gamma_distrib(a,b);
    return gamma_distrib(generator);
}

/* Binomial distribution with trial number h and success probability p.
   Returns the number of successes after h Bernoulli trials. */
unsigned long long int KPS::binomial_distribn(unsigned long long int h, long double p, int seed) {

    static default_random_engine generator(seed);
    if (h<0 || (p>1. && h>0) ) { // || (p<0. && h>0)) {
cout << "h: " << h << " p: " << p << endl; throw exception(); } // quack
    if (h==0 || p==0.)  { return 0;
    } else if (p==1.) { return h; }
    binomial_distribution<unsigned long long int> binom_distrib(h,p);
    return binom_distrib(generator);
}

/* Negative binomial distribution with success number r and success probability p.
   Returns the number of failures before the r-th success. */
unsigned long long int KPS::negbinomial_distribn(unsigned long long int r, long double p, int seed) {

    static default_random_engine generator(seed);
    if (!(r>=0 && (p>0. && p<=1.)) && !(r==0 &p==0.)) { cout << "r: " << r << " p: " << p << endl; throw exception(); }
    if (r==0) return 0;
    negative_binomial_distribution<unsigned long long int> neg_binom_distrib(r,p);
    return neg_binom_distrib(generator);
}

/* Exponential distribution with rate parameter 1./tau */
long double KPS::exp_distribn(long double tau, int seed) {

    static default_random_engine generator(seed);
    exponential_distribution<long double> exp_distrib(1.L/tau);
    return exp_distrib(generator);
}
