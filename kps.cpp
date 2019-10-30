/*
File containing functions relating to kinetic path sampling
*/

#include "kmc_methods.h"

using namespace std;

KPS::KPS(Network& ktn, int nelim, double tau) {

    N_max=nelim;
    basin_ids.reserve(ktn.n_nodes);
    ktn.get_tmtx_lin(tau);
}

KPS::~KPS() {}

/* Reset data of previous kPS iteration and find the microstates of the current trapping basin */
void KPS::setup_basin_sets(Network &ktn) {

    N_c=0; N=0; alpha=epsilon; epsilon=nullptr;
    for (int i=0;i<ktn.n_nodes;i++) basin_ids[i]=0;
}

/* Iterative reverse randomisation procedure to stochastically sample the hopping matrix
   H^(0) corresponding to T^(0), given H^(N) and the {T^(n)} for 0 <= n <= N.
   Return a sampled time for the stochastic escape trajectory. */
double KPS::iterative_reverse_randomisation() {

    for (int i=N;i<=0;i--) {

    }
    double t_esc = KPS::gamma_distribn(0,0.); // time for escape trajectory
    return t_esc;
}

/* Sample a node at the absorbing boundary of the current trapping basin, by the
   categorical sampling procedure based on T^(0) and T^(N) */
Node * KPS::sample_absorbing_node(Network &ktn) {

    return nullptr;
}

/* Graph transformation to eliminate up to N nodes of the current trapping basin.
   Calculates the set of N-1 transition probability matrices {T^(n)} for 0 < n <= N.
   The transition network input to this function is T^(0) */
void KPS::graph_transformation(Network &ktn) {

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
