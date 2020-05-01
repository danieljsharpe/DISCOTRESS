/*
File containing functions to perform the simulation when no enhanced sampling "wrapper" method is used to maintain the ensemble
of trajectories

This file is a part of DISCOTRESS, a software package to simulate the dynamics on arbitrary continuous time Markov chains (CTMCs).
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

using namespace std;

/* main loop to drive a standard kMC simulation (no enhanced sampling wrapper method) */
void STD_KMC::run_enhanced_kmc(const Network &ktn, Traj_Method &traj_method_obj) {

    cout << "\nstd_kmc> beginning kMC simulation with no enhanced sampling wrapper method" << endl;
/*
    long double dummy_randno = Wrapper_Method::rand_unif_met(seed); // seed generator
    Node *dummy_node = Wrapper_Method::get_initial_node(ktn,walker,seed);
    walker.dump_walker_info(false,true,tintvl>=0.);
    n_ab=0; n_traj=0; int n_kmcit=1;
    double next_tintvl=tintvl; // next time interval for writing trajectory data
    bool leftb=false; // flag indicates if trajectory has left initial set of nodes yet
    while ((n_ab<maxn_abpaths) && (n_kmcit<=maxit)) { // algo terminates after max no of iterations of the standard kMC algorithm
        (*kmc_func)(walker);
        if (ktn.nbins>0) walker.visited[walker.curr_node->bin_id]=true;
        if (!leftb && walker.curr_node->aorb!=1) leftb=true;
        walker.dump_walker_info(walker.curr_node->aorb==-1,false, \
            (walker.curr_node->aorb==-1 || tintvl==0. || (tintvl>0. && walker.t>=next_tintvl)));
        if (tintvl>0. && walker.t>=next_tintvl) { // reached time interval for dumping trajectory data, calc next interval
            while (walker.t>=next_tintvl) next_tintvl+=tintvl; }
        n_kmcit++;
        if (walker.curr_node->aorb==-1 || (walker.curr_node->aorb==1 && leftb)) {
            // traj has reached absorbing macrostate A or has returned to B
            update_tp_stats(walker,walker.curr_node->aorb==-1,ktn.nbins>0);
            if (walker.curr_node->aorb==1) continue;
            walker.reset_walker_info();
            walker.path_no++;
            dummy_node = Wrapper_Method::get_initial_node(ktn,walker,seed);
            if ((n_ab<maxn_abpaths) && (n_kmcit<maxit) && (tintvl>=0.)) walker.dump_walker_info(false,true);
            leftb=false; next_tintvl=tintvl;
        }
    }
    cout << "std_kmc> standard kMC simulation terminated after " << n_kmcit-1 << " kMC iterations. Simulated " \
         << n_ab << " transition paths" << endl;
    if (ktn.nbins>0) calc_tp_stats(ktn.nbins);
*/

    n_ab=0; n_traj=0; int n_it=1;
    while ((n_ab<maxn_abpaths) && (n_it<=maxit)) { // if using kPS or MCAMC, algo terminates when max no of basin escapes have been simulated
        bool donebklsteps=false;
        traj_method_obj.kmc_iteration(ktn,walker);
        traj_method_obj.dump_traj(walker,walker.curr_node->aorb==-1,false);
        n_it++;
        check_if_endpoint:
            if (walker.curr_node->aorb==-1 || walker.curr_node->aorb==1) { // traj has reached absorbing macrostate A or has returned to B
                update_tp_stats(walker,walker.curr_node->aorb==-1,!adaptivecomms);
                if (walker.curr_node->aorb==-1) { // transition path, reset walker
                    walker.reset_walker_info();
                    walker.path_no++;
                    traj_method_obj.reset_nodeptrs();
                    continue;
                } else if (ktn.ncomms>0) {
                    walker.visited[walker.curr_node->bin_id]=true;
                }
            }
            if (donebklsteps) continue;
        traj_method_obj.do_bkl_steps(ktn,walker);
        donebklsteps=true;
        goto check_if_endpoint;
    }
    cout << "\nstd_kmc> simulation terminated after " << n_it-1 << " iterations. Simulated " \
         << n_ab << " transition paths" << endl;
    if (!adaptivecomms) calc_tp_stats(ktn.nbins); // calc committors and transn path densities for communities and write to file
}
