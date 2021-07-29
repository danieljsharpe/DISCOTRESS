# DISCOTRESS user information

## Compilation and execution

Begin by downloading the latest version of the DISCOTRESS repository, which can be done from the command line:
```bash
git clone https://github.com/danieljsharpe/DISCOTRESS
```

Then navigate to the directory containing the source code with `cd DISCOTRESS/src` and compile using:
```bash
g++ -std=c++17 discotress.cpp kmc_methods.cpp we.cpp ffs.cpp neus.cpp milestoning.cpp rea.cpp kps.cpp mcamc.cpp keywords.cpp network.cpp -o discotress -fopenmp
```

To run the program, simply type the magic word: `discotress`, having provided the necessary input files documented below.

DISCOTRESS is tested using v5.4.0 of the gcc compiler, which supports OpenMP v4.0. These versions are therefore recommended but not required.

Get started with the [tutorials](https://github.com/danieljsharpe/DISCOTRESS_tutorials).

## Input files

The instructions to the program are set by the keywords in the _input.kmc_ file. 

Additional mandatory input files specify the network topology and parameters:

File | Format | Description
---- | ------ | -----------
*stat\_prob.dat* | single col, length **NNODES** | natural log of stationary probabilities for nodes
*edge\_conns.dat* | double col, length **NEDGES** | each line specifies a bidirectional edge connecting nodes with the stated IDs (indexed from 1)
*edge\_weights.dat* | double col, length **NEDGES** | transition probabilities (for DTMC, if **DISCRETETIME**) or natural log of transition rates (for CTMC, default) associated with forward and reverse transitions (consistent with *edge\_conns.dat*)
*nodes.A* | single col, length set by **NODESAFILE** | list of IDs for nodes in the absorbing set &#120068; (indexed from 1)
*nodes.B* | single col, length set by **NODESBFILE** | list of IDs for nodes in the initial set &#120069; (indexed from 1)

Additional optional input files, corresponding to particular optional keywords, relate to the simulation parameters and to the recording of some statistics:

File | Format | Keyword | Description
---- | ------ | ------- | -----------
*communities.dat* | single col, length **NNODES** | **COMMSFILE** | community IDs for nodes (indexed from 0), specifies network partitioning for use in enhanced sampling methods
*bins.dat* | single col, length **NNODES** | **BINSFILE** | bin IDs for nodes (indexed from 0), used for calculating transition path statistics (defaults to *communities.dat*)
*initcond.dat* | single col, length cf. **NODESBFILE** | **INITCONDFILE** | initial occupation probability distribution for nodes of the initial set &#120069;. Defaults to be proportional to a local equilibrium within &#120069;
*ntrajs.dat* | single col, length cf. **COMMSFILE** | **DIMREDUCTION** | numbers of short trajectories to run, initialised from each community

## Output files

The following output files are written as the result of a computation to simulate trajectories. For a description of the output files when instead instructing the program to perform an exact computation, see the documentation for the keywords in the *input.kmc* file that relate to state reduction procedures.

File | Description | Format (columns)
---- | ----------- | ----------------
*fpp\_properties.dat* | properties of simulated &#120068; &#8592; &#120069; paths, together yielding numerical estimates of the probability distributions for path properties in the first passage path ensemble | path no. / path time / path length / ln of path probability / path entropy flow
*tp\_stats.dat* | bin statistics for the &#120068; &#8592; &#120069; transition path ensemble, written if communities were specified | bin ID / no. of reactive (direct &#120068; &#8592; &#120069;) paths for which bin is visited / no. of paths for which bin is visited and trajectory returned to initial set &#120069; / reactive visitation probability / committor probability
*walker.x.y.dat* | trajectory information dumped at the specified time intervals (or when a trajectory escapes from a community, depending on options). *x* is the walker ID, *y* is the path number | node ID / community ID / path time / path length / path action (negative ln of path probability) / path entropy flow

## Main keywords

The following is a list of keywords that must always appear in the _input.kmc_ file.

----

**NNODES** `int`  
  mandatory, number of nodes in the Markov chain.

**NEDGES** `int`   
  mandatory, number of (bidirectional) edges in the Markov chain.

**WRAPPER** `str`  
  mandatory, employ a 'wrapping' enhanced sampling strategy for accelerating the observation of rare events. Further detail on wrapper methods is given below. Options:  
-    **BTOA**    \- simulate paths initialised in state &#120069; and terminating at state &#120068;, no enhanced sampling strategy employed
-    **FIXEDT**  \- simulate paths of fixed total time, no enhanced sampling strategy employed
-    **DIMREDN** \- special class to simulate many short nonequilibrium trajectories starting from each community in turn, used for estimation of a coarse-grained Markov chain. See **DIMREDUCTION** keyword for more detail 
-    **WE**      \- weighted ensemble sampling
-    **FFS**     \- forward flux sampling
-    **NEUS**    \- non-equilibrium umbrella sampling
-    **MILES**   \- milestoning
-    **REA**     \- recursive enumeration algorithm, a special wrapper method to calculate the highest-probability paths

**TRAJ** `str`  
  mandatory, method for propagating individual trajectories. Options:  
-    **BKL**     \- Bortz-Kalos-Lebowitz \(aka n-fold way\) algorithm, i.e. standard kinetic Monte Carlo (kMC)
-    **KPS**     \- kinetic path sampling algorithm
-    **MCAMC**   \- Monte Carlo with absorbing Markov chains algorithm  

**NODESAFILE** `str` `int`  
  mandatory if not **WRAPPER DIMREDN**, name of the file containing the node IDs (indexed from 1) belonging to the &#120068; (absorbing) set, and number of nodes in the &#120068; set.

**NODESBFILE** `str` `int`  
  mandatory if not **WRAPPER DIMREDN**, name of the file containing the node IDs (indexed from 1) belonging to the &#120069; (initial) set, and number of nodes in the &#120069; set.

----

## Additional detail on wrapper method options

The **WRAPPER** method handles a set of walkers (independent trajectories) that are each propagated by the chosen **TRAJ** method. The following provides further detail on the requirements and considerations relating to each **WRAPPER** method option.

----

**BTOA**  
  the standard wrapper method to straightforwardly simulate &#120068; &#8592; &#120069; paths using the chosen **TRAJ** method. The first passage path properties and the transition path bin statistics, printed to the output files *fpp\_properties.dat* and *tp\_stats.dat*, respectively, correspond to the _nonequilibrium_ path ensembles (i.e. standard first hitting problem).

**FIXEDT**  
  instructs the program to simulate a number of trajectories (equal to **NABPATHS**) of fixed time (equal to **TRAJT**), initialized at state &#120069;. Trajectories are *not* terminated when the absorbing state &#120068; is hit, and the output file *fpp\_properties.dat* is not written. Using the **STEADYSTATE** keyword, steady-state transition path bin statistics, as well as an estimate for the steady-state mean first passage time, are calculated for the _equilibrium_ &#120068; &#8592; &#120069; path ensemble, which exists if the Markov chain is irreducible. Simulation of the steady-state &#120068; &#8592; &#120069; path ensemble is best achieved using a small number of long-timescale trajectories.

**DIMREDN**  
  instructs the program to simulate many short trajectories (numbers specified via the **DIMREDUCTION** keyword) of fixed total time (specified via the **TRAJT** keyword) initialised from each community in turn. These trajectories are printed to files _walker.x.y.dat_, where _x_ is the ID of the community, and _y_ is the iteration number for that community. Trajectory information is written to files whenever a trajectory transitions to a new community. **DUMPINTVLS** must be set so that appropriate trajectory data is output. The simulation is parallelised, using a number of threads equal to **NTHREADS**. This calculation is compatible with two algorithms to propagate individual trajectories, namely, **TRAJ KPS**, and **TRAJ MCAMC** (without **MEANRATE**). The communities of nodes must be specified (**COMMSFILE** keyword). **NABPATHS**, **MAXIT**, and **BINSFILE** keywords are ignored. This setup is incompatible with specification of an initial condition via the **INITCONDFILE** keyword, and with the **NODESAFILE** and **NODESBFILE** keywords. Instead, a local equilibrium within the starting community is assumed as the initial probability distribution for each macrostate. A script to estimate a coarse-grained discrete- or continuous-time Markov chain from the relevant trajectory information (namely, the times at which communities are occupied) is available [here](https://github.com/danieljsharpe/DISCOTRESS_tools).

**WE**  
  the weighted ensemble method accelerates the sampling of &#120068; &#8592; &#120069; paths by using a splitting and culling procedure to maintain a specified number of walkers in each of the specified communties (given via the **COMMSTARGFILE** keyword if known *a priori*, otherwise the **ADAPTIVECOMMS** keyword must be set). Can only be used with **TRAJ BKL**.

**FFS**  
  the forward flux sampling method accelerates the sampling of &#120068; &#8592; &#120069; paths by ratcheting across nested interfaces.

**NEUS**  
  the non-equilibrium umbrella sampling method accelerates the sampling of &#120068; &#8592; &#120069; steady state paths by simulating walkers confined to cells. Can only be used with **TRAJ BKL**.

**MILES**  
  the milestoning method accelerates the sampling of &#120068; &#8592; &#120069; steady state paths by simulating walkers initialised at milestones (interfaces between macrostates) hitting adjacent milestones.

**REA**  
  the recursive enumeration algorithm (REA) determines the highest-probability &#120068; &#8592; &#120069; paths using a *k* shortest paths algorithm wherein the edge costs are given by the contributions of individual transitions to the total path action. There must be only a single initial (source) node and a single absorbing (sink) node (*cf*. the **NODESAFILE** and **NODESBFILE** keywords). The choice of **TRAJ** method option is arbitrary since an explicit simulation is not performed. **NABPATHS** is interpreted as the number of highest-probability paths to be computed (i.e. = *k*). If the **REANOTIRRED** keyword is specified, then the Markov chain is taken to be reducible, and the REA will not throw an error in the case that no candidate paths to a node exist (the default behaviour, suitable for irreducible Markov chains, is to throw an error in this circumstance). If no candidate paths to the target node can be found and the **REANOTIRRED** keyword is specified, then the program will exit the REA loop and print the set of paths that have been determined (which is then the complete set of A<-B paths). If the **WRITEREA** keyword is specified, then trajectory data for the *k* highest-probability paths are written to the files *shortest_path.k.dat* in the usual *walker.x.y.dat* format (see above), except that the paths are printed backwards. The output file *fpp_properties.dat* lists the properties of the dominant *k* first passage paths from the source to the sink node, stated in order of decreasing probability (increasing path action). For a DTMC (keyword **DISCRETETIME**), **NOLOOP** must be set, and for a CTMC (default), **BRANCHPROBS** must be set, so that shortest paths do not contain self-loop transitions for nodes. Hence, the entropy flow along shortest paths is not computed for DTMCs.

----

## Optional keywords relating to simulation parameters and output

The following is a list of keywords that specify generic simulation parameters and keywords that control the recording of statistics or other aspects of the output.

----

**BINSFILE** `str` `int`  
  optional. Default to be the same as **COMMSFILE**, if specified.
  Name of the file containing the bin IDs (indexed from 0) for nodes, and number of bins. The bins are used to collect statistics associated with nodes (or groups thereof) for the &#120068; &#8592; &#120069; transition path ensemble, namely committor and visitation probabilities.

**COMMSFILE** `str` `int`  
  mandatory if **WRAPPER** is **DIMREDN**, **FFS**, **NEUS**, or **MILES**. Also mandatory if **WRAPPER** is **WE**, or if **TRAJ** is **KPS** or **MCAMC**, and **ADAPTIVECOMMS** is not specified.
  Name of the file containing the definitions of communities (single-column, indexed from zero, number of entries equal to the number of nodes **NNODES** in the network) and no. of communities. Is overridden by **ADAPTIVECOMMS**. For both **WRAPPER** and **TRAJ** enhanced sampling methods, except **TRAJ BKL**, the communities are used to divide the state space (eg the communities define the trapping basins in **KPS**, or the communities for resampling in **WE**), and for certain algorithms may dictate the resolution at which the transition path statistics (see **BINSFILE** keyword) can be calculated. The specification of communities must be consistent with the definition of the &#120068; and &#120069; sets. An exception is if the number of communities is 2, in which case the initial set &#120069; can be a subset of the relevant community. Note that if this is chosen to be the case, then re-hitting &#120069; is not detected, and committor and transient visitation probabilities for the bins will be incorrect.

**DUMPINTVLS**  
  if set, then trajectory data is dumped at precisely the time intervals specified by **TINTVL** (which must therefore be >0.). Otherwise, trajectory data written when the next time interval is exceeded is precisely for the current time of the walker. Path probability and entropy flow are not written in the walker files if this keyword is set, but are still dumped to the *fpp_properties.dat* file. This keyword is required with **WRAPPER DIMREDN**, since the trajectory information required to construct a coarse-grained Markov chain is otherwise not printed.

**INITCONDFILE** `str`  
  optional. Name of the file containing initial occupation probabilities for nodes in &#120069; that are alternative to the stationary probabilities. The number of entries is assumed to be the same as the specified number of nodes in &#120069;, and the specified order is assumed to be the same also \(cf **NODESBFILE**\). The values must sum to unity.

**MAXIT** `int`  
  default is inf. The maximum number of iterations of the relevant algorithm to run before the simulation is terminated (if the target number of &#120068; &#8592; &#120069; paths to simulate is not reached). The interpretation of this option depends on the chosen enhanced sampling method. e.g. with **WRAPPER WE**, **MAXIT** is the number of iterations of the resampling procedure. With **WRAPPER BTOA** and **TRAJ KPS** or **TRAJ MCAMC**, **MAXIT** is the number of basin escape trajectories simulated.

**NABPATHS** `int`  
  mandatory if not **WRAPPER DIMREDN** and if none of the state reduction keywords are specified. The simulation is terminated when this number of &#120068; &#8592; &#120069; paths have been successfully sampled. If **WRAPPER FIXEDT**, then this number is the number of paths of fixed total time to be simulated (not necessarily conditioned on the endpoint &#120068; and &#120069; states).

**TINTVL** `double`  
  time interval for dumping trajectory information. Negative value (default) indicates that trajectory data is not written (i.e. files _walker.0.y.dat_ are not output). Zero value specifies that all trajectory information is written. An explicit non-negative value must be set if **WRAPPER DIMREDN**. The exact value of **TINTVL** is ignored if **TRAJ KPS** (in which case trajectory data is written after every basin escape).

----

## Optional keywords relating to enhanced sampling methods

The following is a list of keywords that specify simulation parameters pertaining to particular enhanced sampling methods (selected by the **WRAPPER** and **TRAJ** keywords).

----

**ADAPTIVECOMMS** `double`  
  mandatory if **WRAPPER WE**, **TRAJ KPS**, or **TRAJ MCAMC**, and **COMMSFILE** is not specified. Default _False_.
  Set the partitioning of the state space leveraged in **WE**, **KPS** or **MCAMC** to be defined on-the-fly by a breadth-first search procedure. The argument is the minimum transition rate for a node to be included in the community being built up. If set with **TRAJ** as **KPS** or **MCAMC**, **KPSKMCSTEPS** is ignored.

**COMMSTARGFILE** `str`  
  mandatory if **WRAPPER WE** and not **ADAPTIVECOMMS**.
  Name of the file containing the target number of trajectories in each bin (single-column, number of entries equal to the number of communities in the network).

**DIMREDUCTION** `str`  
  mandatory if **WRAPPER DIMREDN**, which initialises a special wrapper class that does not perform the usual code function, which is to simulate &#120068; &#8592; &#120069; transition paths, and instead instructs the program to simulate many short trajectories starting from each community, each of length in time equal to **TRAJT**. The total number of trajectories that are to be simulated starting from each community is listed in the file given as the string arg (single-column format, length equal to number of communities, set via the **COMMSFILE** keyword).

**KPSKMCSTEPS** `int`  
  optional. If **TRAJ** is **KPS** or **MCAMC**, specifies the number of standard BKL steps to be performed after a kPS or MCAMC escape from a trapping basin. Default is 0 (pure kPS (or MCAMC), no kMC steps). However, this is not the recommended value. If using **TRAJ KPS** or **TRAJ MCAMC**, for most systems, great gains in simulation efficiency will be achieved by setting **KPSKMCSTEPS** to an appropriate nonzero value. This is because many metastable systems will feature transition regions between metastable states. Therefore, after each basin escape, the trajectory will likely flicker between the two basins. Rather than simulate expensive kPS or MCAMC basin escape iterations for these trivial recrossings, it is much more efficient to perform standard BKL steps. Note that this keyword does not require **BRANCHPROBS** to be set, and can also be used with **DISCRETETIME**. Ignored if **ADAPTIVECOMMS**.

**MEANRATE**  
  optional. If **TRAJ MCAMC**, the calculation uses the approximate mean rate method, as opposed to the default exact first passage time analysis (FPTA) method. Default false.

**NELIM** `int`  
  mandatory if **TRAJ KPS**. The maximum number of nodes that are to be eliminated from the current trapping basin. If **NELIM** exceeds the number of nodes in the largest community, then all states of any trapping basin are always eliminated. Note that **NELIM** determines the number of transition matrices stored for the active subnetwork, and therefore the choice of this keyword (along with the sizes of communities) can strongly affect memory usage.

**NWALKERS** `int`  
  mandatory if **WRAPPER** is **WE**, **FFS**, **NEUS**, or **MILES**. Specifies the number of walkers (independent trajectories) on the network, which are simulated in parallel (see **NTHREADS**). This keyword is ignored (and therefore does not need to be explicitly set) if **WRAPPER** is **BTOA** or **DIMREDN**, in which case the number of walkers is set to **NTHREADS**.

**REANOTIRRED**  
  if **WRAPPER REA**, specifies that candidate paths to nodes may not necessarily exist (this situation may occur when the Markov chain is not irreducible). Hence, errors are not thrown in this circumstance (unlike the default behaviour), and the main loop of the REA is exited in the event that no more paths to the target node exist. Default false.

**STEADYSTATE** `double`  
  optional. If **WRAPPER FIXEDT**, indicates that a small number of trajectories (equal to **NTHREADS**) of fixed total time are to be ran, from which statistics for the &#120068; &#8592; &#120069; *equilibrium* (steady state) TPE are to be computed. The argument associated with this keyword specifies the time threshold after which the trajectory is considered to have equilibriated and recording of steady state path statistics begins. The default value for this argument is 0., but this value should be altered to an appropriate finite value. To ensure that the simulation estimates of these steady state properties are unbiased and accurate, the total fixed time of trajectories (set by **TRAJT**) should be long, to ensure that sufficient statistics are obtained, and statistics should be recorded after a suitably long time period has passed (several times the average mixing time [Kemeny constant] of the Markov chain), to ensure that the trajectories have equilibriated prior to recording steady state path statistics.

**TAURE** `double`  
  mandatory if **WRAPPER WE**. The time between resampling trajectories.

**TRAJT** `long double`  
  mandatory if **WRAPPER** is **FIXEDT** or **DIMREDN**. The maximum time for trajectories when simulating paths of fixed total time.

**WRITEREA**  
  if **WRAPPER REA**, write output trajectory files *shortest_path.k.dat*, in the usual *walker.x.y.dat* format (see above) except backwards, for each of the *k* shortest paths. Default false.

----

## Optional keywords related to exact numerical analysis of the dynamics by state reduction methods

The following keywords are used in combination with the keywords **WRAPPER BTOA** and **TRAJ KPS**. Use of any of the following keywords overrides the default functionality of DISCOTRESS, which is to simulate dynamical paths, and instead instructs the program to perform a state reduction procedure to exactly compute one or more dynamical quantities associated with nodes, in a numerically stable manner. **NABPATHS** must be set to 1. The _communities.dat_ file must specify precisely two communities; namely, nodes in the target set &#120068; and nodes not in &#120068;. The **COMMITTOR**, **ABSORPTION**, **MFPT**, and **GTH** keywords can be used together in any combination. The computations performed with the **FUNDAMENTALRED** and **FUNDAMENTALIRRED** keywords are standalone operations.

The memory costs of state reduction computations can be reduced without consequence by setting the type of the `h` members (which represent the numbers of kMC steps for transitions, when using the **KPS** algorithm) of the `Node` and `Edge` structures (defined in the file *ktn.h*) to `int`, since these members are not used in the state reduction methods.

----

**ABSORPTION**  
  specifies that a state reduction procedure is performed to compute the absorption probabilities. The probabilities B\_ij that a trajectory initialised from the non-absorbing node _i_ is absorbed at node _j_ are written to the file *absorption.dat* in the format "_i_ / _j_ / *B\_ij*". For the initial occupation probability distribution (which, by default, is assumed to be a local equilibrium within the initial set &#120069;), the absorption (hitting) probabilities for each absorbing node are printed to the file *hitting\_probs.dat*.

**COMMITTOR**  
  specifies that a state reduction procedure is performed to compute the committor probabilities for all nodes. The committor probabilities are written to the files *committor\_AB.dat* and *committor\_BA.dat* (for &#120068; &#8592; &#120069; and &#120069; &#8592; &#120068; directions, respectively). Note that the committor probabilities determined by this method are for each node, and the calculation is exact and deterministic (unlike calculation of the committor probabilities for the bins from simulation data, cf. the **BINSFILE** keyword).

**FUNDAMENTALIRRED**  
  specifies that a state reduction algorithm is used to compute the fundamental matrix of an irreducible Markov chain. The (_i_,_j_) edge weights *Z\_ij* of the renormalised network resulting from this procedure are the elements of the fundamental matrix, the trace of which gives the Kemeny constant (average mixing time) for the Markov chain. These values are written to the file *fundamental.dat* in the format "_i_ / _j_ / *Z\_ij*", and the Kemeny constant is printed in the output.

**FUNDAMENTALRED**  
  specifies that a state reduction algorithm is used to compute the fundamental matrix of an absorbing (i.e. reducible) Markov chain. The (_i_,_j_) edge weights *n\_ij* of the renormalised network resulting from this procedure are the expected numbers of times that the _j_-th node is visited along first passage paths initialised from the _i_-th node. These values are written to the file *transient\_visits.dat* in the format "_i_ / _j_ / *n\_ij*". The node visitation probabilities can be computed from this information if the committor probabilities are also known. For the initial occupation probability distribution (which, by default, is assumed to be a local equilibrium within the initial set &#120069;), the expected numbers of times that non-absorbing nodes are visited along first passage paths are printed to the file *node\_visits.dat*.

**GTH**  
  specifies that the stationary probability distribution (which exists if the Markov chain is irreducible) is computed using the Grassmann-Taksar-Heyman (GTH) algorithm. Can only be used when the target set &#120068; contains a single node. The input file *stat\_prob.dat* must be provided, but its contents are not used. The stationary probabilities determined by the GTH algorithm are written to the file *stat\_prob\_gth.dat*.

**MFPT**  
  specifies that a state reduction procedure is performed to compute mean first passage times (MFPT). The MFPTs *m\_i*&#120068; for transitions from non-absorbing nodes _i_ to the set of absorbing nodes &#120068; are written to the file *mfpt.dat* in the format "_i_ / *m\_i*&#120068;". Given an initial occupation probability distribution (which, by default, is assumed to be a local equilibrium within the initial set &#120069;), the &#120068; &#8592; &#120069; MFPT is printed in the output. If the initial mean waiting times of nodes are set to the initial mean number of steps to exit (i.e. equal to unity for all nodes), then the MFPTs are in fact the mean first passage path lengths.

**PATHLENGTHS**  
  when used in conjunction with **MFPT**, specifies that mean first passage path lengths (instead of times) are computed. This is achieved by overriding the mean waiting times to instead represent the mean numbers of steps to exit, which are initially equal to unity for all nodes. Is used in conjunction with **BRANCHPROBS**, in which case each transition represents a move to a different node.

----

## Other optional keywords

----

**ACCUMPROBS**  
  if **TRAJ BKL**, the edges for transitions from each node are ordered according to decreasing transition probability. This optimizes the performance of the BKL algorithm, so is generally recommended, but the path entropy flow is then not output. Default false.

**BRANCHPROBS**  
  when simulating a CTMC, this keyword indicates that the transition probabilities used internally in the program are the branching probabilities. In this case, there are no self-loops and the mean waiting times are uniform. Otherwise, the linearised transition probability matrix is used, and **TAU** must be set. The **TRAJ BKL** and **TRAJ KPS** methods are more efficient when the branching probabilities are used, so this keyword is generally recommended. This keyword is ignored if **TRAJ MCAMC**. This keyword is not compatible with **DISCRETETIME**.

**DEBUG**  
  enable extra printing and tests to aid debugging. Default false.

**DISCRETETIME**  
  the edge weights read from *edge\_weights.dat* are taken to be the transition probabilities of a DTMC. This overrides the default behaviour, which is to assume that *edge\_weights.dat* is a list of (log) transition rates parameterising a continuous-time Markov chain (CTMC). Only the probabilities for transitions between *different* nodes need to be specified in the *edge\_weights.dat* file. The self-loop transition probabilities for nodes are inferred to be the difference of the sum of probabilities for outgoing transitions from unity. For a DTMC, **TAU** must be provided, and is interpreted as the (fixed) lag time (i.e. all transitions are associated with a constant time step **TAU**, rather than an exponential distribution with mean **TAU**). Note that this keyword is compatible with all **TRAJ** options and state reduction procedures.

**DUMPWAITTIMES**  
  dump the mean waiting times for nodes to the file _meanwaitingtimes.dat_.

**NOLOOP**  
  if **DISCRETETIME**, the average numbers of self-loop transitions for nodes are accounted for implicitly by renormalization of outgoing transition probabilities and the lag time. Thus the lag time for transitions from nodes becomes node-dependent, and represents an *expectation* with respect to the numbers of self-loop transitions before escape from a node. The length of a path then represents the number of transitions between *different* nodes (as is the case for a CTMC parameterized by a branching probability matrix), often referred to as the *dynamical activity*. When using **TRAJ BKL**, this keyword will increase the efficiency of the simulation, since the self-loop transitions for nodes are not explicitly taken. However, when using this option, only the mean of the simulated first passage time distribution is meaningful. Not compatible with **TRAJ MCAMC**. Default false.

**NTHREADS** `int`  
  number of threads to use in parallel calculations. Defaults to max. no. of threads available. Keyword is overridden and set equal to one when performing a state reduction computation.

**SEED** `int`  
  seed for the random number generators (default 19).

**SEMIMARKOV**  
  indicates that the waiting time distributions for internode transitions in a continuous-time model are not exponential distributions but are instead Weibull distributions. The two Weibull distribution parameters are read from input files. The first parameter overrides the `t_esc` member of the `node` class, which otherwise represents the mean waiting time for a node in a CTMC (or the lag time for a node in a DTMC). Recall that the exponential distribution has the memoryless property, and therefore defines a CTMC. A continuous-time process for which the transition probabilities depend only on the current node, and for which the waiting time distributions are non-exponential, is a semi-Markov process. DISCOTRESS can be used to simulate an arbitrary finite semi-Markov chain by replacing the function `weibull_distribn()` representing the Weibull distribution with any probability distribution of choice. This keyword is not compatible with **TRAJ MCAMC** or **DISCRETETIME**, and is not compatible with any state reduction procedures. [This keyword is not yet implemented].

**TAU** `long double`  
  mandatory if not **BRANCHPROBS**. If **DISCRETETIME**, **TAU** is the lag time at which the DTMC is parameterised. Otherwise, if **BRANCHPROBS** is not provided, then the CTMC is parameterised by a linearised transition probability matrix with **TAU** the uniform mean waiting time.

----

