# DISCOTRESS user information

## Compilation
Compile with:
g++ -std=c++17 discotress.cpp kmc\_methods.cpp we\_kmc.cpp ffs\_kmc.cpp neus\_kmc.cpp milestoning.cpp kps.cpp mcamc.cpp keywords.cpp ktn.cpp -o discotress -fopenmp

## Input files 

mandatory input files containing the network topology information:
-  *stat_prob.dat*   \- single col, length **NNODES**: natural log of stationary probabilities for nodes
-  *ts_conns.dat*    \- double col, length **NEDGES**: edge connecting node IDs (each is assumed to be bidirectional)
-  *ts_weights.dat*  \- single col, length 2\***NEDGES**: natural log of transition rate for forward and reverse transitions (two lines per line in ts\_conns.dat)
-  *nodes.A*         \- single col, see **NODESAFILE** keyword: node IDs in the absorbing set A (indexed from 1)
-  *nodes.B*         \- single col, see **NODESBFILE** keyword: node IDs in the initial set B (indexed from 1)

optional input files relating to the simulation and calculating statistics:
-  *communities.dat* \- single col, length **NNODES**, see **COMMSFILE** keyword: community IDs (indexed from 0) for nodes
-  *bins.dat*        \- single col, length **NNODES**: bin IDs for nodes (indexed from 0), used for calculating transition path statistics (defaults to communities.dat)
-  *initcond.dat*    \- single col, length cf **NODESBFILE**, see **INITCONDFILE** keyword: initial probability distribution for nodes of the initial set B. Defaults to be proportional to a local equilibrium within B
-  *ntrajs.dat*      \- single col, length cf **COMMSFILE**, see **DIMREDUCTION** keyword: numbers of short trajectories to run from each community

keyword file:
-  *input.kmc*       \- see below for a list of keywords

## Output files
-  *walker.x.y.dat*  \- trajectory information dumped at the specified time intervals (or when a trajectory escapes from a community, depending on options). x is the walker ID, y is the path number.  
Format: node ID (indexed from 1) / community ID (indexed from 0) / path length (no. of kMC moves) / path time / natural log of path probability / path entropy flow
-  *tp_distribns.dat* \- properties of A<-B transition paths, from which the distributions can be analysed.  
Format: transition path no. / path length / path time / natural log of path probability / path entropy flow
-  *tp_stats.dat*    \- bin statistics, written if communities were specified.  
Format: bin ID (indexed from 0) / no. of direct A<-B paths for which bin is visited / no. of paths for which bin is visited and trajectory returned to initial set B / transition path probability / committor probability

## Main keywords

**NNODES** _int_  
  mandatory, number of nodes in the transition network

**NEDGES** _int_   
  mandatory, number of (bidirectional) edges in the transition network

**WRAPPER** _str_  
  mandatory, employ a \`wrapping' enhanced sampling strategy for accelerating the observation of rare events. Options:  
-    **NONE**    \- no enhanced sampling strategy employed
-    **DIMREDN** \- special class to simulate many short nonequilibrium trajectories starting from each community in turn
               used for estimation of a coarse-grained transition network. See **DIMREDUCTION** keyword for more detail 
-    **WE**      \- weighted ensemble sampling
-    **FFS**     \- forward flux sampling
-    **NEUS**    \- non-equilibrium umbrella sampling
-    **MILES**   \- milestoning

**TRAJ** _str_  
  mandatory, method for propagating individual kMC trajectories. Options:  
-    **BKL**     \- Bortz-Kalos-Lebowitz \(aka n-fold way\) rejection-free algorithm \(ie standard kMC\)
-    **KPS**     \- kinetic path sampling algorithm
-    **MCAMC**   \- Monte Carlo with Absorbing Markov Chains algorithm
  Note the following excluded combinations of WRAPPER / TRAJ options: WE / KPS, WE / MCAMC, NEUS / KPS, NEUS / MCAMC, DIMREDN / BKL

**NODESAFILE** _str_ _int_  
  mandatory, name of the file containing the node ID's (indexed from 1) belonging to the A (absorbing) set and no. of nodes in A

**NODESBFILE** _str_ _int_  
  mandatory, name of the file containing the node ID's (indexed from 1) belonging to the B (initial) set and no. of nodes in B

**NABPATHS** _int_  
  mandatory, the simulation is terminated when this number of A<-B paths have been successfully sampled

## Optional keywords relating to enhanced sampling methods

**INITCONDFILE** _str_  
  optional. Name of the file containing initial occupation probabilities for nodes in B that are alternative to the stationary probabilities. The number of entries is assumed to be the same as the specified number of nodes in B, and the specified order is assumed to be the same also \(cf **NODESBFILE**\). Values must sum to unity

**COMMSFILE** _str_ _int_  
  mandatory if **WRAPPER** is DIMREDN, FFS, NEUS, or MILES. Also mandatory if **WRAPPER** is WE, or if **TRAJ** is KPS or MCAMC, and **ADAPTIVECOMMS** is not specified.
  Name of the file containing the definitions of communities (single-column, indexed from zero, number of entries equal to the number of nodes **NNODES** in the network) and no. of communities. Is overridden by **ADAPTIVECOMMS**. For both **WRAPPER** and **TRAJ** enhanced sampling methods, except **TRAJ BKL**, the communities are used to divide the state space (eg the communities define the trapping basins in KPS, or the communities for resampling in WE), and for certain algorithms may dictate the resolution at which the transition path statistics (see **BINFILE** keyword) can be calculated. The specification of communities must be consistent with the definition of the A and B sets. An exception is if the number of communities is 2, in which case the initial set B can be a subset of the relevant community. Note that if this is chosen to be the case, then re-hitting B is not detected, and committor probabilities and transition path probability density will be incorrect.

**ADAPTIVECOMMS** _double_  
  mandatory if **WRAPPER WE**, **TRAJ KPS**, or **TRAJ MCAMC**, and **COMMSFILE** is not specified. Default False.
  Set the partitioning of the state space leveraged in WE-kMC, kPS or MCAMC to be defined on-the-fly by a breadth-first search procedure. The argument is the minimum transition rate for a node to be included in the community being built up. If set with KPS or MCAMC, **KPSKMCSTEPS** is ignored.

**COMMSTARGFILE** _str_  
  mandatory if **WRAPPER WE** and not **ADAPTIVECOMMS**.
  Name of the file containing the target number of trajectories in each bin (single-column, number of entries equal to the number of communities in the network).

**BINFILE** _str_ _int_  
  optional. Default to be the same as **COMMSFILE**, if specified.
  Name of the file containing the definition of bins (indexed from 0) for calculating transition path statistics (ie committor functions and transition path probability densities), and no. of bins.

**TAU** _long double_  
  mandatory if **WRAPPER WE** or if **TRAJ KPS** and not **BRANCHPROBS**.
  If WE, tau is the time between resampling trajectories. If KPS, tau is the lag time at which the linear transition probability matrices are estimated (or provided, if **TRANSNPROBS**).

**KPSKMCSTEPS** _int_  
  optional. If **TRAJ** is **KPS** or **MCAMC**, specifies the number of standard BKL steps to be performed after a kPS or MCAMC escape from a trapping basin. Default is 0 (pure kPS (or MCAMC), no kMC steps). Requires **BRANCHPROBS** to be set. Ignored if **ADAPTIVECOMMS**.

**MEANRATE**  
  optional. If **TRAJ MCAMC**, the calculation uses the approximate mean rate method, as opposed to the default exact first passage time analysis (FPTA) method

**NELIM** _int_  
  mandatory if **TRAJ KPS**. The maximum number of nodes that are to be eliminated from the current trapping basin. If **NELIM** exceeds the number of nodes in the largest community, then all states of any trapping basin are always eliminated. Note that **NELIM** determines the number of transition matrices stored for the active sub-network, and therefore the choice of this keyword (along with the sizes of communities) can strongly affect memory usage.

**DIMREDUCTION** _str_ _double_  
  mandatory if **WRAPPER DIMREDN**, which initialises a special wrapper class that does not perform the usual code function, which is to simulate A<-B transition paths, and instead instructs the program to simulate many short trajectories starting from each community, each of length in time equal to the float argument. These trajectories are printed to files _walker.x.y.dat_, where _x_ is the ID of the community, and _y_ is the iteration number for that community. Trajectory information is written to files whenever a trajectory transitions to a new community. The total number of trajectories that are to be simulated starting from each community is listed in the file given as the string arg. This calculation is parallelised, using a number of threads equal to **NTHREADS** (defaults to maximum number of threads available). This calculation is compatible with two algorithms to propagate individual trajectories: **TRAJ KPS**, and **TRAJ MCAMC** (without **MEANRATE**). The communities of nodes must be specified (**COMMSFILE** keyword). **NODESAFILE**, **NODESBFILE**, **NABPATHS**, **MAXIT**, and **BINFILE** keywords are ignored. This setup is incompatible with specification of an initial condition via the **INITCONDFILE** keyword. Instead, a local equilibrium within the starting community is assumed as the initial probability distribution for each trajectory.

## Other optional keywords

**MAXIT** _int_  
  default is inf. The maximum number of iterations of the relevant algorithm to run before the simulation is terminated (if the target number of A<-B paths to simulate is not reached). The interpretation of this option depends on the chosen enhanced sampling method. eg with WRAPPER WE, **MAXIT** is the number of iterations of the resampling procedure. With WRAPPER NONE and TRAJ KPS or TRAJ MCAMC, **MAXIT** is the number of basin escape trajectories simulated.

**NTHREADS** _int_  
  number of threads to use in parallel calculations. Defaults to max. no. of threads available.

**BRANCHPROBS**  
  the transition probabilities are calculated as the branching probabilities. Mandatory if **TRAJ** is **BKL**, or is **KPS** or **MCAMC** and **KPSKMCSTEPS** is specified.

**TRANSNPROBS**  
  the edge weights read in from the file *ts_weights.dat* are transition probabilities, not rates (compatible with **TRAJ KPS** only). The self-loop transition probabilities for nodes are inferred to be the remainders from unity for transitions from each node.

**PFOLD**  
  when used in conjunction with TRAJ KPS, specifies that a committor function calculation is performed instead of a kPS simulation. The _communities.dat_ file must specify precisely three communities; the A set, the B set, and the set of all other nodes ("I"). The committor functions are written to the files "committor\_AB.dat" and "committor\_BA.dat" (for A<-B and B<-A directions, respectively). Note that the committor functions determined by this method are for each node, and the calculation is _exact_ and _deterministic_ (unlike calculation of the committor probabilities for the bins from simulation data, cf the **BINSFILE** keyword). **NABPATHS** must be set to some arbitrary number >0.

**TINTVL** _double_  
  ignored if **TRAJ KPS** (in which case trajectory data is written after every basin escape). Time interval for dumping trajectory information. Negative value (default) indicates that trajectory data is not written. In this case, start and end nodes are still written to files walker.0.y.dat. Zero value specifies that all trajectory information is written.

**SEED** _int_  
  seed for the random number generators (default 19).

**DEBUG**  
  enable extra printing and tests to aid debugging

**DUMPWAITTIMES**  
  dump the mean waiting times for nodes to the file _meanwaitingtimes.dat_

## Keywords relating to nonlinear master equations

**DISCOTRESS** can be used to interpret a set of rules relating to a nonlinear master equation and map the system to a linear transition network.
Additional files are needed for this purpose:

  _changevecs.dat_  - ...
