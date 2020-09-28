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
-  *input.kmc*       \- see below for a list of keywords. Comment lines begin with an exclamation mark (!)

## Output files
-  *walker.x.y.dat*  \- trajectory information dumped at the specified time intervals (or when a trajectory escapes from a community, depending on options). x is the walker ID, y is the path number.  
Format: node ID (indexed from 1) / community ID (indexed from 0) / path length (no. of kMC moves) / path time / natural log of path probability / path entropy flow
-  *tp_distribns.dat* \- properties of A<-B transition paths, from which the distributions can be analysed.  
Format: transition path no. / path length / path time / natural log of path probability / path entropy flow
-  *tp_stats.dat*    \- bin statistics, written if communities were specified.  
Format: bin ID (indexed from 0) / no. of direct A<-B paths for which bin is visited / no. of paths for which bin is visited and trajectory returned to initial set B / transition path probability / committor probability

## Main keywords

**NNODES** _int_  
  mandatory, number of nodes in the Markov chain.

**NEDGES** _int_   
  mandatory, number of (bidirectional) edges in the Markov chain.

**WRAPPER** _str_  
  mandatory, employ a \`wrapping' enhanced sampling strategy for accelerating the observation of rare events. Options:  
-    **NONE**    \- no enhanced sampling strategy employed
-    **DIMREDN** \- special class to simulate many short nonequilibrium trajectories starting from each community in turn
               used for estimation of a coarse-grained Markov chain. See **DIMREDUCTION** keyword for more detail 
-    **WE**      \- weighted ensemble sampling
-    **FFS**     \- forward flux sampling
-    **NEUS**    \- non-equilibrium umbrella sampling
-    **MILES**   \- milestoning

**TRAJ** _str_  
  mandatory, method for propagating individual kMC trajectories. Options:  
-    **BKL**     \- Bortz-Kalos-Lebowitz \(aka n-fold way\) rejection-free algorithm \(ie standard kMC\)
-    **KPS**     \- kinetic path sampling algorithm
-    **MCAMC**   \- Monte Carlo with Absorbing Markov Chains algorithm  
  Note the following excluded combinations of WRAPPER / TRAJ options: WE / KPS, WE / MCAMC, NEUS / KPS, NEUS / MCAMC, DIMREDN / BKL.

**NODESAFILE** _str_ _int_  
  mandatory, name of the file containing the node ID's (indexed from 1) belonging to the A (absorbing) set and no. of nodes in A.

**NODESBFILE** _str_ _int_  
  mandatory, name of the file containing the node ID's (indexed from 1) belonging to the B (initial) set and no. of nodes in B.

**NABPATHS** _int_  
  mandatory, the simulation is terminated when this number of A<-B paths have been successfully sampled.

## Optional keywords relating to enhanced sampling methods

**INITCONDFILE** _str_  
  optional. Name of the file containing initial occupation probabilities for nodes in B that are alternative to the stationary probabilities. The number of entries is assumed to be the same as the specified number of nodes in B, and the specified order is assumed to be the same also \(cf **NODESBFILE**\). The values must sum to unity.

**COMMSFILE** _str_ _int_  
  mandatory if **WRAPPER** is **DIMREDN**, **FFS**, **NEUS**, or **MILES**. Also mandatory if **WRAPPER** is **WE**, or if **TRAJ** is **KPS** or **MCAMC**, and **ADAPTIVECOMMS** is not specified.
  Name of the file containing the definitions of communities (single-column, indexed from zero, number of entries equal to the number of nodes **NNODES** in the network) and no. of communities. Is overridden by **ADAPTIVECOMMS**. For both **WRAPPER** and **TRAJ** enhanced sampling methods, except **TRAJ BKL**, the communities are used to divide the state space (eg the communities define the trapping basins in **KPS**, or the communities for resampling in **WE**), and for certain algorithms may dictate the resolution at which the transition path statistics (see **BINFILE** keyword) can be calculated. The specification of communities must be consistent with the definition of the A and B sets. An exception is if the number of communities is 2, in which case the initial set B can be a subset of the relevant community. Note that if this is chosen to be the case, then re-hitting B is not detected, and committor and transient visitation probabilities for the bins will be incorrect.

**ADAPTIVECOMMS** _double_  
  mandatory if **WRAPPER WE**, **TRAJ KPS**, or **TRAJ MCAMC**, and **COMMSFILE** is not specified. Default _False_.
  Set the partitioning of the state space leveraged in **WE**, **KPS** or **MCAMC** to be defined on-the-fly by a breadth-first search procedure. The argument is the minimum transition rate for a node to be included in the community being built up. If set with **TRAJ** as **KPS** or **MCAMC**, **KPSKMCSTEPS** is ignored.

**COMMSTARGFILE** _str_  
  mandatory if **WRAPPER WE** and not **ADAPTIVECOMMS**.
  Name of the file containing the target number of trajectories in each bin (single-column, number of entries equal to the number of communities in the network).

**BINFILE** _str_ _int_  
  optional. Default to be the same as **COMMSFILE**, if specified.
  Name of the file containing the definition of bins (indexed from 0) for calculating transition path statistics (i.e. committor and transient visitation probabilities), and no. of bins.

**TAU** _long double_  
  mandatory if not **BRANCHPROBS**.
  **TAU** is the lag time at which the linear transition probability matrices are estimated (or provided, if **TRANSNPROBS**). If **TRANSPROBS** and **DISCRETETIME** are specified, then **TAU** is the lag time of the discrete-time Markov chain.

**KPSKMCSTEPS** _int_  
  optional. If **TRAJ** is **KPS** or **MCAMC**, specifies the number of standard BKL steps to be performed after a kPS or MCAMC escape from a trapping basin. Default is 0 (pure kPS (or MCAMC), no kMC steps). However, this is not the recommended value. If using **TRAJ KPS** or **TRAJ MCAMC**, for most systems, great gains in simulation efficiency will be achieved by setting **KPSKMCSTEPS** to an appropriate nonzero value. This is because many metastable systems will feature transition regions between metastable states. Therefore, after each basin escape, the trajectory will likely flicker between the two basins. Rather than simulate expensive kPS or MCAMC basin escape iterations for these trivial recrossings, it is much more efficient to perform standard BKL steps. Note that this keyword does not require **BRANCHPROBS** to be set, and can also be used with **DISCRETETIME**. Ignored if **ADAPTIVECOMMS**.

**MEANRATE**  
  optional. If **TRAJ MCAMC**, the calculation uses the approximate mean rate method, as opposed to the default exact first passage time analysis (FPTA) method.

**NELIM** _int_  
  mandatory if **TRAJ KPS**. The maximum number of nodes that are to be eliminated from the current trapping basin. If **NELIM** exceeds the number of nodes in the largest community, then all states of any trapping basin are always eliminated. Note that **NELIM** determines the number of transition matrices stored for the active subnetwork, and therefore the choice of this keyword (along with the sizes of communities) can strongly affect memory usage.

**TAURE** _double_
  mandatory if **WRAPPER WE**, the time between resampling trajectories

**DIMREDUCTION** _str_ _long double_  
  mandatory if **WRAPPER DIMREDN**, which initialises a special wrapper class that does not perform the usual code function, which is to simulate A<-B transition paths, and instead instructs the program to simulate many short trajectories starting from each community, each of length in time equal to the float argument. These trajectories are printed to files _walker.x.y.dat_, where _x_ is the ID of the community, and _y_ is the iteration number for that community. Trajectory information is written to files whenever a trajectory transitions to a new community. The total number of trajectories that are to be simulated starting from each community is listed in the file given as the string arg. This calculation is parallelised, using a number of threads equal to **NTHREADS** (defaults to maximum number of threads available). This calculation is compatible with two algorithms to propagate individual trajectories: **TRAJ KPS**, and **TRAJ MCAMC** (without **MEANRATE**). The communities of nodes must be specified (**COMMSFILE** keyword). **NABPATHS**, **MAXIT**, and **BINFILE** keywords are ignored. This setup is incompatible with specification of an initial condition via the **INITCONDFILE** keyword, and with the **NODESAFILE** and **NODESBFILE** keywords. Instead, a local equilibrium within the starting community is assumed as the initial probability distribution for each macrostate. Note that **WRAPPER DIMREDN** requires both this keyword and **DUMPINTVLS** to be set.

## Optional keywords related to exact numerical analysis of the dynamics by state reduction methods

The following keywords are used in combination with the keywords **WRAPPER NONE** and **TRAJ KPS**. Use of any of the following keywords overrides the default functionality of DISCOTRESS, which is to simulate dynamical paths, and instead instructs the program to perform a state reduction procedure to exactly compute one or more dynamical quantities associated with nodes, in a numerically stable manner. **NABPATHS** must be set to some arbitrary number >0. The _communities.dat_ file must specify precisely two communities; namely, nodes in the target set A and nodes not in A. The **COMMITTOR**, **ABSORPTION**, **MFPT**, and **GTH** keywords can be used together in any combination. The computations performed with the **FUNDAMENTALRED** and **FUNDAMENTALIRRED** keywords are standalone operations.

The memory costs of state reduction computations can be reduced by setting the type of the "h" members (which represent the numbers of kMC steps for transitions, when using the **KPS** algorithm) of the Node and Edge structures (defined in the file *ktn.h*) to _int_, since these members are not used in the state reduction methods.

**COMMITTOR**  
  specifies that a state reduction procedure is performed to compute the committor probabilities for all nodes. The committor probabilities are written to the files *committor\_AB.dat* and *committor\_BA.dat* (for A<-B and B<-A directions, respectively). Note that the committor probabilities determined by this method are for each node, and the calculation is exact and deterministic (unlike calculation of the committor probabilities for the bins from simulation data, cf. the **BINFILE** keyword).

**ABSORPTION**  
  specifies that a state reduction procedure is performed to compute the absorption probabilities. The probabilities b\_ij that a trajectory initialised from the non-absorbing node i is absorbed at node j are written to the file *absorption.dat* in the format "i / j / b\_ij". For the initial occupation probability distribution (which, by default, is assumed to be a local equilibrium within the initial set B), the absorption (hitting) probabilities for each absorbing node are printed to the file *hitting\_probs.dat*.

**FUNDAMENTALRED**  
  specifies that a state reduction algorithm is used to compute the fundamental matrix of an absorbing (i.e. reducible) Markov chain. The (i,j) edge weights n\_ij of the renormalised network resulting from this procedure are the expected numbers of times that the j-th node is visited along first passage paths initialised from the i-th node. These values are written to the file *transient\_visits.dat* in the format "i / j / n\_ij". The node visitation probabilities can be computed from this information if the committor probabilities are also known. For the initial occupation probability distribution (which, by default, is assumed to be a local equilibrium within the initial set B), the expected numbers of times that non-absorbing nodes are visited along first passage paths are printed to the file *node\_visits.dat*.

**FUNDAMENTALIRRED**  
  specifies that a state reduction algorithm is used to compute the fundamental matrix of an irreducible Markov chain. The (i,j) edge weights z\_ij of the renormalised network resulting from this procedure are the elements of the fundamental matrix, the trace of which gives the Kemeny constant (average mixing time) for the Markov chain. These values are written to the file *fundamental.dat* in the format "i / j / z\_ij", and the Kemeny constant is printed in the output.

**MFPT**  
  specifies that a state reduction procedure is performed to compute mean first passage times (MFPT). The MFPTs m\_iA for transitions from non-absorbing nodes i to the set of absorbing nodes A are written to the file *mfpt.dat* in the format "i / m\_iA". Given an initial occupation probability distribution (which, by default, is assumed to be a local equilibrium within the initial set B), the A<-B MFPT is printed in the output. If the initial mean waiting times of nodes are set to the initial mean number of steps to exit (i.e. equal to unity for all nodes), then the MFPTs are in fact the mean first passage path lengths.

**GTH**  
  specifies that the stationary probability distribution (which exists if the Markov chain is irreducible) is computed using the Grassmann-Taksar-Heyman (GTH) algorithm. Can only be used when the target set A contains a single node. The input file *stat\_prob.dat* must be provided, but its contents are not used. The stationary probabilities determined by the GTH algorithm are written to the file *stat\_prob\_gth.dat*.

**PATHLENGTHS**  
  when used in conjunction with **MFPT**, specifies that mean first passage path lengths (instead of times) are computed. This is achieved by overriding the mean waiting times to instead represent the mean numbers of steps to exit, which are initially equal to unity for all nodes. Is used in conjunction with **BRANCHPROBS**, in which case each transition represents a move to a different node.

## Other optional keywords

**MAXIT** _int_  
  default is inf. The maximum number of iterations of the relevant algorithm to run before the simulation is terminated (if the target number of A<-B paths to simulate is not reached). The interpretation of this option depends on the chosen enhanced sampling method. e.g. with **WRAPPER WE**, **MAXIT** is the number of iterations of the resampling procedure. With **WRAPPER NONE** and **TRAJ KPS** or **TRAJ MCAMC**, **MAXIT** is the number of basin escape trajectories simulated.

**NTHREADS** _int_  
  number of threads to use in parallel calculations. Defaults to max. no. of threads available.

**BRANCHPROBS**  
  indicates that the transition probabilities in **TRAJ BKL** or **TRAJ KPS** are branching probabilities (ie continuous-time, and with non-uniform waiting times for nodes) calculated from the transition rates. Not compatible with **TRAJ MCAMC**. If **BRANCHPROBS** is not specified, then the default setup is that the linearised transition probability matrix at lag time **TAU** will be calculated (ie continuous-time, with uniform waiting times for all nodes). Alternatively, transition probabilities can be read in with the **TRANSNPROBS** keyword.

**TRANSNPROBS**  
  the edge weights read in from the file *ts_weights.dat* are transition probabilities at lag time **TAU** (ie the waiting times of all nodes are uniform), not rates. Note that this keyword is compatible with all **TRAJ** options. The self-loop transition probabilities for nodes are inferred to be the remainders from unity for transitions from each node. The transition probability matrix is taken to be linearised (ie continuous-time, with uniform waiting times [equal to **TAU**] for all nodes) by default. To simulate a discrete-time Markov chain, specify the **DISCRETETIME** keyword.

**DISCRETETIME**  
  when set with the **TRANSNPROBS** keyword, the edge weights read from *ts_weights.dat* are taken to be the transition probabilities of a discrete-time Markov chain (by default, the transition probabilities are assumed to be continuous-time). **TAU** is then the constant lag time (ie all transitions are associated with time **TAU**, rather than sampling from an exponential distribution). Note that this keyword is compatible with all **TRAJ** options.

**TINTVL** _double_  
  time interval for dumping trajectory information. Negative value (default) indicates that trajectory data is not written (i.e. files _walker.0.y.dat_ are not output). Zero value specifies that all trajectory information is written. An explicit non-negative value must be set if **WRAPPER DIMREDN**. The exact value of **TINTVL** is ignored if **TRAJ KPS** (in which case trajectory data is written after every basin escape).

**DUMPINTVLS**  
  if set, then trajectory data is dumped at precisely the time intervals specified by **TINTVL** (which must therefore be >0.). Otherwise, trajectory data written when the next time interval is exceeded is precisely for the current time of the walker. Path probability and entropy flow are not written in the walker files if this keyword is set, but are still dumped to the *tp_distribns.dat* file. This keyword is required with **WRAPPER DIMREDN**, since the trajectory information required to construct a coarse-grained Markov chain is otherwise not printed.

**SEED** _int_  
  seed for the random number generators (default 19).

**DEBUG**  
  enable extra printing and tests to aid debugging.

**DUMPWAITTIMES**  
  dump the mean waiting times for nodes to the file _meanwaitingtimes.dat_
