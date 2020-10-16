# DISCOTRESS user information

## Compilation and execution

Compile with:
```bash
g++ -std=c++17 discotress.cpp kmc_methods.cpp we_kmc.cpp ffs_kmc.cpp neus_kmc.cpp milestoning.cpp kps.cpp mcamc.cpp keywords.cpp ktn.cpp -o discotress -fopenmp
```

To run the code, simply type the magic word: `discotress`

Get started with the [tutorials](https://github.com/danieljsharpe/DISCOTRESS_tutorials).

## Input files

The instructions to the program are set by the keywords in the _input.kmc_ file. 

Additional mandatory input files specify the network topology and parameters:

File | Format | Description
---- | ------ | -----------
*stat\_prob.dat* | single col, length **NNODES** | natural log of stationary probabilities for nodes
*edge\_conns.dat* | double col, length **NEDGES** | each line specifies a bidirectional edge connecting nodes with the stated IDs (indexed from 1)
*edge\_weights.dat* | double col, length **NEDGES** | transition probability (for DTMC) or natural log of transition rate (for CTMC) associated with forward and reverse transitions (consistent with *edge\_conns.dat*)
*nodes.A* | single col, length set by **NODESAFILE** | list of IDs for nodes in the absorbing set _A_ (indexed from 1)
*nodes.B* | single col, length set by **NODESBFILE** | list of IDs for nodes in the initial set _B_ (indexed from 1)

Additional optional input files, corresponding to particular optional keywords, relate to the simulation parameters and to the recording of some statistics:

File | Format | Keyword | Description
---- | ------ | ------- | -----------
*communities.dat* | single col, length **NNODES** | **COMMSFILE** | community IDs for nodes (indexed from 0), specifies network partitioning for use in enhanced sampling methods
*bins.dat* | single col, length **NNODES** | **BINSFILE** | bin IDs for nodes (indexed from 0), used for calculating transition path statistics (defaults to *communities.dat*)
*initcond.dat* | single col, length cf. **NODESBFILE** | **INITCONDFILE** | initial occupation probability distribution for nodes of the initial set _B_. Defaults to be proportional to a local equilibrium within _B_
*ntrajs.dat* | single col, length cf. **COMMSFILE** | **DIMREDUCTION** | numbers of short trajectories to run, initialised from each community

## Output files

The following output files are written as the result of a computation to simulate trajectories. For a description of the output files when instead instructing the program to perform an exact computation, see the documentation for the keywords in the *input.kmc* file that relate to state reduction procedures.

File | Description | Format
---- | ----------- | ------
*tp\_distribns.dat* | properties of _A_ &#8592 _B_ paths, from which the probability distributions for first passage  can be analysed | path no. / path length / path time / ln of path probability / path entropy flow
*tp\_stats.dat* | bin statistics, written if communities were specified | bin ID / no. of reactive (direct _A_ &#8592 _B_) paths for which bin is visited / no. of paths for which bin is visited and trajectory returned to initial set _B_ / reactive visitation probability / committor probability
*walker.x.y.dat* | trajectory information dumped at the specified time intervals (or when a trajectory escapes from a community, depending on options). *x* is the walker ID, *y* is the path number | node ID / community ID / path length / path time / ln of path probability / path entropy flow

## Main keywords

The following is a list of keywords that must always appear in the _input.kmc_ file.

----

**NNODES** `int`  
  mandatory, number of nodes in the Markov chain.

**NEDGES** `int`   
  mandatory, number of (bidirectional) edges in the Markov chain.

**WRAPPER** `str`  
  mandatory, employ a \`wrapping' enhanced sampling strategy for accelerating the observation of rare events. Options:  
-    **NONE**    \- no enhanced sampling strategy employed
-    **DIMREDN** \- special class to simulate many short nonequilibrium trajectories starting from each community in turn, used for estimation of a coarse-grained Markov chain. See **DIMREDUCTION** keyword for more detail 
-    **WE**      \- weighted ensemble sampling
-    **FFS**     \- forward flux sampling
-    **NEUS**    \- non-equilibrium umbrella sampling
-    **MILES**   \- milestoning

**TRAJ** `str`  
  mandatory, method for propagating individual trajectories. Options:  
-    **BKL**     \- Bortz-Kalos-Lebowitz \(aka n-fold way\), i.e. standard kinetic Monte Carlo (kMC)
-    **KPS**     \- kinetic path sampling algorithm
-    **MCAMC**   \- Monte Carlo with absorbing Markov chains algorithm  
  Note the following excluded combinations of WRAPPER / TRAJ options: WE / KPS, WE / MCAMC, NEUS / KPS, NEUS / MCAMC, DIMREDN / BKL.

**NODESAFILE** `str` `int`  
  mandatory if not **WRAPPER DIMREDN**, name of the file containing the node IDs (indexed from 1) belonging to the _A_ (absorbing) set, and number of nodes in the _A_ set.

**NODESBFILE** `str` `int`  
  mandatory if not **WRAPPER DIMREDN**, name of the file containing the node IDs (indexed from 1) belonging to the _B_ (initial) set, and number of nodes in the _B_ set.

----

## Optional keywords relating to simulation parameters and output

The following is a list of keywords that specify generic simulation parameters and keywords that control the recording of statistics or other aspects of the output.

----

**BINSFILE** `str` `int`
  optional. Default to be the same as **COMMSFILE**, if specified.
  Name of the file containing the bin IDs (indexed from 0) for nodes, and number of bins. The bins are used to collect statistics associated with nodes (or groups thereof) for the _A_ &#8592 _B_ transition path ensemble, namely committor and visitation probabilities.

**COMMSFILE** `str` `int`
  mandatory if **WRAPPER** is **DIMREDN**, **FFS**, **NEUS**, or **MILES**. Also mandatory if **WRAPPER** is **WE**, or if **TRAJ** is **KPS** or **MCAMC**, and **ADAPTIVECOMMS** is not specified.
  Name of the file containing the definitions of communities (single-column, indexed from zero, number of entries equal to the number of nodes **NNODES** in the network) and no. of communities. Is overridden by **ADAPTIVECOMMS**. For both **WRAPPER** and **TRAJ** enhanced sampling methods, except **TRAJ BKL**, the communities are used to divide the state space (eg the communities define the trapping basins in **KPS**, or the communities for resampling in **WE**), and for certain algorithms may dictate the resolution at which the transition path statistics (see **BINSFILE** keyword) can be calculated. The specification of communities must be consistent with the definition of the _A_ and _B_ sets. An exception is if the number of communities is 2, in which case the initial set _B_ can be a subset of the relevant community. Note that if this is chosen to be the case, then re-hitting _B_ is not detected, and committor and transient visitation probabilities for the bins will be incorrect.

**DUMPINTVLS**
  if set, then trajectory data is dumped at precisely the time intervals specified by **TINTVL** (which must therefore be >0.). Otherwise, trajectory data written when the next time interval is exceeded is precisely for the current time of the walker. Path probability and entropy flow are not written in the walker files if this keyword is set, but are still dumped to the *tp_distribns.dat* file. This keyword is required with **WRAPPER DIMREDN**, since the trajectory information required to construct a coarse-grained Markov chain is otherwise not printed.

**INITCONDFILE** `str`
  optional. Name of the file containing initial occupation probabilities for nodes in _B_ that are alternative to the stationary probabilities. The number of entries is assumed to be the same as the specified number of nodes in _B_, and the specified order is assumed to be the same also \(cf **NODESBFILE**\). The values must sum to unity.

**MAXIT** `int`
  default is inf. The maximum number of iterations of the relevant algorithm to run before the simulation is terminated (if the target number of _A_ &#8592 _B_ paths to simulate is not reached). The interpretation of this option depends on the chosen enhanced sampling method. e.g. with **WRAPPER WE**, **MAXIT** is the number of iterations of the resampling procedure. With **WRAPPER NONE** and **TRAJ KPS** or **TRAJ MCAMC**, **MAXIT** is the number of basin escape trajectories simulated.

**NABPATHS** `int`
  mandatory if not **WRAPPER DIMREDN** and if none of the state reduction keywords are specified. The simulation is terminated when this number of _A_ &#8592 _B_ paths have been successfully sampled.

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

**DIMREDUCTION** `str` `long double`
  mandatory if **WRAPPER DIMREDN**, which initialises a special wrapper class that does not perform the usual code function, which is to simulate _A_ &#8592 _B_ transition paths, and instead instructs the program to simulate many short trajectories starting from each community, each of length in time equal to the float argument. These trajectories are printed to files _walker.x.y.dat_, where _x_ is the ID of the community, and _y_ is the iteration number for that community. Trajectory information is written to files whenever a trajectory transitions to a new community. The total number of trajectories that are to be simulated starting from each community is listed in the file given as the string arg. This calculation is parallelised, using a number of threads equal to **NTHREADS** (defaults to maximum number of threads available). This calculation is compatible with two algorithms to propagate individual trajectories: **TRAJ KPS**, and **TRAJ MCAMC** (without **MEANRATE**). The communities of nodes must be specified (**COMMSFILE** keyword). **NABPATHS**, **MAXIT**, and **BINSFILE** keywords are ignored. This setup is incompatible with specification of an initial condition via the **INITCONDFILE** keyword, and with the **NODESAFILE** and **NODESBFILE** keywords. Instead, a local equilibrium within the starting community is assumed as the initial probability distribution for each macrostate. Note that **WRAPPER DIMREDN** requires both this keyword and **DUMPINTVLS** to be set.

**KPSKMCSTEPS** `int`  
  optional. If **TRAJ** is **KPS** or **MCAMC**, specifies the number of standard BKL steps to be performed after a kPS or MCAMC escape from a trapping basin. Default is 0 (pure kPS (or MCAMC), no kMC steps). However, this is not the recommended value. If using **TRAJ KPS** or **TRAJ MCAMC**, for most systems, great gains in simulation efficiency will be achieved by setting **KPSKMCSTEPS** to an appropriate nonzero value. This is because many metastable systems will feature transition regions between metastable states. Therefore, after each basin escape, the trajectory will likely flicker between the two basins. Rather than simulate expensive kPS or MCAMC basin escape iterations for these trivial recrossings, it is much more efficient to perform standard BKL steps. Note that this keyword does not require **BRANCHPROBS** to be set, and can also be used with **DISCRETETIME**. Ignored if **ADAPTIVECOMMS**.

**MEANRATE**  
  optional. If **TRAJ MCAMC**, the calculation uses the approximate mean rate method, as opposed to the default exact first passage time analysis (FPTA) method.

**NELIM** `int`  
  mandatory if **TRAJ KPS**. The maximum number of nodes that are to be eliminated from the current trapping basin. If **NELIM** exceeds the number of nodes in the largest community, then all states of any trapping basin are always eliminated. Note that **NELIM** determines the number of transition matrices stored for the active subnetwork, and therefore the choice of this keyword (along with the sizes of communities) can strongly affect memory usage.

**TAURE** `double`
  mandatory if **WRAPPER WE**, the time between resampling trajectories

----

## Optional keywords related to exact numerical analysis of the dynamics by state reduction methods

The following keywords are used in combination with the keywords **WRAPPER NONE** and **TRAJ KPS**. Use of any of the following keywords overrides the default functionality of DISCOTRESS, which is to simulate dynamical paths, and instead instructs the program to perform a state reduction procedure to exactly compute one or more dynamical quantities associated with nodes, in a numerically stable manner. **NABPATHS** must be set to some arbitrary number >0. The _communities.dat_ file must specify precisely two communities; namely, nodes in the target set _A_ and nodes not in _A_. The **COMMITTOR**, **ABSORPTION**, **MFPT**, and **GTH** keywords can be used together in any combination. The computations performed with the **FUNDAMENTALRED** and **FUNDAMENTALIRRED** keywords are standalone operations.

The memory costs of state reduction computations can be reduced without consequence by setting the type of the `h` members (which represent the numbers of kMC steps for transitions, when using the **KPS** algorithm) of the `Node` and `Edge` structures (defined in the file *ktn.h*) to `int`, since these members are not used in the state reduction methods.

----

**ABSORPTION**
  specifies that a state reduction procedure is performed to compute the absorption probabilities. The probabilities b\_ij that a trajectory initialised from the non-absorbing node i is absorbed at node j are written to the file *absorption.dat* in the format "i / j / b\_ij". For the initial occupation probability distribution (which, by default, is assumed to be a local equilibrium within the initial set _B_), the absorption (hitting) probabilities for each absorbing node are printed to the file *hitting\_probs.dat*.

**COMMITTOR**  
  specifies that a state reduction procedure is performed to compute the committor probabilities for all nodes. The committor probabilities are written to the files *committor\_AB.dat* and *committor\_BA.dat* (for _A_ &#8592 _B_ and _B_ &#8592 _A_ directions, respectively). Note that the committor probabilities determined by this method are for each node, and the calculation is exact and deterministic (unlike calculation of the committor probabilities for the bins from simulation data, cf. the **BINSFILE** keyword).

**FUNDAMENTALIRRED**  
  specifies that a state reduction algorithm is used to compute the fundamental matrix of an irreducible Markov chain. The (i,j) edge weights z\_ij of the renormalised network resulting from this procedure are the elements of the fundamental matrix, the trace of which gives the Kemeny constant (average mixing time) for the Markov chain. These values are written to the file *fundamental.dat* in the format "i / j / z\_ij", and the Kemeny constant is printed in the output.

**FUNDAMENTALRED**
  specifies that a state reduction algorithm is used to compute the fundamental matrix of an absorbing (i.e. reducible) Markov chain. The (i,j) edge weights n\_ij of the renormalised network resulting from this procedure are the expected numbers of times that the j-th node is visited along first passage paths initialised from the i-th node. These values are written to the file *transient\_visits.dat* in the format "i / j / n\_ij". The node visitation probabilities can be computed from this information if the committor probabilities are also known. For the initial occupation probability distribution (which, by default, is assumed to be a local equilibrium within the initial set _B_), the expected numbers of times that non-absorbing nodes are visited along first passage paths are printed to the file *node\_visits.dat*.

**GTH**
  specifies that the stationary probability distribution (which exists if the Markov chain is irreducible) is computed using the Grassmann-Taksar-Heyman (GTH) algorithm. Can only be used when the target set _A_ contains a single node. The input file *stat\_prob.dat* must be provided, but its contents are not used. The stationary probabilities determined by the GTH algorithm are written to the file *stat\_prob\_gth.dat*.

**MFPT**  
  specifies that a state reduction procedure is performed to compute mean first passage times (MFPT). The MFPTs m\_iA for transitions from non-absorbing nodes i to the set of absorbing nodes _A_ are written to the file *mfpt.dat* in the format "i / m\_iA". Given an initial occupation probability distribution (which, by default, is assumed to be a local equilibrium within the initial set _B_), the _A_ &#8592 _B_ MFPT is printed in the output. If the initial mean waiting times of nodes are set to the initial mean number of steps to exit (i.e. equal to unity for all nodes), then the MFPTs are in fact the mean first passage path lengths.

**PATHLENGTHS**  
  when used in conjunction with **MFPT**, specifies that mean first passage path lengths (instead of times) are computed. This is achieved by overriding the mean waiting times to instead represent the mean numbers of steps to exit, which are initially equal to unity for all nodes. Is used in conjunction with **BRANCHPROBS**, in which case each transition represents a move to a different node.

----

## Other optional keywords

----

**BRANCHPROBS**  
  indicates that the transition probabilities in **TRAJ BKL** or **TRAJ KPS** are branching probabilities (ie continuous-time, and with non-uniform waiting times for nodes) calculated from the transition rates. Not compatible with **TRAJ MCAMC**. If **BRANCHPROBS** is not specified, then the default setup is that the linearised transition probability matrix at lag time **TAU** will be calculated (ie continuous-time, with uniform waiting times for all nodes). Alternatively, transition probabilities can be read in with the **TRANSNPROBS** keyword.

**DEBUG**
  enable extra printing and tests to aid debugging.

**DISCRETETIME**
  when set with the **TRANSNPROBS** keyword, the edge weights read from *edge\_weights.dat* are taken to be the transition probabilities of a discrete-time Markov chain (by default, the transition probabilities are assumed to be continuous-time). **TAU** is then the constant lag time (ie all transitions are associated with time **TAU**, rather than sampling from an exponential distribution). Note that this keyword is compatible with all **TRAJ** options.

**DUMPWAITTIMES**
  dump the mean waiting times for nodes to the file _meanwaitingtimes.dat_

**NTHREADS** `int`
  number of threads to use in parallel calculations. Defaults to max. no. of threads available.

**SEED** `int`  
  seed for the random number generators (default 19).

**TAU** `long double`
  mandatory if not **BRANCHPROBS**.
  **TAU** is the lag time at which the linear transition probability matrices are estimated (or provided, if **TRANSNPROBS**). If **TRANSPROBS** and **DISCRETETIME** are specified, then **TAU** is the lag time of the discrete-time Markov chain.

**TRANSNPROBS**
  the edge weights read in from the file *edge\_weights.dat* are transition probabilities at lag time **TAU** (ie the waiting times of all nodes are uniform), not rates. Note that this keyword is compatible with all **TRAJ** options. The self-loop transition probabilities for nodes are inferred to be the remainders from unity for transitions from each node. The transition probability matrix is taken to be linearised (ie continuous-time, with uniform waiting times [equal to **TAU**] for all nodes) by default. To simulate a discrete-time Markov chain, specify the **DISCRETETIME** keyword.

----

