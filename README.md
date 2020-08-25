# DISCOTRESS

Welcome to the **DISCOTRESS** (DIscrete State COntinuous Time Rare Event Simulation Suite) program, developed by Daniel J. Sharpe.

![Getting from B to A in a Markov chain](https://github.com/danieljsharpe/danieljsharpe/blob/master/discotress_network_annotated.png?raw=true)

**DISCOTRESS** is a software package to simulate and analyse the dynamics on arbitrary Markov chains. **DISCOTRESS** is designed to enable simulation and analysis of the dynamics for discrete- and continuous-time finite Markov chains (DTMCs and CTMCs, respectively) that exhibit strong metastability (i.e. rare event dynamics). In the rare event regime, standard simulation methods are severely inefficient, and linear algebra methods for the exact computation of dynamical quantities fail to converge.

Discrete-state DTMCs and CTMCs are widely applied models for the stochastic dynamics of many processes. Markov chains are commonly used to represent populations of species in and animal movement within an ecosystem, epidemic spread, financial markets, social groups, gene regulatory and other chemical reaction networks, as well as the dynamics of many-particle systems in condensed matter and biological physics, and more!

## What can I do with DISCOTRESS?

**DISCOTRESS** can be used to:
- simulate the ensemble of A<-B first passage paths (the FPPE) from an initial set of nodes in the Markov chain, denoted B, to an absorbing set of nodes, denoted A [1,3]
- estimate A<-B FPPE statistics by simulation, including the mean first first passage time (MFPT) and FPT distribution, the path probability and entropy flow distributions, and committor probabilities [1,3]
- obtain dynamical quantities characterising the FPPE exactly, including MFPTs, path lengths, committor and absorption probabilities, as well as the expected numbers of node visits and node visitation probabilities, using numerically stable state reduction methods [2,3]
- obtain dynamical quantities characterising the dynamics in the infinite-time limit, namely the stationary distribution and the average mixing time, using numerically stable state reduction methods [2,3]
- estimate and validate a coarse-grained Markov chain constructed from multiple short nonequilibrium trajectories [4]

A host of simulation algorithms for the above purposes are built in to the software, including various state reduction methods, standard kinetic Monte Carlo (BKL), kinetic path sampling (kPS), Monte Carlo with absorbing Markov chains (MCAMC), milestoning, and weighted ensemble (WE) sampling.

## How do I get started?

**Requirements:** C++17  
**Dependencies:** OpenMP

Please see the "documentation" file for a list of keywords that may be included in the input file "input.kmc", a description of the other input/output files and their formats, and the compilation command.

Try running the example in the "tutorials" folder to perform your first kPS simulation!

More tutorials are forthcoming.

## Citations

If you use the DISCOTRESS software in your publication, please cite the following articles:
- [1] D. J. Sharpe and D. J. Wales, _Efficient and exact sampling of transition path ensembles on Markovian networks_, J. Chem. Phys. 153, 024121.
- [2] D. J. Sharpe and D. J. Wales, _Numerical analysis of first passage processes in metastable Markov chains_, (in preparation).
- [3] D. J. Sharpe and D. J. Wales, _Pathways and dynamical observables from absorbing Markov chains_, (in preparation).

Please cite relevant articles describing particular functionality of DISCOTRESS if you use these features:
- [4] D. J. Sharpe and D. J. Wales, _Dimensionality reduction of Markov chains using efficient dynamical simulations_, (in preparation).
- [5] D. Kannan\*, D. J. Sharpe\*, T. D. Swinburne and D. J. Wales, _Optimal dimensionality reduction of Markov chains using graph transformation_, (submitted).
- [6] T. D. Swinburne, D. Kannan, D. J. Sharpe and D. J. Wales, _Rare events and first passage time statistics from the energy landscape_, J. Chem. Phys. (accepted).
- [7] D. Kannan, D. J. Sharpe, T. D. Swinburne and D. J. Wales, _Coarse-graining continuous-time Markov chains with graph transformation_, (in preparation).

More publications for DISCOTRESS are forthcoming.

For references of the various simulation and state reduction algorithms implemented in DISCOTRESS please see the individual .cpp files in the source code.

## Technical FAQs

### 

### Which sampling method should I use for my simulation?

DISCOTRESS divides its dynamical simulation methods into two classes; "wrapper" methods handle an ensemble of so-called walkers (trajectories) via a division of the state space (e.g. WE sampling, milestoning, or no special method), and "trajectory" methods deal with propagating an individual trajectory (BKL, KPS, MCAMC).

For most Markov chains representing a realistic dynamical process of interest, the standard simulation algorithm (BKL) used alone will most likely be too inefficient to simulate the dynamics, because of "flickering" within long-lived macrostates. It becomes necessary to use methods for the simulation of trajectories that are unaffected by metastability (KPS, MCAMC), and/or to employ an enhanced sampling methodology to handle an ensemble of walkers simulated in parallel. Trajectory segments from the walkers can then be stitched together, with appropriate weighting, to yield complete A<-B paths. To choose the appropriate enhanced sampling algorithm, there are many factors to consider, including the information that is desired from the simulation. For instance, milestoning and an adaptation of WE sampling are used to simulate the equilibrium (steady state) ensemble of A<-B transition paths. The other methods, and the state reduction algorithms, compute dynamical quantities for the nonequilibrium ensemble of A<-B paths (i.e. first hitting problem with respect to the initial probability distribution). The characteristics of the topology and dynamics of the Markov chain will also influence which enhanced sampling method is the best choice. For a thorough discussion concerning the choice of sampling method, see Ref. [1].

In general, kPS (see Refs. [1,3]) provides a highly efficient method to simulate the dynamics for metastable DTMCs and CTMCs.

### How should I obtain a predefined partitioning of the Markov chain for use with the enhanced sampling simulation algorithms?

All of the enhanced sampling methods described in Ref. [1] are based on a partitioning of the network into communities, which must accurately characterise the metastable macrostates of the system. Be aware that the choice of this community structure strongly affects the efficiency of the simulation. The question of how to obtain this partitioning is therefore of critical importance.

Even for small networks, spectral methods (such as the original and robust Perron cluster cluster analysis methods, PCCA and PCCA+, respectively) will fail owing to numerical instability if the Markov chain is metastable. Hence, for such networks, BACE (the Bayesian Agglomerative Clustering Engine) is a favourable alternative, but is not scalable. Unfortunately, many state-of-the-art community detection algorithms are based on optimisation of the modularity objective, and therefore _are liable to misrepresent the dynamics_. Multi-level regularised Markov clustering (MLR-MCL), an implementation of which is available at github.com/danieljsharpe/mlr\_mcl, provides a suitable alternative.

### I want my simulation to go faster!

To achieve efficient dynamical simulations using enhanced sampling methods, the community structure must accurately characterise the metastable sets of nodes. It is therefore often worth exploring different algorithms and varying choices of parameters to obtain a suitable partitioning of the Markov chain. There are various rigorous measures to determine the quality of a partitioning [4]. A common mistake when using the kPS and MCAMC algorithms is to not run a fixed number of standard BKL steps after basin escape iterations.

It is also possible to pre-process the Markov chain to reduce dimensionality and eliminate flickering, while maintaining an accurate representation of the slow dynamics for the original Markov chain, see Refs. [6,7].

## Contact the author

github.com/danieljsharpe
daniel.j.sharpe1995@gmail.com
