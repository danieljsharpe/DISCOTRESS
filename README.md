# DISCOTRESS

Welcome to the **DISCOTRESS** (DIscrete State COntinuous Time Rare Event Simulation Suite) program, developed by Daniel J. Sharpe.

![Getting from B to A in a Markov chain](https://github.com/danieljsharpe/danieljsharpe/blob/master/discotress_network_annotated.png?raw=true)

**DISCOTRESS** is a software package to simulate and analyse the dynamics for arbitrary Markov chains. **DISCOTRESS** is designed to enable simulation and analysis of the dynamics for discrete- and continuous-time finite Markov chains (DTMCs and CTMCs, respectively) that are nearly reducible. Such Markov chains exhibit strong metastability; that is, there exists a comparatively slow (i.e. low probability) process. In this rare event regime, standard simulation algorithms are unfeasible owing to their inefficiency, and linear algebra methods for the exact computation of dynamical quantities fail to converge or lead to a severe propagation of numerical error. The advanced methods implemented in DISCOTRESS circumvent these problems.

Discrete-state DTMCs and CTMCs are widely applied models for the stochastic dynamics of many processes. Markov chains are commonly used to represent populations of species in and animal movement within an ecosystem :parrot: :palm_tree:, epidemic spread :microbe: :mosquito:, financial markets :money_with_wings: :chart_with_upwards_trend:, social groups :family: :couple:, gene regulatory and other chemical reaction networks :dna: :test_tube:, as well as the dynamics of many-particle systems in condensed matter and biological physics :atom_symbol: :petri_dish:, and more!

## What can I do with DISCOTRESS?

**DISCOTRESS** can be used to:
- sample the ensemble of &#120068; &#8592; &#120069; first passage paths (the FPPE) from an initial set of nodes in the Markov chain, denoted &#120069;, to an absorbing set of nodes, denoted &#120068; [1,2].
- simulate probability distributions of path properties for the &#120068; &#8592; &#120069; FPPE, including the first passage time (FPT) [1,2], path probability, and entropy flow distributions of paths.
- simulate statistics associated with nodes (or groups thereof) for the ensemble of direct &#120068; &#8592; &#120069; transition paths (the TPE), namely committor and visitation probabilities [1,2]. These quantities can be calculated for both the nonequilbrium and equilibrium (steady state) TPEs [3].
- obtain dynamical quantities characterising the &#120068; &#8592; &#120069; nonequilibrium and equilibrium FPPE and TPE exactly, including MFPTs, committor and absorption probabilities, expected numbers of node visits, and node visitation probabilities, using numerically stable state reduction methods [2,3].
- obtain dynamical quantities characterising the dynamics in the infinite-time limit, namely the stationary distribution and the average mixing time, using numerically stable state reduction methods [2,3].
- estimate and validate a coarse-grained Markov chain constructed from multiple short nonequilibrium trajectories [4].

A host of simulation algorithms to sample paths are built in to the software [1,2], including standard kinetic Monte Carlo (BKL), kinetic path sampling (kPS), Monte Carlo with absorbing Markov chains (MCAMC), milestoning, and weighted ensemble (WE) sampling.

## How do I get started?

**Requirements:** C++17  
**Dependencies:** OpenMP

Refer to the [documentation](https://github.com/danieljsharpe/DISCOTRESS/blob/master/documentation.md) for a list of keywords that may be included in the file *input.kmc*, a description of the other input and output files and their formats, and the compilation command.

Try running the [tutorial examples](https://github.com/danieljsharpe/DISCOTRESS_tutorials) to get started!

Use the available [analysis scripts](https://github.com/danieljsharpe/DISCOTRESS_tools) to find out what it all means.

## Citations

If you use the DISCOTRESS software in your publication, please cite the following articles:
- [1] D. J. Sharpe and D. J. Wales, [Efficient and exact sampling of transition path ensembles on Markovian networks](https://doi.org/10.1063/5.0012128), J. Chem. Phys. 153, 024121.
- [2] D. J. Sharpe and D. J. Wales, Pathways and dynamical observables in metastable Markov chains, (in preparation).

Please cite relevant articles describing particular functionality of DISCOTRESS if you use these features:
- [3] D. J. Sharpe and D. J. Wales, Numerical analysis of first passage processes in metastable Markov chains, (in preparation).
- [4] D. J. Sharpe and D. J. Wales, Dimensionality reduction of Markov chains using efficient dynamical simulations, (in preparation).
- [5] D. Kannan\*, D. J. Sharpe\*, T. D. Swinburne and D. J. Wales, Optimal dimensionality reduction of Markov chains using graph transformation, (submitted).
- [6] D. Kannan, D. J. Sharpe, T. D. Swinburne and D. J. Wales, Coarse-graining continuous-time Markov chains with graph transformation, (in preparation).

Further example applications can be found in the following publications:
- [7] T. D. Swinburne, D. Kannan, D. J. Sharpe and D. J. Wales, [Rare events and first passage time statistics from the energy landscape](https://doi.org/10.1063/5.0016244), J. Chem. Phys. 153, 134115.

More publications for DISCOTRESS are forthcoming.

For additional citations related to the various simulation and state reduction algorithms implemented in DISCOTRESS, please see the individual .cpp files in the source code.

## Technical FAQs

### 

### Which sampling method should I use for my simulation?

DISCOTRESS divides its dynamical simulation methods into two classes; "wrapper" methods handle an ensemble of so-called walkers (independent trajectories) via a division of the state space (e.g. WE sampling, milestoning, or no special method), and "trajectory" methods deal with propagating an individual trajectory (BKL, KPS, MCAMC).

For Markov chains that are nearly reducible, which commonly arise when constructing a network that represents a realistic dynamical process of interest, the standard simulation algorithm (BKL) used alone will most likely be too inefficient to simulate the dynamics, because of "flickering" within long-lived (i.e. metastable) macrostates. It becomes necessary to use methods for the simulation of trajectories that are unaffected by metastability (KPS, MCAMC), and/or to employ an enhanced sampling methodology to handle an ensemble of walkers simulated in parallel. Trajectory segments from the walkers can then be stitched together, with appropriate weighting, to yield complete &#120068; &#8592; &#120069; paths. To choose the appropriate enhanced sampling algorithm, there are many factors to consider, including the information that is desired from the simulation. For instance, milestoning and an adaptation of WE sampling are used to simulate the equilibrium (steady state) ensemble of &#120068; &#8592; &#120069; transition paths [1], which can also be achieved by simulating a very long trajectory that continually transitions between the &#120068; and &#120069; states [2,3]. The other simulation methods compute dynamical quantities for the nonequilibrium ensemble of &#120068; &#8592; &#120069; paths (i.e. first hitting problem with respect to the initial probability distribution). The characteristics of the topology and dynamics of the Markov chain will also influence which enhanced sampling method is the best choice. For a thorough discussion concerning the choice of sampling method, see [1].

In general, kPS [1,2] provides a highly efficient method to simulate the dynamics for metastable DTMCs and CTMCs.

The state reduction algorithms for the exact computation of dynamical properties also determine quantities associated with the nonequilibrium TPE, although equilibrium TPE properties can be subsequently derived [3]. The state reduction algorithms can be used to characterise the &#120068; &#8592; &#120069; process at both a microscopic and macroscopic level of detail, but are typically feasible only for Markov chains of no more than several thousand nodes. Simulation offers a more scalable approach to estimate the same properties via sampling of trajectories.

### How should I obtain a predefined partitioning of the Markov chain for use with the enhanced sampling simulation algorithms?

All of the enhanced sampling methods described in [1] are based on a partitioning of the network into communities, which must accurately characterise the metastable macrostates of the system. Be aware that the choice of this community structure strongly affects the efficiency of the simulation. The question of how to obtain this partitioning is therefore of critical importance.

Even for small networks, spectral methods (such as the original and robust Perron cluster cluster analysis methods, PCCA and PCCA+, respectively) will fail owing to numerical instability if the Markov chain is metastable. Hence, for such networks, BACE (the Bayesian Agglomerative Clustering Engine) is a favourable alternative, but is not scalable. More scalable community detection algorithms are usually stochastic, but such algorithms may misclassify nodes close to the boundary between metastable macrostates. This problem can be attenuated by variational refinement of an initial clustering [4].  Unfortunately, many state-of-the-art community detection algorithms are based on optimisation of the modularity objective, and therefore _are liable to misrepresent the dynamics_. Multi-level regularised Markov clustering (MLR-MCL), an implementation of which is [here](https://github.com/danieljsharpe/mlr_mcl), provides a suitable alternative [1].

### I want my simulation to go faster!

To achieve efficient dynamical simulations using enhanced sampling methods, the community structure must accurately characterise the metastable sets of nodes. It is therefore often worth exploring different algorithms and varying choices of parameters to obtain a suitable partitioning of the Markov chain. There are various rigorous measures to determine the quality of a partitioning [4]. A common mistake when using the kPS and MCAMC algorithms is to not run a fixed number of standard BKL steps after basin escape iterations.

It is also possible to pre-process the Markov chain to reduce dimensionality and eliminate flickering, while maintaining an accurate representation of the slow dynamics for the original Markov chain. A simple method for this purpose is to recursively regroup states according to a specified transition rate or probability threshold [7]. A more complex approach is to renormalise the Markov chain following the elimination of appropriately chosen nodes, which has the advantage of preserving the mean and variance of the &#120068; &#8592; &#120069; FPT distribution if only nodes outside of the &#120068; and &#120069; sets are eliminated [6].

## Contact the author

Would you like to request a new feature in DISCOTRESS? Or simply have a question? Get in touch:

daniel.j.sharpe1995@gmail.com
