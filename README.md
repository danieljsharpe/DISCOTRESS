# DISCOTRESS

Welcome to the **DISCOTRESS** (DIscrete State COntinuous Time Rare Event Simulation Suite) program, developed by Daniel J. Sharpe.

![Getting from B to A in a Markov chain](https://github.com/danieljsharpe/danieljsharpe/blob/master/discotress_network_annotated.png?raw=true)

**DISCOTRESS** is a software package to simulate and analyse the dynamics for arbitrary Markov chains. **DISCOTRESS** is designed to enable simulation and analysis of the dynamics for discrete- and continuous-time finite Markov chains (DTMCs and CTMCs, respectively) that are nearly reducible. Such Markov chains exhibit strong metastability; that is, there exists a comparatively slow (i.e. low probability) process. In this rare event regime, which is frequently encountered in realistic modeling tasks, standard simulation algorithms are unfeasible owing to their inefficiency, and linear algebra methods for the exact computation of dynamical quantities fail to converge or lead to a severe propagation of numerical error. The advanced methods implemented in DISCOTRESS circumvent these problems.

Discrete-state DTMCs and CTMCs are widely applied models for the stochastic dynamics of many processes. Markov chains are commonly used to represent populations of species in and animal movement within an ecosystem :parrot: :palm_tree:, epidemic spread :microbe: :mosquito:, financial markets :money_with_wings: :chart_with_upwards_trend:, climate dynamics :sun_behind_rain_cloud: :tornado:, gene regulatory and other chemical reaction networks :dna: :test_tube:, as well as the dynamics of many-particle systems in condensed matter and biological physics :atom_symbol: :petri_dish:, and more!

## What can I do with DISCOTRESS?

A host of simulation algorithms to sample paths are built in to **DISCOTRESS** [1,2], including standard kinetic Monte Carlo (BKL), kinetic path sampling (kPS), Monte Carlo with absorbing Markov chains (MCAMC), milestoning, and weighted ensemble (WE) sampling. The software also includes algorithms to perform exact numerical analysis. The flexibility of the software and the relative advantages of the different algorithms allow for a variety of problems to be treated. **DISCOTRESS** can be used to:

- sample the ensemble of &#120068; &#8592; &#120069; first passage paths (the FPPE) from an initial set of nodes in the Markov chain, denoted &#120069;, to an absorbing (target) set of nodes, denoted &#120068; [1,2].
- simulate probability distributions of path properties for the &#120068; &#8592; &#120069; FPPE, including the first passage time (FPT) [1,2], path probability [4], and entropy flow distributions of paths.
- simulate statistics associated with nodes (or groups thereof) for the ensemble of direct &#120068; &#8592; &#120069; transition paths (the TPE), namely committor and visitation probabilities [1,2]. These quantities can be calculated for both the nonequilbrium and equilibrium (steady state) TPEs [3].
- obtain dynamical quantities characterising the &#120068; &#8592; &#120069; nonequilibrium and equilibrium FPPE and TPE exactly, including MFPTs, committor and absorption probabilities, expected numbers of node visits, and node visitation probabilities, using numerically stable state reduction methods [2,3].
- obtain dynamical quantities characterising the dynamics in the infinite-time limit, namely the stationary distribution and the average mixing time, using numerically stable state reduction methods [2,3].
- determine the set of &#120068; &#8592; &#120069; first passage paths with the highest probabilities, using a *k* shortest paths algorithm [4].
- estimate and validate a coarse-grained Markov chain constructed from multiple short nonequilibrium trajectories [5].

## How do I get started?

**Requirements:** C++17  
**Dependencies:** OpenMP

Refer to the [documentation](https://github.com/danieljsharpe/DISCOTRESS/blob/master/documentation.md) for a list of keywords that may be included in the file *input.kmc*, a description of the other input and output files and their formats, and the compilation command.

Try running the [tutorial examples](https://github.com/danieljsharpe/DISCOTRESS_tutorials) to get started!

Use the available [analysis scripts](https://github.com/danieljsharpe/DISCOTRESS_tools) to find out what it all means.

Read the [FAQs](https://github.com/danieljsharpe/DISCOTRESS/blob/master/documentation.md) for hints and tips.

## Citations

If you use the DISCOTRESS software in your publication, please cite the following articles:
- [1] D. J. Sharpe and D. J. Wales, [Efficient and exact sampling of transition path ensembles on Markovian networks](https://doi.org/10.1063/5.0012128), J. Chem. Phys. 153, 024121.
- [2] D. J. Sharpe and D. J. Wales, Pathways and dynamical observables in metastable Markov chains, (in preparation).

Please cite relevant articles describing particular functionality of DISCOTRESS if you use these features:
- [3] D. J. Sharpe and D. J. Wales, Numerical analysis of first passage processes in metastable finite Markov chains, (in preparation).
- [4] D. J. Sharpe and D. J. Wales, Graph transformation algorithm for the expectation of first passage path properties in finite Markov chains, (in preparation).
- [5] D. J. Sharpe and D. J. Wales, Dimensionality reduction of Markov chains using efficient dynamical simulations, (in preparation).
- [6] D. Kannan\*, D. J. Sharpe\*, T. D. Swinburne and D. J. Wales, Optimal dimensionality reduction of Markov chains using graph transformation, (in press).
- [7] D. Kannan, D. J. Sharpe, T. D. Swinburne and D. J. Wales, Coarse-graining continuous-time Markov chains with graph transformation, (in preparation).

Further example applications can be found in the following publications:
- [8] T. D. Swinburne, D. Kannan, D. J. Sharpe and D. J. Wales, [Rare events and first passage time statistics from the energy landscape](https://doi.org/10.1063/5.0016244), J. Chem. Phys. 153, 134115.

More publications for DISCOTRESS are forthcoming.

For additional citations related to the various simulation and state reduction algorithms implemented in DISCOTRESS, please see the individual .cpp files in the source code.

## Contact the author

Would you like to request a new feature in DISCOTRESS? Or simply have a question? Get in touch:

daniel.j.sharpe1995@gmail.com
