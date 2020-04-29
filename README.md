# DISCOTRESS

Welcome to the **DISCOTRESS** (DIscrete State COntinuous Time Rare Event Simulation Suite) program, developed by Daniel J. Sharpe.

**DISCOTRESS** is a software package to simulate the dynamics on arbitrary continuous time Markov chains (CTMCs). **DISCOTRESS** is designed to enable simulation of the dynamics even for CTMCs that exhibit strong metastability (i.e. rare event dynamics), where standard simulation methods fail.

CTMCs are widely applied models for the stochastic dynamics of many processes. CTMCs are commonly used to represent populations of species in an ecosystem, epidemic spread, financial markets, and gene regulatory networks, as well as the dynamics of many-particle systems in condensed matter and biological physics.

## What can I do with DISCOTRESS?

**DISCOTRESS** can be used to:
- simulate the ensemble of A<-B transition paths (the TPE) from an initial set of nodes in the CTMC, denoted B, to an absorbing set of nodes, denoted A
- obtain A<-B TPE statistics including the mean first first passage time (MFPT) and FPT distribution, the path probability and entropy flow distributions, and committor probabilities
- estimate and validate a coarse-grained Markov chain constructed from multiple short nonequilibrium trajectories

A host of simulation algorithms for the above purposes are built in to the software, including standard kinetic Monte Carlo (BKL), kinetic path sampling (kPS), Monte Carlo with absorbing Markov chains (MCAMC), milestoning, and weighted ensemble (WE) sampling.

## How do I get started?

**Requirements:** C++17
**Dependencies:** None!

Please see the "documentation" file for a list of keywords that may be included in the input file "input.kmc", a description of the other input/output files and their formats, and the compilation command.

Try running the example in the "tutorials" folder to perform your first kPS simulation!

Tutorials are forthcoming.

## Citations

If you use the DISCOTRESS software in your publication, please cite the following articles:
- [1] D. J. Sharpe and D. J. Wales, _Efficient and exact sampling of transition path ensembles on Markovian networks_, J. Chem. Phys. (submitted).

More publications for DISCOTRESS are forthcoming.

For references of the individual simulation algorithms please see the individual .cpp files in the source code.

## Technical FAQs

### Which sampling method should I use for my simulation?

For most CTMCs representing a realistic dynamical process of interest, the standard simulation algorithm (BKL) will most likely be too inefficient to simulate the dynamics. To choose the appropriate enhanced sampling algorithm, there are many factors to consider, including the information that is desired from the simulation. The characteristics of the topology and dynamics of the CTMC will also influence which enhanced sampling method is optimal. For a thorough discussion concerning the choice of sampling method, see Ref. [1].

### How should I obtain a predefined partitioning of the Markov chain for use with the enhanced sampling simulation algorithms?

The partitioning of the network into communities must accurately characterise the metastable states of the system. Even for small networks, spectral methods (such as the original and robust Perron cluster cluster analysis methods, PCCA and PCCA+, respectively) will fail owing to numerical instability if the CTMC is metastable. Hence, for such networks, BACE (the Bayesian Agglomerative Clustering Engine) is a favourable alternative, but is not scalable. Unfortunately, many state-of-the-art community detection algorithms are based on optimisation of the modularity objective, and therefore _are liable to misrepresent the dynamics_. Multi-level regularised Markov clustering (MLR-MCL), an implementation of which is available at github.com/danieljsharpe/mlr_mcl, provides a suitable alternative.

## Contact the author

github.com/danieljsharpe
daniel.j.sharpe1995@gmail.com
