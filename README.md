# DISCOTRESS

Welcome to the **DISCOTRESS** (DIscrete State COntinuous Time Rare Event Simulation Suite) program, developed by Daniel J. Sharpe.

**DISCOTRESS** is a software package to simulate the dynamics on arbitrary Markov chains. **DISCOTRESS** is designed to enable simulation of the dynamics for continuous- and discrete-time Markov chains (CTMCs and DTMCs, respectively) that exhibit strong metastability (i.e. rare event dynamics), where standard simulation methods fail.

CTMCs and DTMCs are widely applied models for the stochastic dynamics of many processes. Markov chains are commonly used to represent populations of species in an ecosystem, epidemic spread, financial markets, social groups, and gene regulatory networks, as well as the dynamics of many-particle systems in condensed matter and biological physics, and more!

## What can I do with DISCOTRESS?

**DISCOTRESS** can be used to:
- simulate the ensemble of A<-B transition paths (the TPE) from an initial set of nodes in the CTMC or DTMC, denoted B, to an absorbing set of nodes, denoted A [1]
- obtain A<-B TPE statistics including the mean first first passage time (MFPT) and FPT distribution, the path probability and entropy flow distributions, and committor probabilities [1]
- estimate and validate a coarse-grained Markov chain constructed from multiple short nonequilibrium trajectories [2]

A host of simulation algorithms for the above purposes are built in to the software, including standard kinetic Monte Carlo (BKL), kinetic path sampling (kPS), Monte Carlo with absorbing Markov chains (MCAMC), milestoning, and weighted ensemble (WE) sampling.

## How do I get started?

**Requirements:** C++17  
**Dependencies:** OpenMP

Please see the "documentation" file for a list of keywords that may be included in the input file "input.kmc", a description of the other input/output files and their formats, and the compilation command.

Try running the example in the "tutorials" folder to perform your first kPS simulation!

Tutorials are forthcoming.

## Citations

If you use the DISCOTRESS software in your publication, please cite the following articles:
- [1] D. J. Sharpe and D. J. Wales, _Efficient and exact sampling of transition path ensembles on Markovian networks_, J. Chem. Phys. (submitted).

Please cite relevant articles describing particular functionality of DISCOTRESS if you use these features:
- [2] D. J. Sharpe and D. J. Wales, _Dimensionality reduction of Markovian transition networks using efficient dynamical simulations_, (in preparation).
- [3] D. Kannan\*, D. J. Sharpe\*, T. D. Swinburne and D. J. Wales, _Computation of inter-community first passage times in Markov chains_, (in preparation).
- [4] T. D. Swinburne, D. Kannan, D. J. Sharpe and D. J. Wales, _Rare events and first passage time statistics from the energy landscape_, J. Chem. Phys. (submitted).
- [5] D. Kannan, D. J. Sharpe, T. D. Swinburne and D. J. Wales, _Coarse-graining continuous-time Markov chains with graph transformation_, (in preparation).

More publications for DISCOTRESS are forthcoming.

For references of the individual simulation algorithms please see the individual .cpp files in the source code.

## Technical FAQs

### 

### Which sampling method should I use for my simulation?

For most Markov chains representing a realistic dynamical process of interest, the standard simulation algorithm (BKL) will most likely be too inefficient to simulate the dynamics. To choose the appropriate enhanced sampling algorithm, there are many factors to consider, including the information that is desired from the simulation. The characteristics of the topology and dynamics of the Markov chain will also influence which enhanced sampling method is optimal. For a thorough discussion concerning the choice of sampling method, see Ref. [1].

In general, kPS (see Refs. [1]-[3]) provides a highly efficient method to simulate the dynamics for metastable CTMCs or DTMCs.

### How should I obtain a predefined partitioning of the Markov chain for use with the enhanced sampling simulation algorithms?

The enhanced sampling methods described in Ref. [1] are based on a partitioning of the network into communities, which must accurately characterise the metastable states of the system. Be aware that the choice of this community structure strongly affects the efficiency of the simulation. The question of how to obtain this partitioning is therefore of critical importance.

Even for small networks, spectral methods (such as the original and robust Perron cluster cluster analysis methods, PCCA and PCCA+, respectively) will fail owing to numerical instability if the CTMC is metastable. Hence, for such networks, BACE (the Bayesian Agglomerative Clustering Engine) is a favourable alternative, but is not scalable. Unfortunately, many state-of-the-art community detection algorithms are based on optimisation of the modularity objective, and therefore _are liable to misrepresent the dynamics_. Multi-level regularised Markov clustering (MLR-MCL), an implementation of which is available at github.com/danieljsharpe/mlr_mcl, provides a suitable alternative.

### I want my simulation to go faster!

To achieve efficient dynamical simulations using enhanced sampling methods, the community structure must accurately characterise the metastable sets of nodes. It is therefore often worth exploring different algorithms and varying choices of parameters to obtain a suitable partitioning of the Markov chain. A common mistake when using the kPS and MCAMC algorithms is to not run a fixed number of standard BKL steps after basin escape iterations.

It is also possible to pre-process the Markov chain to reduce dimensionality and eliminate flickering, while maintaining an accurate representation of the original Markov chain, see Refs. [4] and [5].

## Contact the author

github.com/danieljsharpe
daniel.j.sharpe1995@gmail.com
