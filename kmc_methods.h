/*
Classes and functions for handling enhanced kinetic Monte Carlo simulations and propagating the trajectories
*/

#ifndef __KMC_METHODS_H_INCLUDED__
#define __KMC_METHODS_H_INCLUDED__

/* class containing functions to handle enhanced sampling KMC methods */
class KMC_Enhanced_Methods {

    public:

    KMC_Enhanced_Methods();
    ~KMC_Enhanced_Methods();
};

/* Weighted ensemble KMC */
class WE_KMC : KMC_Enhanced_Methods {

    public:

    WE_KMC();
    ~WE_KMC();

};

/* Forward flux sampling KMC */
class FFS_KMC : KMC_Enhanced_Methods {

    public:

    FFS_KMC();
    ~FFS_KMC();
};

/* class containing functions to propagate KMC trajectories */
class KMC_Standard_Methods {

    public:

    KMC_Standard_Methods();
    ~KMC_Standard_Methods();
    static void Gillespie(); // the stochastic simulation algorithm (SSA) of Gillespie - slow but exact
    static void Tau_Leaping(); // the approximate tau-leaping algorithm
    static void BKL(); // rejection-free algorithm of Bortz, Kalos and Lebowitz
};

#endif
