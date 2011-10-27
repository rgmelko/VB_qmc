// 
// TFIM_proj.cpp: the main body of the code
// A program to calculate the Renyi entropy in the transverse field Ising model
// Roger Melko, October 2011, University of Vermont
//
#include "head_proj.h"
#include "simparam.h"
#include "measure.h"
#include "basis.h"


int main(){


    PARAMS param; //read parameter file
    //param.printBst();

    MTRand mrand(param.SEED_); //random number for metropolis

    Basis Proj(mrand);
    //Proj.printBasis();

    for (int i=0; i<param.EQL_; i++){
        Proj.DiagonalUpdate(mrand);
        Proj.LinkedList();
        Proj.ClusterUpdate(mrand);
    }

    Measure observ;
    //observ.zero();
    for (int i=0; i<param.MCS_; i++){
        Proj.DiagonalUpdate(mrand);
        Proj.LinkedList();
        Proj.ClusterUpdate(mrand);
        observ.measure_E(Proj);
    }
    observ.output();


    return 0;

}

