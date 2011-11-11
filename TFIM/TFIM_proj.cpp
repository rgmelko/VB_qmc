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
    param.printBst();

    MTRand mrand(param.SEED_); //random number for metropolis

    Basis Proj(mrand);
    //Proj.printBasis();

    int Loopsize2; //loopsize squared for M^2 measurement

    for (int i=0; i<param.EQL_; i++){
        Proj.DiagonalUpdate(mrand);
        Proj.LinkedList();
        //Proj.printLinkedList();
        Proj.ClusterUpdate(mrand,Loopsize2);
    }


    vector<int> inA(param.numSpin/2,0); //size of physical spins
    inA[0] = 1;

    Measure observ;

    for (int j=0; j< param.nBin_; j++){
        observ.zero();
        for (int i=0; i<param.MCS_; i++){
            Proj.DiagonalUpdate(mrand);
            Proj.LinkedList();
            Proj.printBasis();
            Proj.printLinkedList();
            Proj.ClusterUpdate(mrand,Loopsize2);
            Proj.printBasis();
            cout<<"swapped "<<Proj.SWAP(inA)<<endl;
            return 0;
            observ.measure_E(Proj);
            observ.measure_M(Proj, Loopsize2);
        }
        observ.output();
    }


    return 0;

}

