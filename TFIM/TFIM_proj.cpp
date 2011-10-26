// 
// TFIM_proj.cpp: the main body of the code
// A program to calculate the Renyi entropy in the transverse field Ising model
// Roger Melko, October 2011, University of Vermont
//
#include "head_proj.h"
#include "simparam.h"
#include "basis.h"


int main(){


    PARAMS param; //read parameter file
    //param.printBst();

    MTRand mrand(param.SEED_); //random number for metropolis

    Basis Proj(mrand);
    Proj.printBasis();

    Proj.DiagonalUpdate(mrand);
    Proj.printBasis();

    Proj.LinkedList();
    Proj.printLinkedList();
    Proj.printBasis();


    return 0;

}

