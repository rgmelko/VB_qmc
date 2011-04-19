// 
// RVB_oop.cpp: the main body of the program
// A OOP VB basis projector code for the RVB wavefunction
// Roger Melko, April 2011, Santa Barbara Station Q
//
#include "head_proj.h"
#include "simparam.h"
#include "basis.h"
#include "spinstate.h"

int main(){

    PARAMS param;
    MTRand mrand(param.SEED_); //random number for metropolis

    Basis Vbeta(param);   //bra
    Basis Valpha(param);  //ket

    Vbeta.print();
    Valpha.print();

    SpinState Z1(param); //the Sz basis state
    Z1.print();

    //int temp;
    //temp = Vbeta|Valpha;
    //cout<<temp<<endl;




	return 0;
};
