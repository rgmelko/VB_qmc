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

    //param.printBst();
    //param.printPst();

    Basis Vbeta(param);   //bra
    Basis Valpha(param);  //ket

    //Vbeta.print();
    //Valpha.print();

    SpinState Z1(param); //the Sz basis state

    int temp;

	//initialize the spin state
    temp = Z1.SampleRandomState(mrand,Valpha,Vbeta);
 
	for (int i=0; i<10000; i++){
		for (int j=0; j<param.numSpin/2; j++){  //sample VB bonds
			Valpha.TwoBondUpdate(mrand,param,Z1.Sstate);
			Vbeta.TwoBondUpdate(mrand,param,Z1.Sstate);
		}
		temp = Z1.SampleRandomState(mrand,Valpha,Vbeta); //sample spin state
	}

	Z1.print();
    Vbeta.print();
    Valpha.print();

	temp = Valpha|Vbeta;
	cout<<temp<<endl;


	return 0;
};
