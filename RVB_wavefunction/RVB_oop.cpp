// 
// RVB_oop.cpp: the main body of the program
// A OOP VB basis projector code for the RVB wavefunction
// Roger Melko, April 2011, Santa Barbara Station Q
//
#include "head_proj.h"
#include "simparam.h"
#include "basis.h"
#include "spinstate.h"
#include "measure.h"

int main(){

    PARAMS param;
    MTRand mrand(param.SEED_); //random number for metropolis

    Basis Vbeta(param);   //bra
    Basis Valpha(param);  //ket

    SpinState Z1(param); //the Sz basis state

    int temp;

    Measure Observ; //create measurement object
	Observ.zero(param);

	//initialize the spin state
    temp = Z1.SampleRandomState(mrand,Valpha,Vbeta);
 
    //MC sampling
	for (int i=0; i<param.MCS_; i++){

		for (int j=0; j<param.numSpin/2; j++){  //sample VB bonds
			Valpha.TwoBondUpdate(mrand,param,Z1.Sstate);
			Vbeta.TwoBondUpdate(mrand,param,Z1.Sstate);
		}
		temp = Z1.SampleRandomState(mrand,Valpha,Vbeta); //sample spin state

        //measurements
		Observ.measure_Cx(Vbeta, Valpha);
		Observ.record();
	}

	Observ.output(param);

	//Z1.print();
    //Vbeta.print();
    //Valpha.print();
	//temp = Valpha|Vbeta;
	//cout<<temp<<endl;


	return 0;
};
