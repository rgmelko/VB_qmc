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

    int temp;
    PARAMS param; //read parameter file
    MTRand mrand(param.SEED_); //random number for metropolis

    //initialize your VB basis states
    Basis Vbeta(param);   //bra
    Basis Valpha(param);  //ket

	//initialize the spin state
    SpinState Z1(param); //the Sz basis state
    temp = Z1.SampleRandomState(mrand,Valpha,Vbeta);

	Measure Observ; //create measurement object

	cout<<"("<<Valpha.TopoX()<<","<<Valpha.TopoY()<<")"<<endl;
	cout<<"("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<endl;

	//choose your topological sector
	for (int i=0; i<100; i++){
		Valpha.LoopUpdate(mrand,param,Z1.Sstate);
		Vbeta.LoopUpdate(mrand,param,Z1.Sstate);
		temp = Z1.SampleRandomState(mrand,Valpha,Vbeta); //sample spin state
		//if( (Vbeta.TopoX() == param.Wx_) && (Vbeta.TopoY() == param.Wy_)) break;
		cout<<"("<<Valpha.TopoX()<<","<<Valpha.TopoY()<<")"<<endl;
		cout<<"("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<endl;
	}
	//Valpha.print();
	//Vbeta.print();
	Valpha.LoopUpdate(mrand,param,Z1.Sstate);
	Vbeta.LoopUpdate(mrand,param,Z1.Sstate);
	temp = Z1.SampleRandomState(mrand,Valpha,Vbeta); //sample spin state
	//if( (Vbeta.TopoX() == param.Wx_) && (Vbeta.TopoY() == param.Wy_)) break;
	cout<<"("<<Valpha.TopoX()<<","<<Valpha.TopoY()<<")"<<endl;
	cout<<"("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<endl;
	//Valpha.print();
	//Vbeta.print();

	//cout<<"("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<endl;
	//Valpha = Vbeta;
	//cout<<"("<<Valpha.TopoX()<<","<<Valpha.TopoY()<<")"<<endl;
	//temp = Z1.SampleRandomState(mrand,Valpha,Vbeta); //sample spin state
	return 1;

    //Monte Carlo binning loop
	for (int j=0; j<param.nBin_; j++){

		Observ.zero(param);
		for (int i=0; i<param.MCS_; i++){ //MC sampling steps

			for (int j=0; j<param.numSpin/2; j++){  //sample VB bonds
				Valpha.TwoBondUpdate(mrand,param,Z1.Sstate);
				Vbeta.TwoBondUpdate(mrand,param,Z1.Sstate);
				//Vbeta.LoopUpdate(mrand,param,Z1.Sstate);
				//Valpha.LoopUpdate(mrand,param,Z1.Sstate);
			}
			temp = Z1.SampleRandomState(mrand,Valpha,Vbeta); //sample spin state

			//measurements
			Observ.measure_Cx(Vbeta, Valpha);
			Observ.record();
		}//i

		Observ.output(param);
	}//j

	//Vbeta.print();
    //cout<<"("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<endl;
   
    Valpha.filewrite(0);
	Vbeta.filewrite(1);

	//Z1.print();
    //Vbeta.print();
    //Valpha.print();
	//temp = Valpha|Vbeta;
	//cout<<temp<<endl;


	return 0;
};
