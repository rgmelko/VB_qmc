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
    Basis Vbeta_old(param);   //old bra for rejection
    Basis Valpha_old(param);  //old ket for rejection

	//initialize the spin state
    SpinState Z1(param); //the Sz basis state
    temp = Z1.SampleRandomState(mrand,Valpha,Vbeta);

	Measure Observ; //create measurement object

	//cout<<"("<<Valpha.TopoX()<<","<<Valpha.TopoY()<<")"<<endl;
	//cout<<"("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<endl;

	//choose your topological sector
	for (int i=0; i<10000; i++){
		temp = Valpha.LoopUpdate(mrand,param,Z1.Sstate);
		temp = Vbeta.LoopUpdate(mrand,param,Z1.Sstate);
		temp = Z1.SampleRandomState(mrand,Valpha,Vbeta); //sample spin state
		if( (Vbeta.TopoX() == param.Wx_) && (Vbeta.TopoY() == param.Wy_)) break;
	}
	cout<<"("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<endl;
	Valpha = Vbeta;
	cout<<"("<<Valpha.TopoX()<<","<<Valpha.TopoY()<<")"<<endl;
	temp = Z1.SampleRandomState(mrand,Valpha,Vbeta); //sample spin state

	Vbeta_old = Vbeta;
	Valpha_old = Valpha;

	int Walpha, Wbeta;  //detects winding number changes in the loop
	int i;
	for (int j=0; j<param.nBin_; j++){ //Monte Carlo binning loop

		Observ.zero(param);
		i = 0;
		while (i<param.MCS_){

			//for (int j=0; j<param.numSpin/2; j++){  //sample VB bonds
			//	Valpha.TwoBondUpdate(mrand,param,Z1.Sstate);
			//	Vbeta.TwoBondUpdate(mrand,param,Z1.Sstate);
			//}

			Wbeta = Vbeta.LoopUpdate(mrand,param,Z1.Sstate);
			Walpha = Valpha.LoopUpdate(mrand,param,Z1.Sstate);
			
			if (Walpha == 0 && Wbeta == 0){ //no winding number change

				temp = Z1.SampleRandomState(mrand,Valpha,Vbeta); //sample spin state
				//measurements
				Observ.measure_Cx(Vbeta, Valpha);
				Observ.record();
				i++;
				Valpha_old = Valpha;
				Vbeta_old= Vbeta;

			}//no Wnum change
			else{
				if (Walpha != 0) Valpha = Valpha_old;
				if (Wbeta != 0)  Vbeta = Vbeta_old;
			}//rejection based on winding number change

		}//i
		cout<<i<<endl;

		Observ.output(param);
	}//j

	//Vbeta.print();
    //cout<<"("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<endl;
   	cout<<"("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<endl;
	cout<<"("<<Valpha.TopoX()<<","<<Valpha.TopoY()<<")"<<endl;

    Valpha.filewrite(0);
	Vbeta.filewrite(1);

	//Z1.print();
    //Vbeta.print();
    //Valpha.print();
	//temp = Valpha|Vbeta;
	//cout<<temp<<endl;

	return 0;
};
