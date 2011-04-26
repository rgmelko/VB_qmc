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
	int numLoops; //number of loops to perform

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

	//-----choose your topological sector
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
	//---------------------------------------

	Vbeta_old = Vbeta;
	Valpha_old = Valpha;

	int Walpha, Wbeta;  //detects winding number changes in the loop
	int i,k;

	i = 0;
	//********Equilibriation
	while (i<param.EQL_){ 

		Wbeta = Vbeta.LoopUpdate(mrand,param,Z1.Sstate); //oop update
		Walpha = Valpha.LoopUpdate(mrand,param,Z1.Sstate);

		if (Walpha == 0 && Wbeta == 0){ //no winding number change
			temp = Z1.SampleRandomState(mrand,Valpha,Vbeta); //sample spin state
			i++;
			Valpha_old = Valpha;
			Vbeta_old= Vbeta;
		}//no Wnum change
		else{
			if (Walpha != 0) Valpha = Valpha_old;
			if (Wbeta != 0)  Vbeta = Vbeta_old;
		}//rejection based on winding number change

	}
	//*********End equilibriation

	temp = (Valpha.Scount + Vbeta.Scount)/(2*param.EQL_);
	numLoops = param.numSpin/temp;
	if (numLoops < 1) numLoops = 1; //adjust the number of loops

    //--------------Main Monte Carlo data binning loop
	for (int j=0; j<param.nBin_; j++){ 

		Observ.zero(param);
		i = 0;
		for (i=0; i<param.MCS_; i++){

			//for (int j=0; j<param.numSpin/2; j++){  //sample VB bonds
			//	Valpha.TwoBondUpdate(mrand,param,Z1.Sstate);
			//	Vbeta.TwoBondUpdate(mrand,param,Z1.Sstate);
			//}

        
		    k=0;
			while (k<numLoops){
				Wbeta = Vbeta.LoopUpdate(mrand,param,Z1.Sstate);
				Walpha = Valpha.LoopUpdate(mrand,param,Z1.Sstate);

				if (Walpha == 0 && Wbeta == 0){ //no winding number change
					Valpha_old = Valpha;
					Vbeta_old= Vbeta;
					k++;
				}
				else{  //rejection based on winding number change
					if (Walpha != 0) Valpha = Valpha_old;
					if (Wbeta != 0)  Vbeta = Vbeta_old;
				}
			}//k

			Observ.measure_Cx(Vbeta, Valpha);
			Observ.record();

			temp = Z1.SampleRandomState(mrand,Valpha,Vbeta); //sample spin state
		}//i
		cout<<i<<endl;

		Observ.output(param);
		Valpha.filewrite(0);
		Vbeta.filewrite(1);

	}//j
	//------------------------------------

	//Vbeta.print();
    //cout<<"("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<endl;
   	cout<<"("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<endl;
	cout<<"("<<Valpha.TopoX()<<","<<Valpha.TopoY()<<")"<<endl;


	//Z1.print();
    //Vbeta.print();
    //Valpha.print();
	//temp = Valpha|Vbeta;
	//cout<<temp<<endl;

	return 0;
};
