// 
// RVB_oop.cpp: the main body of the program
// A OOP VB basis projector code for the RVB wavefunction
// Roger Melko, April 2011, Santa Barbara Station Q
//
#include "head_proj.h"
#include "simparam.h"
#include "basis.h"
#include "measure.h"
#include "renyi.h"

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
	Vbeta.SWAP(param.ratio_);
    Valpha.SampleSpinState(mrand,Vbeta);
	Vbeta.SWAP(param.ratio_);

	Measure Observ; //create measurement object
	Renyi renyi(param.nX_,param.ratio_);

	//cout<<"("<<Valpha.TopoX()<<","<<Valpha.TopoY()<<")"<<endl;
	//cout<<"("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<endl;

	//-----choose your topological sector
	for (int i=0; i<100000; i++){
		temp = Valpha.LoopUpdate(mrand,param);
		temp = Vbeta.LoopUpdate(mrand,param);
		Vbeta.SWAP(param.ratio_);
		Valpha.SampleSpinState(mrand,Vbeta);
		Vbeta.SWAP(param.ratio_);
		if( (Vbeta.TopoX() == param.Wx_) && (Vbeta.TopoY() == param.Wy_)) break;
	}
	cout<<"beta("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<", ";
	cout<<"("<<Vbeta.TopoXanc()<<","<<Vbeta.TopoYanc()<<")"<<endl;
	Vbeta.CopyTop();
	cout<<"beta("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<", ";
	cout<<"("<<Vbeta.TopoXanc()<<","<<Vbeta.TopoYanc()<<")"<<endl;
	Valpha = Vbeta;
	cout<<"alpha("<<Valpha.TopoX()<<","<<Valpha.TopoY()<<")"<<", ";
	cout<<"("<<Valpha.TopoXanc()<<","<<Valpha.TopoYanc()<<")"<<endl;
	Vbeta.SWAP(param.ratio_);
	Valpha.SampleSpinState(mrand,Vbeta);
	Vbeta.SWAP(param.ratio_);
	//---------------------------------------
	Valpha.filewrite(0);
	Vbeta.filewrite(1); //save configuration file

	Vbeta_old = Vbeta;
	Valpha_old = Valpha;

	int Walpha, Wbeta;  //detects winding number changes in the loop
	int i,k;

	i = 0;
	//********Equilibriation
	while (i<param.EQL_){ 

		//for (int j=0; j<param.numSpin/2; j++){  //sample VB bonds
		//	Valpha.TwoBondUpdate(mrand,param,Z1.Sstate);
		//	Vbeta.TwoBondUpdate(mrand,param,Z1.Sstate);
		//}
		//temp = Z1.SampleRandomState(mrand,Valpha,Vbeta); //sample spin state
		//i++;

		Wbeta = Vbeta.LoopUpdate(mrand,param); //oop update
		Walpha = Valpha.LoopUpdate(mrand,param);

		if (Walpha == 0 && Wbeta == 0){ //no winding number change
			Vbeta.SWAP(param.ratio_);
			Valpha.SampleSpinState(mrand,Vbeta);
			Vbeta.SWAP(param.ratio_);
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
	//cout<<"end eq \n";

	//temp = (Valpha.Scount + Vbeta.Scount)/(2*param.EQL_);
	//numLoops = param.numSpin/temp;
	//if (numLoops < 1) numLoops = 1; //adjust the number of loops
	numLoops = 10;

    //--------------Main Monte Carlo data binning loop
	for (int j=0; j<param.nBin_; j++){ 

		Observ.zero(param);
        renyi.zero();

		for (i=0; i<param.MCS_; i++){

			//for (int j=0; j<param.numSpin/2; j++){  //sample VB bonds
			//	Valpha.TwoBondUpdate(mrand,param,Z1.Sstate);
			//	Vbeta.TwoBondUpdate(mrand,param,Z1.Sstate);
			//}
        
		    k=0;
			while (k<numLoops){
				Wbeta = Vbeta.LoopUpdate(mrand,param);
				Walpha = Valpha.LoopUpdate(mrand,param);

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

			renyi.measure_ratio(Valpha,Vbeta,param.ratio_);      // for ratio
            renyi.record();

			Observ.measure_Cx(Vbeta, Valpha);
			Observ.record();

			Vbeta.SWAP(param.ratio_);
			Valpha.SampleSpinState(mrand,Vbeta); //sample spin state
			Vbeta.SWAP(param.ratio_);
		}//i
		//cout<<i<<endl;

		renyi.output(param);
		Observ.output(param);
		Valpha.filewrite(0);
		Vbeta.filewrite(1);

	}//j
	//------------------------------------

	cout<<"beta("<<Vbeta.TopoX()<<","<<Vbeta.TopoY()<<")"<<", ";
	cout<<"("<<Vbeta.TopoXanc()<<","<<Vbeta.TopoYanc()<<")"<<endl;
	cout<<"alpha("<<Valpha.TopoX()<<","<<Valpha.TopoY()<<")"<<", ";
	cout<<"("<<Valpha.TopoXanc()<<","<<Valpha.TopoYanc()<<")"<<endl;

	return 0;
};
