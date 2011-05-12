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

int RestricTOPO=1;  //restricts the topological sector if 1, lets fluctuate if 0

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
    Valpha.SampleSpinState(mrand,Vbeta);

	Measure Observ; //create measurement object
	Renyi renyi(param.numSpin);

	//-----choose your topological sector
    while(1){
		temp = Valpha.LoopUpdate(mrand,param);
		temp = Vbeta.LoopUpdate(mrand,param);
		Valpha.SampleSpinState(mrand,Vbeta);
		if( (Vbeta.TopoX() == param.Wx_) && (Vbeta.TopoY() == param.Wy_)) break;
	}
	cout<<"beta: "; Vbeta.printTOPO();
    Vbeta.CopyTop();
    cout<<"beta: "; Vbeta.printTOPO();
	Valpha = Vbeta;
    cout<<"alpha: "; Valpha.printTOPO();

	Valpha.SampleSpinState(mrand,Vbeta);
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
		//	Valpha.TwoBondUpdate(mrand,param);
		//	Vbeta.TwoBondUpdate(mrand,param);
		//}
        //Valpha.SampleSpinState(mrand,Vbeta);
		//i++;

		Wbeta = Vbeta.LoopUpdate(mrand,param); //loop update
		Walpha = Valpha.LoopUpdate(mrand,param);

		//check winding number change
		if (RestricTOPO == 0){
			Valpha.SampleSpinState(mrand,Vbeta);
			i++;
		}
		else if (Walpha == 0 && Wbeta == 0){ 
			Valpha.SampleSpinState(mrand,Vbeta);
			Valpha_old = Valpha;
			Vbeta_old= Vbeta;
			i++;
		}//no Wnum change
		else{
			Valpha = Valpha_old;
			Vbeta = Vbeta_old;
		}//rejection based on winding number change

	} //*********End equilibriation
	Vbeta_old = Vbeta;
	Valpha_old = Valpha;

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
			//	Valpha.TwoBondUpdate(mrand,param);
			//	Vbeta.TwoBondUpdate(mrand,param);
			//}
        
			k=0;
			while (k<numLoops){
				Wbeta = Vbeta.LoopUpdate(mrand,param);
				Walpha = Valpha.LoopUpdate(mrand,param);

				if (RestricTOPO == 0){
					k++;
				}
				else if (Walpha == 0 && Wbeta == 0){ //no winding number change
					Valpha_old = Valpha;
					Vbeta_old= Vbeta;
					k++;
				}
				else{  //rejection based on winding number change
					Valpha = Valpha_old;
					Vbeta = Vbeta_old;
				}
			}//k

			//renyi.measure_H2(Valpha,Vbeta);      // for SWAP
			renyi.measure_ratio(Valpha,Vbeta);      // for ratio
            renyi.record();

			Observ.measure_Cx(Vbeta, Valpha);
			Observ.record();

			Valpha.SampleSpinState(mrand,Vbeta); //sample spin state
			Vbeta_old = Vbeta;
			Valpha_old = Valpha;

		}//i
		//cout<<i<<endl;

		renyi.output();
		Observ.output(param);
		Valpha.filewrite(0);
		Vbeta.filewrite(1);

	}//j
	//------------------------------------
    cout<<"beta : "; Vbeta.printTOPO();
    cout<<"alpha: "; Valpha.printTOPO();

	return 0;
};
