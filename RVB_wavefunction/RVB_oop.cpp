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

int RestricTOPO  = 1;  //restricts the topological sector if 1, lets fluctuate if 0
int Bond_Updates = 0; // 1 for bonds, 0 for loops

int main(){


    int temp;
	int numLoops; //number of loops to perform

    PARAMS param; //read parameter file
    MTRand mrand(param.SEED_); //random number for metropolis

	//log what's going on
	if (RestricTOPO == 0) cout<<"Winding numbers are fluctuating \n"; 
	else cout<<"Topological sector restricted to "<<param.Wx_<<" "<<param.Wy_<<endl;
	if (Bond_Updates == 1) cout<<"Bond updates \n";
	else cout<<"Loop Updates \n";

    //initialize your VB basis states
    Basis Vbeta(param);   //bra
    Basis Valpha(param);  //ket

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

	int Walpha, Wbeta;  //detects winding number changes in the loop
	int i,k;

	i = 0;
	//********Equilibriation
	while (i<param.EQL_){ 

		if (Bond_Updates == 1){ //perform bond updates
			for (int j=0; j<param.numSpin/2; j++){  //sample VB bonds
				Valpha.TwoBondUpdate(mrand,param);
				Vbeta.TwoBondUpdate(mrand,param);
			}
			Valpha.SampleSpinState(mrand,Vbeta);
			i++;
		}
		else{ //perform loop updates
			temp = Vbeta.LoopUpdate(mrand,param); //loop update
			temp = Valpha.LoopUpdate(mrand,param);
			Valpha.SampleSpinState(mrand,Vbeta);
			i++;
		}
		
	} //*********End equilibriation

	//temp = (Valpha.Scount + Vbeta.Scount)/(2*param.EQL_);
	//numLoops = param.numSpin/temp;
	//if (numLoops < 1) numLoops = 1; //adjust the number of loops
	numLoops = 10;

    //--------------Main Monte Carlo data binning loop
	for (int j=0; j<param.nBin_; j++){ 

		Observ.zero(param);
        renyi.zero();

		for (i=0; i<param.MCS_; i++){

			//Bond update
			if (Bond_Updates == 1){
				for (int j=0; j<param.numSpin/2; j++){  
					Valpha.TwoBondUpdate(mrand,param);
					Vbeta.TwoBondUpdate(mrand,param);
				}
				renyi.measure_ratio(Valpha,Vbeta);   //measure
				renyi.record();
			}
			else {//Loop update
				k=0;
				while (k<numLoops){
					temp = Vbeta.LoopUpdate(mrand,param);
					temp = Valpha.LoopUpdate(mrand,param);
					k++;
				}//k

				Walpha = Valpha.RightTopoNum(param);
				Wbeta = Vbeta.RightTopoNum(param);
				if ( (RestricTOPO == 0) || (Walpha == 0 && Wbeta == 0) ){ 
					//renyi.measure_H2(Valpha,Vbeta);      // for SWAP
					renyi.measure_ratio(Valpha,Vbeta);      // for ratio
					renyi.record();
				}
			}//end updates

			Observ.measure_Cx(Vbeta, Valpha);  //this could be moved up 
			Observ.record();

			Valpha.SampleSpinState(mrand,Vbeta); //sample spin state

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
