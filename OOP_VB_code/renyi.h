#ifndef RENYI_H
#define RENYI_H

#include "head_proj.h"
#include "basis.h"
#include "projector.h"

class Renyi
{
    private: 
	  double entropy;
	  double TOTAL_H1;

    public: 
      void zero();
      void measure_H1(const Basis &, const Basis &);
      void record();
      void output(PARAMS &);


};//Renyi


void Renyi::zero(){

	TOTAL_H1 = 0; 

}//zero


void Renyi::measure_H1(const Basis & A, const Basis & B){

	int Nloop_num, Nloop_den; //number of loops in numerator and denominator
	Basis Vl(A);   //copy constructors
	Basis Vr(B);
    Nloop_den = Vl|Vr ; 

	int a,b, bond1, bond2;
	int old1, old2, old3, old4;

    int num_Swap = 4;
	for (int i=0; i<num_Swap; i++){
	  a = i;
	  b = i+A.LinX;

	  bond1 = Vr.VBasis[a]; 
	  bond2 = Vr.VBasis[b];

	  if (a != bond2) { //off-diagonal case: SWAP
		Vr.VBasis[a] = b;
		Vr.VBasis[b] = a;
		Vr.VBasis[bond1] = bond2;
		Vr.VBasis[bond2] = bond1;
	  }

	}//num_Swap 

    Nloop_num = Vl|Vr; //new overlap

	//a=0;
	//b=A.LinX;

	//bond1 = Vr.VBasis[a]; 
	//bond2 = Vr.VBasis[b];
	////diagonal operation
	//old1 = Vr.VBasis[a];
	//old2 = Vr.VBasis[b];
	//old3 = Vr.VBasis[bond1];
	//old4 = Vr.VBasis[bond2];
	//if (a == bond2) {
	//	if (b != bond1) cout<<"Measurement connection error \n";
	//	Nloop_num = Vl|Vr;  //Calculate overlap here
	//}
	////off diagonal operation
	//else{
	//	Vr.VBasis[a] = b;
	//	Vr.VBasis[b] = a;
	//	Vr.VBasis[bond1] = bond2;
	//	Vr.VBasis[bond2] = bond1;
	//	Nloop_num = Vl|Vr;  //Calculate overlap here
	//	Vr.VBasis[a] = old1;
	//	Vr.VBasis[b] = old2;
	//	Vr.VBasis[bond1] = old3;
	//	Vr.VBasis[bond2] = old4;
	//}

    entropy = 1.0*pow(2,Nloop_num - Nloop_den);

}//measure_H1

void Renyi::record(){

	TOTAL_H1 += entropy;

}//record

void Renyi::output(PARAMS & p){

	cout<<-log(TOTAL_H1/(2.0*p.MCS_))<<endl;

}//output

#endif
