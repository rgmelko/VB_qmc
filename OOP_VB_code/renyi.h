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
      void measure_H1(const Basis &, const Basis &, const int &);
      void record();
      void output(PARAMS &);


};//Renyi


void Renyi::zero(){

	TOTAL_H1 = 0; 

}//zero


void Renyi::measure_H1(const Basis & A, const Basis & B, const int & num_Swap){

	int Nloop_num, Nloop_den; //number of loops in numerator and denominator
	Basis Vl(A);   //copy constructors
	Basis Vr(B);
    Nloop_den = Vl|Vr ; 

	int a,b, bond1, bond2;
	int old1, old2, old3, old4;

	for (int i=0; i<num_Swap; i++){
	  a = i;
	  b = i+A.LinX; //b in the other layer

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

    entropy = 1.0*pow(2,Nloop_num - Nloop_den);

}//measure_H1

void Renyi::record(){

	TOTAL_H1 += entropy;

}//record

void Renyi::output(PARAMS & p){

	cout<<-log(TOTAL_H1/(2.0*p.MCS_))<<endl;

}//output

#endif
