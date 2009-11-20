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
	  vector<int> inAreg; //inside the "A region
	  int nSwap;

    public: 
	  Renyi(const int &, const int &);
      void zero();
      void measure_H1(const Basis &, const Basis &);
      void record();
      void output(PARAMS &);


};//Renyi

Renyi::Renyi(const int & Lsize, const int & num_Swap){

	nSwap = num_Swap;

    inAreg.assign(2*Lsize,0);
	for (int i=0; i<num_Swap; i++){
		inAreg.at(i)=1;
		inAreg.at(i+Lsize)=1;
	}


};//constructor


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

    //cout<<nSwap<<" ";
    //Vr.print();
	for (int i=0; i<nSwap; i++){
      
	  a = i;
	  b = i+A.LinX; //b in the other layer

	  bond1 = Vr.VBasis[a]; 
	  bond2 = Vr.VBasis[b];

	  if (inAreg.at(bond2) == 1 )
		  Vr.VBasis[a] = bond2 - A.LinX;
      else{
		  Vr.VBasis[a] = bond2;
		  Vr.VBasis[bond2] = a;
	  }

	  if (inAreg.at(bond1) == 1 )
		  Vr.VBasis[b] = bond1 + A.LinX;
      else{
		  Vr.VBasis[b] = bond1;
		  Vr.VBasis[bond1] = b;
	  }

	}//nSwap 
	//Vr.print();
	//cout<<endl;

    Nloop_num = Vl|Vr; //new overlap

    entropy = 1.0*pow(2,Nloop_num - Nloop_den);

}//measure_H1

void Renyi::record(){

	TOTAL_H1 += entropy;

}//record

void Renyi::output(PARAMS & p){

	cout<<-log(TOTAL_H1/(2.0*p.MCS_))<<" ";

}//output

#endif
