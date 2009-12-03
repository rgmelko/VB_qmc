#ifndef RENYI_H
#define RENYI_H

#include "head_proj.h"
#include "basis.h"
#include "projector.h"

class Renyi
{
    private: 
	  vector<double> entropy;
	  vector<double> TOTAL_H2;
	  vector<int> inAreg; //inside the "A region
	  int nSwap;
	  int Nloop_den;   //number of loops in the denominator

    public: 
	  Renyi(const int &);
      void zero();
      void measure_H2(const Basis &, const Basis &);
      int calc_SWAP(const Basis &, const Basis &, const int &);
      void record();
      void output(PARAMS &);


};//Renyi

Renyi::Renyi(const int & Lsize){

	nSwap = Lsize;

	entropy.assign(nSwap-1,0);  //resize and initialize entropy
	TOTAL_H2.assign(nSwap-1,0);  //resize and initialize entropy total

     //********TEMPORARY FIX: SEE calc_SWAP
    //inAreg.assign(2*Lsize,0);
	//for (int i=0; i<nSwap; i++){
	//	inAreg.at(i)=1;
	//	inAreg.at(i+Lsize)=1;
	//}


};//constructor


void Renyi::zero(){

	TOTAL_H2.assign(nSwap-1,0);

}//zero


void Renyi::measure_H2(const Basis & A, const Basis & B){

	Basis Vl(A);   //copy constructors
	Basis Vr(B);
    Nloop_den = Vl|Vr ; 

    int Nloop_num;
	for (int r=1; r<nSwap; r++){
        Nloop_num = calc_SWAP(A,B,r);
		entropy.at(r-1) = 1.0*pow(2,Nloop_num - Nloop_den);
	}

}//measure_H2


int Renyi::calc_SWAP(const Basis & A, const Basis & B, const int & X){

    //  //1D
    //  inAreg.assign(2*A.LinX,0);
	//  for (int i=0; i<X; i++){
	//  	inAreg.at(i)=1;
	//  	inAreg.at(i+A.LinX)=1;
	//  }

	//2D
	inAreg.assign(A.numSpin,0);
	for (int i=0; i<X; i++){
		for (int j=0; j<X; j++){
			inAreg.at(i+j*A.LinX)=1;
			inAreg.at(i+j*A.LinX+A.numSpin/2)=1;
		}
	}

	//for (int i=0; i < inAreg.size(); i++)
	//	cout<<i<<" "<<inAreg.at(i)<<endl;

	int Nloop_num; //number of loops in numerator and denominator
	Basis Vl(A);   //copy constructors
	Basis Vr(B);
    //Nloop_den = Vl|Vr ; 

	int a,b, bond1, bond2;
	int old1, old2, old3, old4;

    int sA; //spin in region "A"
	for (int i=0; i<X; i++){ //X is the maximum *linear* distance to take region "A"
		for (int j=0; j<X; j++){ 

			sA = i+j*A.LinX;
			if (inAreg.at(sA) != 1) cout<<"Renyi error: A reg\n";

			a = sA;
			b = sA + A.numSpin/2; //b in the other layer

			bond1 = Vr.VBasis[a]; 
			bond2 = Vr.VBasis[b];

			if (inAreg.at(bond2) == 1 )
				Vr.VBasis[a] = bond2 - A.numSpin/2;
			else{
				Vr.VBasis[a] = bond2;
				Vr.VBasis[bond2] = a;
			}

			if (inAreg.at(bond1) == 1 )
				Vr.VBasis[b] = bond1 + A.numSpin/2;
			else{
				Vr.VBasis[b] = bond1;
				Vr.VBasis[bond1] = b;
			}

		}//j
	}//i

    Nloop_num = Vl|Vr; //new overlap

    //entropy = 1.0*pow(2,Nloop_num - Nloop_den);
    //return 1.0*pow(2,Nloop_num - Nloop_den);
    return Nloop_num;

}//calc_SWAP


void Renyi::record(){

    for (int i=0; i<TOTAL_H2.size(); i++)
	   TOTAL_H2.at(i) += entropy.at(i);

}//record


void Renyi::output(PARAMS & p){

	ofstream cfout;
	cfout.open("00.data",ios::app);

    for (int i=0; i<TOTAL_H2.size(); i++)
	  cfout<<-log(TOTAL_H2.at(i)/(2.0*p.MCS_))<<" ";
	cfout<<endl;

	cfout.close();

}//output

#endif
