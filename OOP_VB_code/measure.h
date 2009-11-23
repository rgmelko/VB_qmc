#ifndef MEASURE_H
#define MEASURE_H

#include "head_proj.h"
#include "basis.h"
#include "projector.h"

class Measure
{
    private: 
      //observables here
      double E_old; //energy
      double E2_old;  //energy
      double C_old; //C(L/2,L/2)

      void MapLoops(const Basis &, const Basis &, vector<int> & temp); //creates overlap map of numbered loops


    public:
      double TOT_energy;   //energy
      double TOT_energy2;   //energy
      double TOT_cL_2;    //C(L/2,L/2)

      void zero();
      void measure_energy(const Basis &, const Basis &, const PARAMS &);
      void measure_energy2(const Basis &, const Basis &, const PARAMS &);
      void measure_CL2L2(const Basis &, const Basis &);
      void record();
      void output(const PARAMS &);
  
};


void Measure::zero(){

    E_old=0; 
    E2_old=0;
    C_old=0;

    TOT_energy = 0.0;
    TOT_energy2 = 0.0;
    TOT_cL_2 = 0.0;
}

void Measure::measure_energy(const Basis & A, const Basis & B, const PARAMS & p){
//******************** Energy *********************	

	int Nloop_num, Nloop_den; //number of loops in numerator and denominator

	//first, calculate overlap <A|B>
	Basis Vl(A);
	Basis Vr(B);   //copy constructors

	Nloop_den = Vl|Vr ; 
	//cout<<N_loop<<"\n";
	//cout.flush();

	int a,b, bond1, bond2;
	int old1, old2, old3, old4;
	double Enrgy=0;
	int off_d_count;
	for (int i=0; i<Vr.numLattB; i++){ //Propogate Vr

		//a=Vr.Bst.at(i).A;
		//b=Vr.Bst.at(i).B;
		a=p.Bst.at(i).A;
		b=p.Bst.at(i).B;

		bond1 = Vr.VBasis[a]; 
		bond2 = Vr.VBasis[b];
		//diagonal operation
		old1 = Vr.VBasis[a];
		old2 = Vr.VBasis[b];
		old3 = Vr.VBasis[bond1];
		old4 = Vr.VBasis[bond2];
		if (a == bond2) {
			if (b != bond1) cout<<"Measurement connection error \n";
			Nloop_num = Vl|Vr;  //Calculate overlap here
			off_d_count = 0;
		}
		//off diagonal operation
		else{
			Vr.VBasis[a] = b;
			Vr.VBasis[b] = a;
			Vr.VBasis[bond1] = bond2;
			Vr.VBasis[bond2] = bond1;
			Nloop_num = Vl|Vr;  //Calculate overlap here
			Vr.VBasis[a] = old1;
			Vr.VBasis[b] = old2;
			Vr.VBasis[bond1] = old3;
			Vr.VBasis[bond2] = old4;
			off_d_count = 1;
		}

		Enrgy += 1.0*pow(2,Nloop_num - Nloop_den - off_d_count);

	}//i

    //cout<<Enrgy<<endl;
	E_old = Enrgy; 
	//E_new = Enrgy; 

}//measure_energy 

void Measure::measure_energy2(const Basis & A, const Basis & B, const PARAMS & p){
//******************** Ann's Energy *********************	

	vector<int> is_in_loop;  //records whether a spin is counted in a loop 
	MapLoops(A,B,is_in_loop);

    int a,b;
    int m_diff = 0;
	for (int i=0; i<A.numLattB; i++){

		//a=A.Bst.at(i).A;
		//b=A.Bst.at(i).B;
		a=p.Bst.at(i).A;
		b=p.Bst.at(i).B;

		if (is_in_loop.at(a) == 0 || is_in_loop.at(b) == 0) cout<<"Energy 2 error \n";
		else if (is_in_loop.at(a) != is_in_loop.at(b)) m_diff ++;
	}

	//E2_new = 0.75*(m_diff - A.numLattB);
	E2_old = 0.75*(m_diff - A.numLattB);

}//energy2


void Measure::measure_CL2L2(const Basis & A, const Basis & B){
//********************C(L/2,L/2)*********************	

	vector<int> is_in_loop;  //records whether a spin is counted in a loop 
	MapLoops(A,B,is_in_loop);

	if (is_in_loop.at(0) == is_in_loop.at(A.LinX*A.LinX/2+A.LinX/2) ) //fixed the (L/2,L/2) bug
		C_old= 0.75;
	else 
		C_old= 0;


}//measure_CL2L2


void Measure::MapLoops(const Basis & A, const Basis & B, vector<int> & is_in_loop){
//****************creates overlap map of numbered loops ****************	

	is_in_loop.assign(B.VBasis.size(),0);

	int next;
	int Nloop = 0;

	int count = 1;
	for (int i=0; i<B.VBasis.size(); i++){

		if (is_in_loop.at(i) == 0){
			is_in_loop.at(i) = count;
			next = A.VBasis.at(i); //V_A basis
			while (is_in_loop.at(next) == 0){

				if  (is_in_loop.at(next) != 0) cout<<"loop error 1 \n";
				else is_in_loop.at(next) = count;

				next = B.VBasis.at(next);      //V_B basis
				if  (is_in_loop.at(next) == 0) is_in_loop.at(next) = count; 
				else break;

				next = A.VBasis.at(next); //V_A basis 
			}//while

			count++;
			Nloop ++;

		}//if
	}//i

}//MapLoops

void Measure::record(){

    TOT_energy += E_old;
    TOT_energy2 += E2_old;
    TOT_cL_2 += C_old;

}//update

void Measure::output(const PARAMS & p){

    //TOT_energy/= (2.0*p.MCS_);
    //cout<<-TOT_energy/p.numSpin+0.25*p.numLattB/p.numSpin<<" ";
    cout<<TOT_energy2/(2.0*p.MCS_ * p.numSpin)<<"\n";
    //cout<<TOT_cL_2/(2.0*p.MCS_)<<endl; //factor of 2 for 2 projector samples

}//output


#endif
