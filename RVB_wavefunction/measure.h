#ifndef MEASURE_H
#define MEASURE_H

#include "head_proj.h"
#include "basis.h"
#include <fstream>

class Measure
{
    private: 
      //observables here
      double C_old; //C(L/2,L/2)
      vector<double> C_x;
      vector<double> C_x_anc; //the ancilary

      void MapLoops(const Basis &, const Basis &, vector<int> & temp); //creates overlap map of numbered loops


    public:
      double TOT_cL_2;    //C(L/2,L/2)

      void zero(const PARAMS &);
      void measure_CL2L2(const Basis &, const Basis &);
      void measure_Cx(const Basis &, const Basis &);
      void record();
      void output(const PARAMS &);
  
};


void Measure::zero(const PARAMS & p){

    C_old=0;
    TOT_cL_2 = 0.0;

    C_x.clear();
    C_x_anc.clear();
    for (int i = 0; i<p.nX_; i++){
        C_x.push_back(0.0);
        C_x_anc.push_back(0.0);
    }
}//zero

void Measure::measure_Cx(const Basis & A, const Basis & B){
//********************C(0,x)*********************	

    vector<int> is_in_loop;  //records whether a spin is counted in a loop 
	MapLoops(A,B,is_in_loop);

    int half = A.LinX * A.LinX;
    int j;
    for (int i = 0; i<A.LinX; i++){
        if (is_in_loop.at(0) == is_in_loop.at(i) ) {
            if (i%2 == 0) //same sublattice
                C_x.at(i) += 0.75;
            else
                C_x.at(i) -= 0.75; //different sublattice
        }//real lattice

        j = i+half;
        if (is_in_loop.at(half) == is_in_loop.at(j) ) {
            if (i%2 == 0) //same sublattice
                C_x_anc.at(i) += 0.75;
            else
                C_x_anc.at(i) -= 0.75; //different sublattice
        }//ancillary

    }

}//measure_Cx


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

    TOT_cL_2 += C_old;

}//update

void Measure::output(const PARAMS & p){

	ofstream cfout;
	cfout.open("00.data",ios::app);

    //cfout<<p.MCS_<<endl;
    //cfout<<TOT_cL_2/(1.0*p.MCS_)<<endl; 

    for (int i=0; i<p.nX_/2; i++)
        //cfout<<i<<" "<<C_x.at(i)/(1.0*p.MCS_)<<endl;
        cfout<<C_x.at(i)/(1.0*p.MCS_)<<" ";

    cfout<<endl;

	cfout.close();

	cfout.open("01.data",ios::app);

    //cfout<<p.MCS_<<endl;
    //cfout<<TOT_cL_2/(1.0*p.MCS_)<<endl; 

    for (int i=0; i<p.nX_/2; i++)
        //cfout<<i<<" "<<C_x.at(i)/(1.0*p.MCS_)<<endl;
        cfout<<C_x_anc.at(i)/(1.0*p.MCS_)<<" ";

    cfout<<endl;

	cfout.close();

}//output


#endif
