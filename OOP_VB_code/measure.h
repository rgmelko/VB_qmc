#ifndef MEASURE_H
#define MEASURE_H

#include "head_proj.h"
#include "basis.h"
#include "projector.h"

class Measure
{
    private: 
      //observables here
      double E_old, E_new;   //energy
      double C_old, C_new;   //C(L/2,L/2)


    public:
      double TOT_energy;   //energy
      double TOT_cL_2;    //C(L/2,L/2)

      void zero();
      void measure(const Basis &, const Basis &);
      void update();
      void record();
      void output(const int &);
  
};


void Measure::zero(){

    TOT_energy = 0;
    TOT_cL_2 = 0;
}

void Measure::measure(const Basis & A, const Basis & B){

    vector<int> is_in_loop;  //records whether a spin is counted in a loop 
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

    if (is_in_loop.at(0) == is_in_loop.at(A.LinX/2) )
        C_new= 0.75;
    else 
        C_new= 0;


}//measure

void Measure::update(){

    C_old = C_new;

}//update

void Measure::record(){

    TOT_cL_2 += C_old;

}//update

void Measure::output(const int & MCS){

    cout<<TOT_cL_2/(2.0*MCS)<<endl; //factor of 2 for 2 projector samples

}//output


#endif
