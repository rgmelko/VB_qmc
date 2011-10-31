#ifndef MEASURE_H
#define MEASURE_H

#include "head_proj.h"
#include "basis.h"
#include <fstream>

class Measure: public PARAMS
{
    private: 
      //observables here
      double Energy;

    public:

      Measure(){Energy = 0.0;};
      void zero();
      void measure_E(const Basis &);
      void output();
  
};

void Measure::zero(){

    Energy = 0.0;

}//zero

void Measure::measure_E(const Basis & basis){

    int n_0=0;

    for(int i=0; i<basis.OperatorList.size(); i++){

        if (basis.OperatorList[i].A == -1)
            n_0++;
    }//i

    Energy += 1.0*n_0;
    //cout<<Energy<<endl;

}//measure_E


void Measure::output(){

    double one_over_n;
	ofstream cfout;
	cfout.open("00.data",ios::app);

    one_over_n = Energy/(1.0*MCS_);
    cfout<<numSpin*h_x*2.0*m_ / one_over_n<<" ";
    cfout<<-(-1.0-h_x + h_x*2.0*m_ / one_over_n);
    cfout<<endl;

	cfout.close();

}//output


#endif
