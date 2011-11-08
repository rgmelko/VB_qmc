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
	  double Mag1;
	  double Mag2;

    public:

      Measure(){Energy = 0.0; Mag1 = 0.0; Mag2 = 0.0;};
      void zero();
      void measure_E(const Basis &);
      void measure_M(const Basis &, const int &);
      void output();
  
};

void Measure::zero(){

    Energy = 0.0;
    Mag1 = 0.0;
    Mag2 = 0.0;

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


void Measure::measure_M(const Basis & basis, const int & L2){

	int m_0 = 0;

    vector<int> S_prop; 
    S_prop = basis.S_left; 

    for(int i=0; i<basis.OperatorList.size()/2; i++){

        if (basis.OperatorList[i].A == -2) //this is a off-diagonal site operator
            S_prop[basis.OperatorList[i].B] = S_prop[basis.OperatorList[i].B]^1 ; //spin flip

	}

	for (int i=0; i<basis.numSpin; i++)
		m_0 += (2*S_prop[i]-1);


    Mag1 += 1.0*m_0*m_0; //m^2

    Mag2 += 1.0*L2;

}//measure_M


void Measure::output(){

    double one_over_n;
	ofstream cfout;
	cfout.open("00.data",ios::app);

    one_over_n = Energy/(1.0*MCS_);
    cfout<<numSpin*h_x*2.0*m_ / one_over_n<<" ";
    cfout<<-(numSpin*h_x*2.0*m_ / one_over_n - numSpin*h_x - numLattB)/numSpin<<" ";
    //cfout<<-(-1.0-h_x + h_x*2.0*m_ / one_over_n)<<" "; //wrong for OBC?
	cfout<<Mag1/(1.0*MCS_*1.0*numSpin*numSpin)<<" ";
	cfout<<Mag2/(1.0*MCS_*1.0*numSpin*numSpin);
    cfout<<endl;

	cfout.close();

}//output


#endif
