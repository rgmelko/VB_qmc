#ifndef SPINSTATE_H
#define SPINSTATE_H

#include "head_proj.h"
#include "simparam.h"


class SpinState
{

	public:
		int numSpin;

		vector<int>  Sstate;   //VB basis

        SpinState(const PARAMS &);
        void print();

		//void filewrite(const int & num);
		//void fileread(const int & num);

};

//constructor
SpinState::SpinState(const PARAMS &p){

	numSpin = p.numSpin;
    for (int i=0; i<numSpin; i++) //initialize the Sz basis to null
        Sstate.push_back(-1);
}//constructor




//print
void SpinState::print(){

    cout<<"Spin basis: \n";

    for (int i=0; i<numSpin; i++) //initialize the Sz basis to null
      cout<<Sstate.at(i)<<endl;

}//print

#endif
