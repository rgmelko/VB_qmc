#ifndef SPINSTATE_H
#define SPINSTATE_H

#include "head_proj.h"
#include "simparam.h"
#include "basis.h"


class SpinState
{

	public:
		int numSpin;

		vector<int>  Sstate;   //VB basis

        SpinState(const PARAMS &);

        int SampleRandomState(const Basis&, const Basis&);
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


//A function which chooses a random spin state compatible with 
//the two input VB basis states.  See also Basis::operator|
int SpinState::SampleRandomState(const Basis& alpha, const Basis & beta){

    vector<int> is_in_loop;  //records whether a spin is counted in a loop 
    is_in_loop.assign(beta.VBasis.size(),0);

    int next;
    int Nloop = 0;

    int spinval;
    for (int i=0; i<beta.VBasis.size(); i++){
        spinval = i%2; //change to random
        Sstate.at(i) = spinval;
        if (is_in_loop.at(i) == 0){

            is_in_loop.at(i) = 1;
            next = alpha.VBasis.at(i); //V_A basis
            while (is_in_loop.at(next) == 0){
                if  (is_in_loop.at(next) == 1) cout<<"loop error 1 \n";
                else is_in_loop.at(next) = 1;

                spinval = spinval^1;  //bit flip
                Sstate.at(next) = spinval;

                next = beta.VBasis.at(next);      //V_B basis
                if  (is_in_loop.at(next) != 1) is_in_loop.at(next) = 1; 
                else break;

                spinval = spinval^1;  //bit flip
                Sstate.at(next) = spinval;

                next = alpha.VBasis.at(next); //V_A basis 
            }//while

            Nloop ++;

        }//if
    }//i

    return Nloop;

}


//print
void SpinState::print(){

    cout<<"Spin basis: \n";

    for (int i=0; i<numSpin; i++) //initialize the Sz basis to null
      cout<<Sstate.at(i)<<endl;

}//print

#endif
