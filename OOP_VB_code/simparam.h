#ifndef PARAMSIM_H
#define PARAMSIM_H

#include <fstream>
#include "head_proj.h"

//Class to read in the simulation parameters from a file
class PARAMS
{
    public:
        int LinX;  //linear lattice size
        int numSpin; 
        int numLattB;  //number of lattice bonds is 2N
		int EQL_; //the number of Monte Carlo steps
		int MCS_; //the number of Monte Carlo steps
        long SEED_;
        PARAMS();

        //lattice spin coordination numbers
        vector <index2> Bst;

        void printBst();
		
    protected:

        //derived constants
        int numVB; //number of valence bonds = N/2
        long int NN_;  //This is the length of the operator string
        int nX_;     //linear size of square lattice
        int sample_; //number of operators to sample each iteration

}; //PARAMS

PARAMS::PARAMS(){
    //initializes commonly used parameters from a file
    ifstream pfin;
    pfin.open("param.dat");

    pfin >> nX_;
    pfin >> NN_;
    pfin >> sample_;
    pfin >> EQL_;
    pfin >> MCS_;
    pfin >> SEED_;
    pfin.close();

    //initialize constants : 
	//1D linear
	LinX = nX_;
	numSpin = 2*nX_;
	numVB = numSpin/2;
	numLattB = numSpin;
    //2D square
  // LinX = nX_;
  // numSpin = nX_ * nX_;
  // numVB= numSpin/2;
  // numLattB= 2*numSpin;

	//initialize lattice bond array: 1D linear
	int a,b;
	index2 temp;
	for (int i=0; i<numSpin; i++){
		a = i;
		b = i+1;
		if ( b%nX_ == 0) b -= nX_;
		if (i%2 == 0) temp.set(a,b);  //order (A,B) sublattice
	    else temp.set(b,a);
		Bst.push_back(temp);

	}

   // //initialize lattice bond array: 2D square
   // int a, b;
   // int x,y;
   // index2 temp;
   // for (int i=0; i<numSpin; i++){
   //     //horizontal bond
   //     a = i;
   //     b = i+1;
   //     if ( b%nX_ == 0) b -= nX_;
   //     x = a%nX_; y = a/nX_;
   //     if ((x+y)%2==0) temp.set(a,b);  //order (A,B) sublattice
   //     else temp.set(b,a);
   //     Bst.push_back(temp);
   //     //vertical bond
   //     a = i;
   //     b = i+nX_;
   //     if (b>= numSpin) b -= numSpin;
   //     x = a%nX_; y = a/nX_;
   //     if ((x+y)%2==0) temp.set(a,b);  //order (A,B) sublattice
   //     else temp.set(b,a);
   //     Bst.push_back(temp);
   // }

}//constructor

void PARAMS::printBst(){

    for (int i=0; i<Bst.size(); i++)
        Bst[i].print();

}

#endif
