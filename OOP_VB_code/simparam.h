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
        int nX_;     //linear size of square lattice
        int nBin_;     //number of production bins
        long SEED_;
        PARAMS();

        //lattice spin coordination numbers
        vector <index2> Bst;

        void printBst();
		
        //derived constants
        int numVB; //number of valence bonds = N/2
        long int NN_;  //This is the length of the operator string
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
    pfin >> nBin_;
    pfin.close();

	sample_ *= 2*nX_*nX_;  // multiply by number of spins (in 2D)

   //   //initialize constants : 1D linear
   //   LinX = nX_;
   //   numSpin = 2*nX_;
   //   numVB = numSpin/2;
   // 	//initialize lattice bond array: 1D linear
   // 	numLattB = 0;
   // 	int a,b;
   // 	index2 temp;
   // 	for (int i=0; i<numSpin; i++){
   // 		a = i;
   // 		b = i+1;
   // 		if ( b%nX_ == 0) b -= nX_;
   // 		else{
   // 		  if (i%2 == 0) temp.set(a,b);  //order (A,B) sublattice
   // 	      else temp.set(b,a);
   // 		  Bst.push_back(temp);
   // 		  numLattB ++;
   // 		}//OBC
   // 
   // 	}
   //  //cout<<numLattB<<endl;

    //2D square
    LinX = nX_;
    numSpin = 2* nX_ * nX_;
    numVB= numSpin/2;
    numLattB= 2*numSpin;
    //initialize lattice bond array: 2D square
    int a, b;
    int x,y;
    index2 temp;
    for (int i=0; i<numSpin/2; i++){  //REAL SYSTEM (LAYER 1)
        //horizontal bond
        a = i;
        b = i+1;
        if ( b%nX_ == 0) b -= nX_;
        x = a%nX_; y = a/nX_;
        if ((x+y)%2==0) temp.set(a,b);  //order (A,B) sublattice
        else temp.set(b,a);
        Bst.push_back(temp);
        //vertical bond
        a = i;
        b = i+nX_;
        if (b>= LinX*LinX) b -= LinX*LinX;
        x = a%nX_; y = a/nX_;
        if ((x+y)%2==0) temp.set(a,b);  //order (A,B) sublattice
        else temp.set(b,a);
        Bst.push_back(temp);
    }//layer 1
    for (int i=numSpin/2; i<numSpin; i++){  //ANCILLARY SYSTEM (LAYER 2)
        //horizontal bond
        a = i;
        b = i+1;
        if ( b%nX_ == 0) b -= nX_;
        x = a%nX_; y = a/nX_;
        if ((x+y)%2==0) temp.set(a,b);  //order (A,B) sublattice
        else temp.set(b,a);
        Bst.push_back(temp);
        //vertical bond
        a = i;
        b = i+nX_;
        if (b>= numSpin) b -= LinX*LinX;
        x = a%nX_; y = a/nX_;
        if ((x+y)%2==0) temp.set(a,b);  //order (A,B) sublattice
        else temp.set(b,a);
        Bst.push_back(temp);
    }//layer2

}//constructor

void PARAMS::printBst(){

    for (int i=0; i<Bst.size(); i++)
        Bst[i].print();

}

#endif
