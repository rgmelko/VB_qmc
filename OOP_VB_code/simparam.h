#ifndef PARAMSIM_H
#define PARAMSIM_H

#include <fstream>
#include "head_proj.h"

//Class to read in the simulation parameters from a file
class PARAMS
{
    public:
        int numSpin; 
        int numLattB;  //number of lattice bonds is 2N
        long SEED_;
        PARAMS();

    protected:

        //derived constants
        int numVB; //number of valence bonds = N/2

        //lattice spin coordination numbers
        vector <index2> Bst;

        long int NN_;  //This is the length of the operator string
        int nX_;     //linear size of square lattice
        int sample_; //number of operators to sample each iteration
        void printBst();

}; //PARAMS

PARAMS::PARAMS(){
    //initializes commonly used parameters from a file
    ifstream pfin;
    pfin.open("param.dat");

    pfin >> nX_;
    pfin >> NN_;
    pfin >> sample_;
    pfin >> SEED_;
    pfin.close();

    //initialize constants
    numSpin = nX_ * nX_;
    numVB= numSpin/2;
    numLattB= 2*numSpin;

    //initialize lattice bond array
    int a, b;
    int x,y;
    index2 temp;
    for (int i=0; i<numSpin; i++){

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
        if (b>= numSpin) b -= numSpin;
        x = a%nX_; y = a/nX_;
        if ((x+y)%2==0) temp.set(a,b);  //order (A,B) sublattice
        else temp.set(b,a);
        Bst.push_back(temp);

    }

}//constructor

void PARAMS::printBst(){

    for (int i=0; i<Bst.size(); i++)
        Bst[i].print();

}

#endif
