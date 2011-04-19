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
        int ratio_;     //the number for Renyi in the denominator of SWAP
        long SEED_;
        PARAMS();

        //lattice spin coordination numbers
        vector <index2> Bst;

        //lattice plaquette coordination numbers (square)
        vector <index4> Pst;

        void printBst();
        void printPst();
		
        //derived constants
        int numVB; //number of valence bonds = N/2

}; //PARAMS

PARAMS::PARAMS(){
    //initializes commonly used parameters from a file
    ifstream pfin;
    pfin.open("param.dat");

    pfin >> nX_;
    pfin >> EQL_;
    pfin >> MCS_;
    pfin >> nBin_;
    pfin >> SEED_;
    pfin >> ratio_;
    pfin.close();

//	NN_ *= 2*nX_*nX_;  // multiply by number of spins (in 2D)


    //2D square
    LinX = nX_;
    //numSpin = 2* nX_ * nX_;
    numSpin = nX_ * nX_;
    numVB= numSpin/2;
    numLattB= 2*numSpin;
    //initialize lattice bond array: 2D square
    int a, b, c, d;
    int x,y;
    index2 temp;
    index4 temp4;
    for (int i=0; i<numSpin; i++){  //REAL SYSTEM (LAYER 1)
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
        d = i+nX_;
        if (d>= LinX*LinX) d -= LinX*LinX;
        x = a%nX_; y = a/nX_;
        if ((x+y)%2==0) temp.set(a,d);  //order (A,B) sublattice
        else temp.set(d,a);
        Bst.push_back(temp);

        //diagonal bond (for Pst)
        c = i+nX_ + 1;
        if ( c%nX_ == 0) c -= nX_;
        if (c>= LinX*LinX) c -= LinX*LinX;
        temp4.set(a,b,c,d);
        Pst.push_back(temp4);

    }//layer 1


}//constructor

void PARAMS::printBst(){
    cout<<"Bst "<<endl;
    for (int i=0; i<Bst.size(); i++)
        Bst[i].print();
}

void PARAMS::printPst(){
    cout<<"Pst "<<endl;
    for (int i=0; i<Pst.size(); i++)
        Pst[i].print();
}

#endif
