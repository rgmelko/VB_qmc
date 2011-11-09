#ifndef PARAMSIM_H
#define PARAMSIM_H

#include <fstream>
#include "head_proj.h"

//Class to read in the simulation parameters from a file
class PARAMS
{
    public:
        int numSpin; 
        int numLattB; //number of lattice bonds is 2N
		int EQL_;     //the number of Monte Carlo steps
		int MCS_;     //the number of Monte Carlo steps
        int nX_;      //linear size of lattice
        int nY_;      //linear size of lattice
        int nBin_;    //number of production bins
        long SEED_;

        long m_;     //The order of the expansion
        double h_x;  //Transverse field 

        PARAMS();

        //lattice spin coordination numbers
        vector <index2> Bst;

        void printBst();
   
}; //PARAMS

PARAMS::PARAMS(){
    //initializes commonly used parameters from a file
    ifstream pfin;
    pfin.open("param.dat");

    pfin >> nX_;
    pfin >> nY_;
    pfin >> EQL_;
    pfin >> MCS_;
    pfin >> nBin_;
    pfin >> SEED_;
    pfin >> m_;
    pfin >> h_x;
    pfin.close();


    if (nY_ == 0) {    // ---------1D chain
        //derived constants
        numSpin = nX_;
        numLattB = nX_-1; //Periodic BC

        //Initialize lattice bond array
        int a,b;
        index2 temp;
        for (int i=0; i<numLattB; i++){  //REAL SYSTEM (LAYER 1)
            a = i;
            b = i+1;
            if (b == numSpin) b = 0; //PBC
            temp.set(a,b);
            Bst.push_back(temp);
        }//i

    }//1D chain
    else{

        numSpin = nX_*nX_;
        numLattB = 2*nX_*nX_; //Periodic BC for 2D lattice

        //Initialize lattice bond array
        int a,b,d;
        index2 temp;

        for (int i=0; i<numSpin; i++){  
            //horizontal bond
            a = i;
            b = i+1;
            if ( b%nX_ == 0) b -= nX_;
            temp.set(a,b);
            Bst.push_back(temp);
            //vertical bond
            a = i;
            d = i+nX_;
            if (d>= nX_*nX_) d -= nX_*nX_;
            temp.set(a,d);
            Bst.push_back(temp);
        }//i

    }



}//constructor

void PARAMS::printBst(){
    cout<<"Bst "<<endl;
    for (int i=0; i<Bst.size(); i++)
        Bst[i].print();
}

#endif
