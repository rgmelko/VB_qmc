// 
// RVB_oop.cpp: the main body of the program
// A OOP VB basis projector code for the RVB wavefunction
// Roger Melko, April 2011, Santa Barbara Station Q
//
#include "head_proj.h"
#include "simparam.h"
#include "basis.h"
//#include "measure.h"
//#include "renyi.h"

int main(){

    PARAMS param;
    MTRand mrand(param.SEED_); //random number for metropolis

    Basis Vbeta(param); 
    Basis Valpha(param); 

    Vbeta.print();
    Valpha.print();

    int temp;
    temp = Vbeta|Valpha;
    cout<<temp<<endl;




	return 0;
};
