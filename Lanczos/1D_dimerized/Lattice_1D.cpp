#include "GenHam.h"

// 1D dimerized square lattice with periodic boundary conditions
//----------------------------------------------------------
void GENHAM::Bonds_1D(){
//
//   0 - 1 = 2 - 3 = 0
//
//   -J,  +J;
//
  Bond.resize(Nsite,3);

  //     (redundant) index, horizontal, J=0 and J2=1

  for (int i=0; i<(Nsite-1); i++){
    Bond(i,0)=i;
    Bond(i,1)=i+1;
    if (i%2 == 0) Bond(i,2)=0;
    else Bond(i,2)=1; 
  }
  //PBC
  Bond(Nsite-1,0)=Nsite-1;
  Bond(Nsite-1,1)=0;   //FOR OBC USE Bond(Nsite-1,1)=-99;
  Bond(Nsite-1,2)=1;
  

//  cout<<Bond;

  
}//MakeBonds

