#include "GenHam.h"

// 4x4 SQUARE lattice with periodic boundary conditions
//----------------------------------------------------------
void GENHAM::Bonds_16B(){
//   This is all you need for the Heisenberg model
//   Each row indexes two sites associated with ll
//   
//     4
//     |           e.g. site 0 related to site 1 and 4
//   ( 0 ) - 1

  Bond.resize(Nsite,4);

  //     (redundant) index, horizontal, vertical, special (anisotripic) flag

  Bond = 0, 1, 4,     0,
         1, 2, 5,     0,
         2, 3, 6,     0,
         3, 0, 7,     0,
         4, 5, 8,     0,
         5, 6, 9,     0,
         6, 7, 10,    0,
         7, 4, 11,    0,
         8, 9, 12,    0,
         9, 10, 13,   0,
         10, 11, 14,  0,
         11, 8, 15,   0,
         12, 13, 0,   0,
         13, 14, 1,   0,
         14, 15, 2,   0,
         15, 12, 3,   0;

// If you wish to remove any of the bonds, simply make the 2nd or 3rd column above -99
//
// To add anisotropic interactions, read the GenHam code

  
}//MakeBonds

