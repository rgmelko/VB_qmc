// Oct 2, 2009 // New double-projector vb qmc program

#include"header.h"
#include"header_dubproj.h"
#include<iomanip>
using namespace std;

class LATTICE
{

public:
 
  MTRand drand; //drand() gives you a random double precision number
  MTRand_int32 irand; // irand() gives you a random integer

  int width, height, number_of_bondops, iterations;
  int number_of_nnbonds, number_of_bonds, number_of_sites;
  int change_number;  //number of bondops changed per step

