#ifndef ladder_header
#define ladder_header

#include <vector>
#include "mtrand.h"
using namespace std;

class LADDER 
{

 public:

  MTRand drand; //drand() gives you a random double precision number
  MTRand_int32 irand; // irand() gives you a random integer

  int legs, length, number_of_bondops;
  int number_of_nnbonds;
  int offdiag=0;
  int number_of_bonds;
  int change_number;

  vector <int> bonds0;
  vector <int> bonds1;
  vector <int> nnbonds0;
  vector <int> nnbonds1;
  vector <int> bondops;
  vector <bool> nn;

  LADDER(int legs,int length,int number_of_bondops, int change_number);
  
  void nnbondlist();
  void generate_ops();
  void apply_ops();
  
};
#endif


LADDER::LADDER(int a,int b,int c, int d)
{
  legs = a;
  length = b;
  number_of_bondops = c;
  change_number = d;
  number_of_nnbonds = 2*a*b - a - b;
  number_of_bonds = a*b/2;

  nnbonds0.resize (number_of_nnbonds);
  nnbonds1.resize (number_of_nnbonds);
  bonds0.resize (number_of_bonds);
  bonds1.resize (number_of_bonds);
  bondops.resize (number_of_bondops);
  nn.resize (number_of_bonds);
}

void LADDER::nnbondlist()
{
  int bondnum = 0;
  
  //the first bonds are of the form (0,1),(1,2),(2,3) etc
  for(bondnum; bondnum < legs*length-1; bondnum++)
    {
      nnbonds0[bondnum] = bondnum;
      nnbonds1[bondnum] = nnbonds0[bondnum]+1;
    }

  //the rest are more complicated
  while(bondnum < number_of_nnbonds)
    {
      int sitenum = 0;
      for(int legcounter=legs*2-1; legcounter>1; legcounter-=2)
	{
	  for(int i00 = sitenum; i00<(length-1)*legs; i00+=legs)
	    {
	      nnbonds0[bondnum]= i00;
	      nnbonds1[bondnum]= i00+legcounter;
	      bondnum++;
	    }
	  sitenum++;
	}
    }
}

void LADDER::generate_ops()
{
  for(int i01=0; i01<number_of_bondops; i01++)
    {
      bondops[i01] = irand() % number_of_nnbonds;
    }
}

void LADDER::apply_ops()
{


}
