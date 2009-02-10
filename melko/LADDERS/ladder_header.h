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
  int offdiagA, offdiagB;
  int number_of_bonds, number_of_sites;
  int change_number;

  vector <int> bondsA, bondsB;
  vector <int> init;
  vector <int> nnbonds0, nnbonds1;
  vector <int> bondops;
 
  LADDER(int legs,int length,int number_of_bondops, int change_number);
  
  void nnbondlist();
  void generate_ops();
  void apply_ops(vector<int> bonds, int offdi);
  
};
#endif


LADDER::LADDER(int a,int b,int c, int d)
{
  legs = a;
  length = b;
  number_of_bondops = c;
  change_number = d;
  number_of_nnbonds = 2*a*b - a - b;
  number_of_sites = a*b;
  number_of_bonds = number_of_sites/2;
  offdiagA = 0;
  offdiagB = 0;

  nnbonds0.resize (number_of_nnbonds);
  nnbonds1.resize (number_of_nnbonds);
  bondsA.resize (number_of_sites);
  bondsB.resize (number_of_sites);
  init.resize (number_of_sites);
  bondops.resize (number_of_bondops);
 
  for(int i02=0; i02<number_of_sites; i02+=2)
    {
      init[i02] = i02+1;
      init[i02+1] = i02;
    }
  
  //sets the bonds to the initial state
  bondsA = init;
  bondsB = init;

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

void LADDER::apply_ops(vector<int> bonds, int offdi)
{
  int a(0),b(0),c(0),d(0);

  for(int i04=0; i04<number_of_bondops; i04++)
    {      
      a = nnbonds0[bondops[i04]];
      b = nnbonds1[bondops[i04]];
      if(bonds[a] == b)
	{}
      else
	{
	  c = bonds[a];
	  d = bonds[b];
	  
	  bonds[c] = d;
	  bonds[d] = c;
	  bonds[a] = b;
	  bonds[b] = a;

	  offdi++;
	  //add to offdiag or weight
	}
    }

}
