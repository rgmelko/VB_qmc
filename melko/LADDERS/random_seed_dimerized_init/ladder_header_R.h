#ifndef ladder_header_R
#define ladder_header_R

#include <vector>
#include "header.h"
#include "matrix.h"
using namespace std;

class LADDER 
{

 public:

  MTRand drand; //drand() gives you a random double precision number
  MTRand_int32 irand; // irand() gives you a random integer

  int legs, length, number_of_bondops, iterations;
  int number_of_nnbonds;
  int number_of_bonds, number_of_sites;
  int change_number; // number of bond ops changed per step
  
  int offdiagA, offdiagB;/*number of offdiag bond ops for
				 system A and B */
  long long enercounter; //counter of nn bonds for energy
  long double accept;
  long double energy;
  string bondfile;
  
  vector <int> bonds; // the bonds
  vector <int> init; // the inital bond state
  vector <int> bondopsA, bondopsB; // bondops for systems A and B
  vector <long long> entrocounter; /* the entropies where element 0
				      is for zone size 1, etc. */
  vector <long double> entropies;

  iMatrix nnbonds; //list of all nnbonds (used to pick bond ops)
  iMatrix nncheck; //matrix st nncheck(a,b) is 1 if a,b are nn. O otherwise
 
  // Constructor
  LADDER(int legs,int length,int number_of_bondops, int change_number, 
	 int iterations, string bondopfilename, long long randnumbseed);
  
  void nnbondlist();//creates list of nnbonds & nncheck
  void generate_ops();//generates the initial operators
  void apply_ops();//applies ops
  void change_ops();//changes a certain number of operators
  void decide();//decides to keep the new ops or go back to the old ones
  void measure();//measures energy and entropy after a step
  void reinitialize();//reinitializes bonds, offdiag etc for the next step
  void first_step();
  void calculate_stuff();
  void read_bonds();
  void super_initialize();
};

LADDER::LADDER(int a,int b,int c, int d, int e, string f, long long g)
{
  irand.seed(g);
  drand.seed(g);

  legs = a;
  length = b;
  number_of_bondops = c;
  change_number = d;
  iterations = e;
  bondfile = f;
  number_of_nnbonds = 2*a*b - a - b;
  number_of_sites = a*b;
  number_of_bonds = number_of_sites/2;
  offdiagA = 0;
  offdiagB = 0;
  enercounter = 0;
  accept = 0;

  nnbonds.resize (number_of_nnbonds,2);
  nncheck.resize (number_of_sites,number_of_sites);
  bonds.resize (number_of_sites);
  init.resize (number_of_sites);
  bondopsA.resize (number_of_bondops); 
  bondopsB.resize (number_of_bondops);
  entrocounter.resize (number_of_sites-1);
  entropies.resize (number_of_sites-1);

  for(int i09=0; i09<number_of_sites; i09++)
    {
      for(int i10=0; i10<number_of_sites; i10++)
	{
	  nncheck(i09,i10) = 0;
	}
    }

  entrocounter.assign(number_of_sites-1,0);
  
  for(int i02=0; i02<number_of_sites; i02+=2)
    {
      init[i02] = i02+1;
      init[i02+1] = i02;
    }

 
  //sets the bonds to the initial state
  bonds = init;
}

void LADDER::nnbondlist()
{
  int bondnum = 0;
  
  //the first bonds are of the form (0,1),(1,2),(2,3) etc
  for(bondnum; bondnum < legs*length-1; bondnum++)
    {
      nnbonds(bondnum,0) = bondnum;
      nnbonds(bondnum,1) = bondnum+1;

      nncheck(bondnum, bondnum+1) = 1;
      nncheck(bondnum+1, bondnum) = 1;
    }

  //the rest are more complicated
  while(bondnum < number_of_nnbonds)
    {
      int sitenum = 0;
      for(int legcounter=legs*2-1; legcounter>1; legcounter-=2)
	{
	  for(int i00 = sitenum; i00<(length-1)*legs; i00+=legs)
	    {
	      nnbonds(bondnum,0)= i00;
	      nnbonds(bondnum,1)= i00+legcounter;

	      nncheck(i00, i00+legcounter) = 1;
	      nncheck(i00+legcounter, i00) = 1;

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
      bondopsA[i01] = irand() % number_of_nnbonds;
      bondopsB = bondopsA;
    }
}

void LADDER::change_ops()
{
  for(int i04=0; i04<change_number; i04++)
    {
      bondopsB[irand()%number_of_bondops] = irand()%number_of_nnbonds;
    }
}

void LADDER::apply_ops()
{
  int a(0),b(0),c(0),d(0);

  for(int i04=0; i04<number_of_bondops; i04++)
    {      
      a = nnbonds(bondopsB[i04],0);
      b = nnbonds(bondopsB[i04],1);

      //cout << "---------------- " <<a << "," << b <<  " -----------------"<< endl;

      if(bonds[a] != b)
	{
	  c = bonds[a];
	  d = bonds[b];
	  
	  bonds[c] = d;
	  bonds[d] = c;
	  bonds[a] = b;
	  bonds[b] = a;
	  
	  offdiagB++;
	}
      /*
	for(int i=0; i< bonds.size(); i++)
	{
	cout << i << "," << bonds[i] << endl;
	}
	cout << endl;
      */
    }
}

void LADDER::measure()
{
  for(int i05=0; i05<number_of_sites; i05++)
    {
      if(i05<bonds[i05])
	{
	  enercounter += nncheck(i05,bonds[i05]);

	  for(int i06=i05; i06<bonds[i05]; i06++)
	    {
	      entrocounter[i06]++;
	    }
	}
    }
   //cout << enercounter << endl;
}


/*given the system (A or B) that was kept last, this function decides whether
  the the new system is kept, or if we go back to the old system*/
void LADDER::decide()
{
  double prob = pow(2,offdiagA-offdiagB);
  if(drand()<prob)
    {
      bondopsA = bondopsB;
      offdiagA = offdiagB;
      
      accept++;
    }
    
}

void LADDER::reinitialize()
{
  bonds = init;
  offdiagB = 0;
  bondopsB = bondopsA;
  //cout << "STEP COMPLETED!!!!!!!!" << endl;
}

void LADDER::first_step()
{
  generate_ops();
  apply_ops();
  offdiagA = offdiagB;
  bondopsA = bondopsB;
  reinitialize();
}

void LADDER::calculate_stuff()
{
  // cout << enercounter << endl;
  energy = (-0.5*number_of_nnbonds)*((long double)enercounter/((long double)number_of_nnbonds*
							       ((long double)iterations+1)) + 0.5);

  for(int i07=0; i07<number_of_sites-1; i07++)
    {
      entropies[i07] = entrocounter[i07]*log(2)/(iterations);
    }
}

void LADDER::read_bonds()
{
  ifstream fin2(bondfile.c_str());
  
  for(int i11=0; i11<number_of_bondops; i11++)
    {
      fin2 >> bondopsA[i11];
    }
  fin2 >> offdiagA;
  fin2.close();

  bondopsB = bondopsA;
}

void LADDER::super_initialize()
{
  offdiagA = 0;
  offdiagB = 0;
  enercounter = 0;
  accept = 0;
  entrocounter.assign(number_of_sites-1,0);
 }

#endif
