// Nov 4, 2009 // Header for new double-projector vb qmc prog (dubproj.cpp)
// Nov 17,2009 // Fixed energy measurement so it doesn't take so long
// Nov 18,2009 // Putting in bond length measurement along (x,0)
// Nov 20,2009 // Trying to adapt this prog to measure renyi ent

#ifndef header_renyi
#define header_renyi

#include"header.h"
#include"matrix.h"
#include<vector>

using namespace std;

class LATTICE
{

public:
  
  MTRand drand; //drand() gives you a random double precision number
  MTRand_int32 irand; // irand() gives you a random integer

  int xsites, ysites, zsites, dim, number_of_bondops, iterations;
  int number_of_nnbonds;
  int number_of_bonds, number_of_sites;
  int change_number; // number of bond ops changed per step
 
  int offdiag1, offdiag2, offdiag3;
  int newloops, oldloops;
  long double accept;
  long double energy;
  long int energyint, mdiff;
  string bondfile1, bondfile2;
  
  vector <int> bonds1, bonds2, bonds3, nrgbonds; // the bonds
  vector <int> init; // the inital bond state
  vector <int> bondops1, bondops2, bondops3; // bondops for systems A and B
  vector <int> whichloop_new, whichloop_old; // stores which loop each site is in. used for nrg
  vector <long double> entropy;
  
 
  iMatrix nnbonds; //list of all nnbonds (used to pick bond ops)
 
  // Constructor
  LATTICE(int xsites, int ysites, int zsites,
	  int number_of_bondops, int change_number, int iterations, 
	  string bondopfile1, string bondfile2, long long randnumbseed);

  void nnbondlist();//creates list of nnbonds & nncheck
  void generate_ops();//generates the initial operators
  void apply_ops(int stepnum);//applies ops
  void change_ops(int stepnum);//changes a certain number of operators
  void decide(int stepnum);//decides to keep the new ops or go back to the old ones
  void measure();//measures energy and entropy after a step
  void reinitialize(int stepnum);//reinitializes bonds, offdiag etc for the next step
  void first_step();
  void calculate_stuff();
  void read_bonds();
  void super_initialize(); 
  int loop_counter();
  void measure_energy();
  void print_quantities(string filename, long double quantity);
  void print_bondops(string filename, vector <int> bondlist, int offdiag);
  void print_entropies(string filename, vector <long double> entropy);
  void swaperator();
  
};

LATTICE::LATTICE(int x, int y, int z, int bondops, int r, int its, 
		 string bondopfile1, string bondopfile2, long long rseed)
{
  irand.seed(rseed);
  drand.seed(rseed);

  xsites = x;
  ysites = y;
  zsites = z;
  number_of_bondops = 2*bondops;
  change_number = r;
  iterations = its;
  bondfile1 = bondopfile1;
  bondfile2 = bondopfile2;
  number_of_sites = 2*x*y*z;
  number_of_bonds = number_of_sites/2;
  accept = 0;
  offdiag1=0; offdiag2=0; offdiag3=0;
  newloops=0; oldloops=0;
  energy=0, energyint=0, mdiff=0;

  //determining dimension of the system
  if((x==1&&y==1)|(x==1&&z==1)|(y==1&&z==1)){dim=1;}
  else if((x==1)|(y==1)|(z==1)){dim=2;}
  else {dim=3;}

  number_of_nnbonds = dim*number_of_sites; //depends on dimension

  nnbonds.resize (number_of_nnbonds,2);
  bonds1.resize (number_of_sites);
  bonds2.resize (number_of_sites);
  bonds3.resize (number_of_sites);
  init.resize (number_of_sites);
  bondops1.resize (number_of_bondops); 
  bondops2.resize (number_of_bondops);
  bondops3.resize (number_of_bondops);
  whichloop_new.resize (number_of_sites);
  whichloop_old.resize (number_of_sites);
  entropy.assign(xsites-1,0);
  
  //initial state is dimerized..
  for(int i01=0; i01<number_of_sites; i01+=2)
    {
      init[i01] = i01+1;
      init[i01+1] = i01;
    }
 
  //sets the bonds to the initial state
  bonds1 = init; bonds2 = init; bonds3 = init;
}

/******************** BONDLIST ************************************/
void LATTICE::nnbondlist()
{
  int b1(0), b2(0), b3(0), b4(0);
  int counter = 0;

  for(int bz=0; bz < zsites; bz++){
    for(int by=0; by < ysites; by++){
      for(int bx=0; bx < xsites; bx++){
	b1 = bx + by*xsites + bz*xsites*ysites; 
	b2 = (bx+1)%xsites + by*xsites + bz*xsites*ysites;
	b3 = bx + ((by+1)%ysites)*xsites + bz*xsites*ysites;
	b4 = bx + by*xsites + ((bz+1)%zsites)*xsites*ysites;
	
	if(b1!=b2){
	  nnbonds(counter,0) = b1;
	  nnbonds(counter,1) = b2;
	  counter++;
	}
	if(b1!=b3){
	  nnbonds(counter,0) = b1;
	  nnbonds(counter,1) = b3;
	  counter++;
	}
	if(b1!=b4){
	  nnbonds(counter,0) = b1;
	  nnbonds(counter,1) = b4;
	  counter++;
	}
      }
    }
  }
  for(int b=number_of_nnbonds/2; b<number_of_nnbonds; b++){
    int a = number_of_nnbonds/2;
    nnbonds(b,0) = nnbonds(b-a,0)+a;
    nnbonds(b,1) = nnbonds(b-a,1)+a;
  }
}

/*********** GENERATE OPERATORS *****************************************/
void LATTICE::generate_ops()
{
  for(int i=0; i<number_of_bondops; i++)
    {
      bondops1[i] = irand() % number_of_nnbonds;
      bondops2[i] = irand() % number_of_nnbonds;
      bondops3 = bondops1;
    }
}

/*********** APPLY OPERATORS ********************************************/
void LATTICE::apply_ops(int stepnum)
{
  int a(0),b(0),c(0),d(0);

  if(stepnum%2==0){
    for(int i04=0; i04<number_of_bondops; i04++){ 
      a = nnbonds(bondops1[i04],0);
      b = nnbonds(bondops1[i04],1);
      
      if(bonds1[a] != b){
	c = bonds1[a];
	d = bonds1[b];
	
	bonds1[c] = d;
	bonds1[d] = c;
	bonds1[a] = b;
	bonds1[b] = a;
	
	offdiag1++;
      }
    }
  }
  else{
    for(int i04=0; i04<number_of_bondops; i04++){      
      a = nnbonds(bondops2[i04],0);
      b = nnbonds(bondops2[i04],1);
      
      if(bonds2[a] != b){
	c = bonds2[a];
	d = bonds2[b];
	
	bonds2[c] = d;
	bonds2[d] = c;
	bonds2[a] = b;
	bonds2[b] = a;
	  
	offdiag2++;	
      }
    }
  }
  newloops = loop_counter();
}

/************* CHANGE OPERATORS ****************************************/
void LATTICE::change_ops(int stepnum)
{
  if(stepnum%2==0){
    for(int i04=0; i04<change_number; i04++){
      bondops1[irand()%number_of_bondops] = irand()%number_of_nnbonds;
    }
  }
  else{ 
    for(int i04=0; i04<change_number; i04++){
      bondops2[irand()%number_of_bondops] = irand()%number_of_nnbonds;
    }
  }
}

/************* REINITIALIZE ********************************************/
void LATTICE::reinitialize(int stepnum)
{
  newloops = 0;

  if(stepnum%2==0){
    bondops3 = bondops2;
    bonds3 = bonds2;
    bonds2 = init;
    offdiag3 = offdiag2;
    offdiag2 = 0;
    
  }
  else{
    bondops3=bondops1;
    bonds3 = bonds1;
    bonds1 = init;
    offdiag3 = offdiag1;
    offdiag1 = 0;
  }
}

/************ FIRST STEP ************************************************/
void LATTICE::first_step()
{
  generate_ops();
  apply_ops(1);
  apply_ops(2);

  oldloops = newloops;
  whichloop_old = whichloop_new;
  reinitialize(1);
}

/************ READ BOND FILES  *******************************************/
void LATTICE::read_bonds()
{
  ifstream fin2(bondfile1.c_str());
  for(int i11=0; i11<number_of_bondops; i11++){
    fin2 >> bondops1[i11];
  }
  fin2 >> offdiag1;
  fin2.close();
  
  ifstream fin9(bondfile2.c_str());
  for(int i11=0; i11<number_of_bondops; i11++){
    fin9 >> bondops2[i11];
  }
  fin9 >> offdiag2;
  fin9.close();

  apply_ops(1);
  apply_ops(2);

  oldloops = newloops;
  whichloop_old = whichloop_new;
  reinitialize(1);
}

/************ LOOP COUNTER ***********************************************/
int LATTICE::loop_counter()
{
  int counter(0), loopnumber(0), startsite(0), a(-99), which(0);
  vector <int> site(number_of_sites+2,0);

  while(counter < number_of_sites){

    site[counter]=1;
    startsite = counter;
    
    whichloop_new[startsite] = loopnumber+1;

    a = bonds1[counter];
    which = 0;
   
    while(a!=startsite){
      
      site[a] = 1;
      whichloop_new[a] = loopnumber+1;

      if(which==0){
	a = bonds2[a];
	which++;
      }
      else{
	a = bonds1[a];
	which--;
      }
      whichloop_new[a] = loopnumber+1;
    }
    loopnumber++;
    while(site[counter]==1){counter++;}
  }
  return loopnumber;
}


//************ DECIDE TO KEEP OR REVERT ****************************/
void LATTICE::decide(int stepnum)
{
  double prob = 0;
  double weight = 0;
  double brand = drand();

  int loopyloop = newloops - oldloops;
  prob = pow(2,loopyloop);

  if(stepnum%2==0){
    int moffbag = offdiag3-offdiag1;
    weight = pow(2,moffbag);
    prob = weight*prob;
    
    if(brand<prob){
      oldloops = newloops;
      whichloop_old = whichloop_new;
      accept++;
    }
    else{
      bondops1 = bondops3;
      bonds1 = bonds3;
      offdiag1 = offdiag3;
    }
  }
  else{
    int moffbag = offdiag3-offdiag2;
    weight = pow(2,moffbag);
    prob = weight*prob;

    if(brand<prob){
      oldloops = newloops;
      whichloop_old = whichloop_new;
      accept++;
    }
    else{
      bondops2 = bondops3;
      bonds2 = bonds3;
      offdiag2 = offdiag3;
    }
  }   
}

/****************** SUPER INITIALIZE *******************************/
void LATTICE::super_initialize()
{
  offdiag1 = 0;
  offdiag2 = 0;
  offdiag3 = 0;
  accept = 0;
  energy = 0;
  energyint = 0;
  entropy.assign(xsites-1,0);
 }

/********************* ENERGY MEASUREMENT **************************/
 void LATTICE::measure_energy()
{
  mdiff=0;

  for(int k1=0; k1<number_of_nnbonds; k1++){

    int a(0),b(0);
    a = nnbonds(k1,0);
    b = nnbonds(k1,1);
    
    if(whichloop_old[a]!=whichloop_old[b]){mdiff++;}
  }
  energyint += mdiff - number_of_nnbonds;
}


/******************** ENERGY/ENTROPY CALCULATIONS ****************/
void LATTICE::calculate_stuff()
{
  energy = energyint*0.75/iterations;

  for(int rint=0; rint<xsites-1; rint++){
    entropy[rint]/=(1.0*iterations);
    entropy[rint] = -log(entropy[rint]);
  }
}

/***************** PRINT ENERGIES / CORRELATIONS ***********************/
void LATTICE::print_quantities(string filename, long double quantity)
{
  ofstream fouttemp(filename.c_str(),ios::app);
  fouttemp.precision(10);
  fouttemp << quantity << endl;
  fouttemp.close();
}

void LATTICE::print_entropies(string filename, vector <long double> entropy)
{
  ofstream foutent(filename.c_str(),ios::app);
  foutent.precision(10);
  foutent.width(12);
  for(int sint=0; sint<xsites-1; sint++){
    foutent << entropy[sint] << " ";
  }
  foutent << endl;
  foutent.close();
}

/******************* PRINT BONDOP LIST *********************************/
void LATTICE::print_bondops(string filename, vector <int> bondlist, int offdiag)
{
  ofstream foutbond(filename.c_str());
  for(int mint=0; mint<number_of_bondops; mint++){
    foutbond << bondlist[mint] << endl;
  }
  foutbond << offdiag << endl;
  foutbond.close();
}

void LATTICE::swaperator()
{
  //figure out a way to define the zone, or many zones.
  // maybe coordinates in the param file

  //swap bonds within zone for top & bottom systems

  //do interproduct ie count number of loops for swapped system

  //entropy = -ln(<swap>) = -ln(<bonds1|swap|bonds2>/<bonds1|bonds2>)

  //bam

  //extra credit: figure out if there's a way to use old loop#s
  //              to determine loop # for larger zone sizes....

  vector <int> tempbonds, temploops;
  tempbonds = bonds2;
  temploops = whichloop_old;
  int temploopnum = oldloops;
  int loopnummax = oldloops;
  int a,b,c,d;
  int ollld, neeew;


  for(int lint=0; lint<xsites-1; lint++){
    a = lint;
    d = lint+xsites;
    b = tempbonds[d];
    c = tempbonds[a];

    tempbonds[c] = d;
    tempbonds[d] = c;
    tempbonds[a] = b;
    tempbonds[b] = a;

    if(temploops[a]==temploops[b]){  
      temploopnum++;                    
      temploops[a]=loopnummax;
      temploops[b]=loopnummax;
      ollld = bonds1[b];
      while(ollld!=a){    
	temploops[ollld]=loopnummax;
	ollld=tempbonds[ollld];
	temploops[ollld]=loopnummax;
	ollld = bonds1[ollld];
	}
      loopnummax++;
    }
    else{                    //If sites are in a different loop, number of loops decreases.
      temploopnum--;          // Relabel loop 2 so it's the same as loop 1.
      ollld=temploops[b];
      neeew=temploops[a];
      for(int nint=0; nint<number_of_sites; nint++){
	if(temploops[nint]==ollld){temploops[nint]=neeew;}
      } 
    }
    int loopdiff = temploopnum - oldloops;
    entropy[lint] += pow(2,loopdiff);
  }
}

#endif


