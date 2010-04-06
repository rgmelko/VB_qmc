// Nov 4, 2009 // Header for new double-projector vb qmc prog (dubproj.cpp)
// Nov 17,2009 // Fixed energy measurement so it doesn't take so long
// Nov 18,2009 // Putting in bond length measurement along (x,0)
// Nov 20,2009 // Trying to adapt this prog to measure renyi ent
// Apr 1, 2010 // adapting to use ratio trick

#ifndef header_renyi
#define header_renyi

#include"header.h"
#include"matrix.h"
#include<vector>

class LATTICE
{

public:
  
  MTRand drand; //drand() gives you a random double precision number
  MTRand_int32 irand; // irand() gives you a random integer

  int xsites, ysites, zsites, dim, number_of_bondops, iterations;
  int number_of_nnbonds;
  int number_of_bonds, number_of_sites;
  int change_number; // number of bond ops changed per step
  int its_per_meas;
 
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
	  string bondopfile1, string bondfile2, int its_per_measurement,
	  long long randnumbseed);

  void nnbondlist();//creates list of nnbonds
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
  int  loop_counter(vector <int> leftside, vector <int> rightside, int superwhich);
  void measure_energy();
  void print_quantities(string filename, long double quantity);
  void print_bondops(string filename, vector <int> bondlist, int offdiag);
  void print_entropies(string filename, vector <long double> entropy);
  double swaperator(int linsize);
  void measure_swap();
  
};

LATTICE::LATTICE(int x, int y, int z, int bondops, int r, int its, 
		 string bondopfile1, string bondopfile2, int its_per_measurement,
		 long long rseed)
{
  irand.seed(rseed);
  drand.seed(rseed);

  xsites = x;
  ysites = y;
  zsites = z;
  number_of_bondops = 2*bondops;
  change_number = r;
  iterations = its;
  its_per_meas = its_per_measurement;
  bondfile1 = bondopfile1;
  bondfile2 = bondopfile2;
  number_of_sites = 2*xsites*ysites;
  number_of_bonds = xsites*ysites;
  accept = 0;
  offdiag1=0; offdiag2=0; offdiag3=0;
  newloops=0; oldloops=0;
  energy=0, energyint=0, mdiff=0;

  //determining dimension of the system
  if((x==1&&y==1)|(x==1&&z==1)|(y==1&&z==1)){dim=1;}
  else if((x==1)|(y==1)|(z==1)){dim=2;}
  else {dim=3;}

  number_of_nnbonds = 4*xsites*ysites; //depends on dimension

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
  int L = xsites;

  for(int i=0; i<L; i++){
    for(int j=0; j<L; j++){
      nnbonds(i+j*L,0) = i + j*L;
      nnbonds(i+j*L,1) = j*L + (i+1)%L;
    }
  }

  for(int i=0; i<L; i++){
    for(int j=0; j<L; j++){
      nnbonds(L*L+i+j*L,0) = i + j*L;
      nnbonds(L*L+i+j*L,1) = i + ((j+1)%L)*L;
    }
  }


  for(int b1=2*L*L; b1<4*L*L; b1++){
    int a1 = 2*L*L;
    nnbonds(b1,0) = nnbonds(b1-a1,0)+L*L;
    nnbonds(b1,1) = nnbonds(b1-a1,1)+L*L;
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

  // for(int q =0; q<number_of_sites ; q++){ cout << q << "," << bonds1[q] << "     " << q << "," <<  bonds2[q] << endl;}
  // cout << endl;

  newloops = loop_counter(bonds1, bonds2, 1);
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
  offdiag1=0;
  fin2.close();
  
  ifstream fin9(bondfile2.c_str());
  for(int i11=0; i11<number_of_bondops; i11++){
    fin9 >> bondops2[i11];
  }
  offdiag2=0;
  fin9.close();

  apply_ops(1);
  apply_ops(2);

  oldloops = newloops;
  whichloop_old = whichloop_new;
  reinitialize(1);
}

/************ LOOP COUNTER ***********************************************/
int LATTICE::loop_counter(vector <int> leftside, vector <int> rightside, int superwhich)
{
  int counter(0), loopnumber(0), startsite(0), AA(-99), which(0);
  vector <int> site(number_of_sites+2,0);

  while(counter < number_of_sites){

    site[counter]=1;
    startsite = counter;
    
    if(superwhich){whichloop_new[startsite] = loopnumber+1;}

    AA = leftside[counter];
    which = 0;
   
    while(AA!=startsite){
      site[AA] = 1;
      if(superwhich){whichloop_new[AA] = loopnumber+1;}
      
      if(which==0){
	AA = rightside[AA];
	which++;
      }
      else{
	AA = leftside[AA];
	which--;
      }
      whichloop_new[AA] = loopnumber+1;
      site[AA] = 1;
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

    int a3(0),b3(0);
    a3 = nnbonds(k1,0);
    b3 = nnbonds(k1,1);
    
    if(whichloop_old[a3]!=whichloop_old[b3]){mdiff++;}
  }
  energyint += mdiff - number_of_nnbonds;
}


/******************** ENERGY/ENTROPY CALCULATIONS ****************/
void LATTICE::calculate_stuff()
{
  energy = energyint*0.75/iterations*its_per_meas;

  for(int rint=0; rint<xsites-1; rint++){
    entropy[rint]/=(1.0*iterations/its_per_meas);
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

/******************* SWAPERATOR ****************************************/

double LATTICE::swaperator(int linsize)
{  

  vector <int> tempbonds;
  tempbonds = bonds2;

  vector <int> VLflip = bonds1, VLuse = bonds1;
  int a,b,c,d;
 
  for(int A=0; A<linsize; A++){
    for(int B=0; B<linsize; B++){
    	  
    a = A + B*xsites;
    d = a + number_of_sites/2;
    b = tempbonds[d];
    c = tempbonds[a];
    
    tempbonds[a] = b;
    tempbonds[b] = a;
    tempbonds[d] = c;
    tempbonds[c] = d;

    }
  }
  
  int temploopnum = loop_counter(VLuse, tempbonds, 0);
  int loopdiff = temploopnum - oldloops;
  return pow(2,loopdiff);
}

void LATTICE::measure_swap()
{  
  for(int i=0; i<xsites-1; i++){
    entropy[i] += swaperator(i+1);
  }
}
#endif


