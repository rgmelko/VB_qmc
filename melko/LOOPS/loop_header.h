//Jan 18, 2010 --- starting loop code

#ifndef loop_header
#define loop_header

#include "header.h"
#include "matrix.h"
#include <vector>

class LOOPS
{
 public:

  MTRand drand; //drand() gives you a random double precision number
  MTRand_int32 irand; // irand() gives you a random integer

  int dim1, dim2, number_of_bondops, number_of_sites, number_of_nnbonds,
    vlegs;
  bool OBC;

  vector <int> Vlinks, Hlinks; // the vertical and horizontal links for the LL
  vector <int> init_bonds; //the initial bonds
  vector <int> antipar, init_antipar, whichbonds, init_whichbonds; //keeping track of antiparallelness

  iMatrix nnbonds;//list of all possible nnbonds
  iMatrix nn_mat; /*matrix of the nnbonds. indices are sites
		    and contents are bond numbers.*/
  iMatrix Nnnbonds;//list of all nearest nnbonds
  iMatrix bops; // list of bond operators
  iMatrix superbops;

  //CONSTRUCTOR
  LOOPS::LOOPS(int xsites, int ysites, int bondops, bool ob, long long rseed);

  void nnbondlist();
  void Nnnbondlist();
  void generate_ops();
  void create_Vlinks();
  void create__Hlinks();
  void make_flip_loops();
  void change_operators();
};

//*************** CONSTRUCTOR ******************************************
LOOPS::LOOPS(int xsites, int ysites, int bondops, bool ob, long long rseed)
{
  irand.seed(rseed);
  drand.seed(rseed);
  
  dim1 = xsites;
  dim2 = ysites;
  OBC = ob;
  number_of_sites = dim1*dim2;
  number_of_bondops = 2*bondops;
  vlegs = 4*number_of_sites + 4*number_of_bondops;

  Vlinks.assign(vlegs, -99);
  Hlinks.assign(vlegs, -99);
  for(int i=0; i<vlegs; i++){ Hlinks[i]=i; }

  bops.resize(number_of_bondops,2);
  init_bonds.assign(number_of_sites, -99);  
}
//******** END of constructor *******************************************

/* nnbondlist()
   Create the list of nnbonds and a matrix called nn_mat for 2D,
   which gives the bond number for a pair of sites, and -99 if 
   they're not nearest neighbour sites.
*/
void LOOPS::nnbondlist()
{
  if(dim2==1){
    number_of_nnbonds = number_of_sites;
    if(OBC){number_of_nnbonds--;}//if we have open BCs there's one less bond
    
    nnbonds.resize (number_of_nnbonds,2);//make nnbonds the proper size
    nn_mat.resize(number_of_sites, number_of_sites);
    
    for(int i=0; i<number_of_nnbonds; i++){
      nnbonds(i,0) = i;
      nnbonds(i,1) = i+1;
      nn_mat(i,i+1) = i;
      nn_mat(i+1,i) = i;
    }
    if(!OBC){ //Add the periodic bond for this (1D) system
      nnbonds(number_of_nnbonds-1,0)=number_of_nnbonds-1;
      nnbonds(number_of_nnbonds-1,1)=0;
      nn_mat(number_of_nnbonds-1,0)=number_of_nnbonds-1;
      nn_mat(0,number_of_nnbonds-1)=number_of_nnbonds-1;
    }
    
  }
  else{
    number_of_nnbonds = 2*number_of_sites;
    if(OBC){number_of_nnbonds -= (dim1+dim2);}
    
    nnbonds.resize (number_of_nnbonds,2);

    //resize and initialize the matrix of nnbonds
    nn_mat.resize(number_of_sites, number_of_sites);
    for(int i=0; i<number_of_sites; i++){
      for(int j=0; j<number_of_sites; j++){
	nn_mat(i,j) = -99;
      }}

    int counter = 0;
    for(int i=0; i<dim1-1; i++){
      for(int j=0; j<dim2; j++){
	nnbonds(counter,0) = i+dim1*j;
	nnbonds(counter,1) = i+dim1*j+1;
	nn_mat(nnbonds(counter,0),nnbonds(counter,1))=counter;
	nn_mat(nnbonds(counter,1),nnbonds(counter,0))=counter;
	counter++;
      }
    }

    for(int i=0; i<dim1; i++){
      for(int j=0; j<dim2-1; j++){
	nnbonds(counter,0) = i+dim1*j;
	nnbonds(counter,1) = i+dim1*j + dim1;
	nn_mat(nnbonds(counter,0),nnbonds(counter,1))=counter;
	nn_mat(nnbonds(counter,1),nnbonds(counter,0))=counter;
	counter++;
      }
    }

    if(!OBC){
      for(int i=0; i<dim2; i++){
	nnbonds(counter,0) = i*dim1 + (dim1-1);
	nnbonds(counter,1) = i*dim1;
	nn_mat(nnbonds(counter,0),nnbonds(counter,1))=counter;
	nn_mat(nnbonds(counter,1),nnbonds(counter,0))=counter;
	counter++;
      }
      for(int j=0; j<dim1; j++){
	nnbonds(counter,0) = dim1*(dim2-1) + j;
	nnbonds(counter,1) = j;
	nn_mat(nnbonds(counter,0),nnbonds(counter,1))=counter;
	nn_mat(nnbonds(counter,1),nnbonds(counter,0))=counter;
	counter++;
      }
    }
    // make sure we have the proper number of nnbonds
    if(counter!=number_of_nnbonds){cout << "supererror" << endl;}
  }  

  /*  Check for nnbonds
      for(int i=0; i<number_of_nnbonds; i++){
      cout << nnbonds(i,0) << "," << nnbonds(i,1) << endl;
      }*/

  // resizing antiparallelness vectors
  init_antipar.assign(number_of_nnbonds, 0);
  init_whichbonds.resize(number_of_nnbonds);
  for(int i=0; i<number_of_nnbonds; i++){init_antipar[i]=i; init_whichbonds[i]=i;}
  antipar.resize(number_of_nnbonds);
  whichbonds.resize(number_of_nnbonds);
}
//******** END of creating NN bond list *********************************

/* Nnnbondlist() 
   Creates a list of the neighbouring bonds for a given bond.
   Works for OBC, PBC, 1D, and 2D.
   For 2D, uses the matrix nn_mat, which is created in nnbondlist().
 */
void LOOPS::Nnnbondlist()
{
  if(dim2==1){  // 1D case
    Nnnbonds.resize(number_of_nnbonds,2);

    // generate nearest nnbonds
    for(int i=0; i<number_of_nnbonds; i++){
      Nnnbonds(i,0)=(i+1)%number_of_nnbonds;
      Nnnbonds(i,1)=(i+number_of_nnbonds-1)%number_of_nnbonds;
    }
    if(OBC){  // For open BCs first and last bonds only have 1 neighbour
      Nnnbonds(0,1) = -99;
      Nnnbonds(number_of_nnbonds-1,0) = -99;
    }
  }
  
  else{  // the 2D case... more complicated...

    // Resize and initialize Nnnbonds
    Nnnbonds.resize(number_of_nnbonds,6);
    for(int i=0; i<number_of_nnbonds; i++){
      for(int j=0; j<6; j++){
	Nnnbonds(i,j)=-99;
      }
    }
    
    for(int i=0; i<number_of_nnbonds; i++){
      int counter=0;
      for(int j=1; j<number_of_sites; j++){
	
	//going through the nn matrix
	int bnum = nn_mat(nnbonds(i,0), (nnbonds(i,1) + j)%number_of_sites);
	//if the bond is a nn of the initial bond;
	if(bnum != -99){Nnnbonds(i,counter) = bnum; counter++;}
	

	//gdoing the same but moving down instead of across
	bnum = nn_mat((nnbonds(i,0) + j)%number_of_sites, nnbonds(i,1));
	if(bnum != -99){Nnnbonds(i,counter) = bnum; counter++;}
      }
    } 

    /*checking Nnnbonds
      for(int i=0; i<number_of_nnbonds; i++){
      for(int j=0; j<6; j++){
      cout << Nnnbonds(i,j) << ", ";
      }
      cout << endl;
      }
    */
  } 
}
//***** END of Nnnbondlist() **********************************************

void LOOPS::generate_ops()
{
  for(int i=0; i<number_of_bondops; i++)
    {
      bops(i,0) = irand() % number_of_nnbonds; //the bond being operated on
      bops(i,1) = 0; //0 = diagonal, 1 = off-diagonal
    }

  superbops.resize(number_of_bondops+number_of_sites,2);
  
    //initial state is dimerized..
    for(int i01=0; i01<number_of_sites; i01+=2){
      init_bonds[i01] = i01+1;
      init_bonds[i01+1] = i01;
      superbops(i01/2, 0) = nn_mat(i01,i01+1); 
      superbops(i01/2, 1) = 0;
      superbops(number_of_sites/2+number_of_bondops+i01/2,0) = nn_mat(i01,i01+1); 
      superbops(number_of_sites/2+number_of_bondops+i01/2,1) = 0;
    }

    for(int i=0; i<number_of_bondops; i++){
      superbops(number_of_sites/2+i,0)=bops(i,0);
      superbops(number_of_sites/2+i,1)=bops(i,1);
    }
}

void LOOPS::create_Vlinks()
{
  for(int i=0; i<vlegs; i+=2){
    Vlinks[i] = i+1;
    Vlinks[i+1] = i;
  }

  
}

void LOOPS::create__Hlinks()
{
  vector <int> last (number_of_sites,-99);
  for(int i=0; i<number_of_sites; i+=2){ 
    last[i]=i+2*(i/2+1); 
    last[i+1]=i+1+2*(i/2+1); 
  }
  
  int legnum = 0;
  // iterate through bond operators and create horizontal links
  // ******definitely fix the *order* problem if the probabilities change******
  for(int bopnum=number_of_sites/2; bopnum<number_of_bondops+number_of_sites; bopnum++){
    legnum = 4*bopnum;

    Hlinks[legnum] = last[nnbonds(superbops(bopnum,0),0)]; //does it matter if I screw
    Hlinks[last[nnbonds(superbops(bopnum,0),0)]] = legnum; //up the order??? because 
    last[nnbonds(superbops(bopnum,0),0)] = legnum + 2; // I am

    Hlinks[legnum+1] = last[nnbonds(superbops(bopnum,0),1)]; 
    Hlinks[last[nnbonds(superbops(bopnum,0),1)]] = legnum+1;
    last[nnbonds(superbops(bopnum,0),1)] = legnum + 3;
  }

  for(int i=0; i<number_of_bondops+number_of_sites; i++){
    cout << superbops(i,0) << endl;
  }
  cout << "------" << endl;
  /*
    for(int i=0; i<Hlinks.size(); i++){
    cout << i << ", " << Hlinks[i] << endl;
    }
    cout << "--------" << endl;
  */
}

void LOOPS::make_flip_loops()
{

  vector <int> crossloops(), loopnums(vlegs,-99);
  int loopnum(1), site(0), startsite(0), counter(0);
  bool which(0), flip=0;

  while(counter < vlegs-2){

    if(drand()<0.5){flip=1;}
    else{flip=0;}
    
    startsite = counter;
    site = Hlinks[counter];
    
    while((site == Hlinks[site])|(loopnums[counter]>0)){
      counter++;
      site = Hlinks[counter];
      startsite = counter;
    }
    loopnums[counter] = loopnum;
    loopnums[site] = loopnum;
    which = 0;

    if(counter > vlegs-2){break;}
    
    while(site!=startsite){

      loopnums[site] = loopnum;
      if(!which){
	site = Vlinks[site];
	if(flip){superbops(site/4,1) = (superbops(site/4,1)+1)%2;}
      }
      else{
	site = Hlinks[site];
      }
      which = !which;
      loopnums[site] = loopnum;
    }
    loopnum++;
    while(loopnums[counter]>0){counter++;}
  }

  for(int i=0; i<vlegs/4; i++){
    cout.width(4);
    cout.fill(' ');

    cout << left << i ;
    cout.width(4);
    cout << left << superbops(i,0); 
    cout.width(4);
    cout << left << superbops(i,1); 
    cout << endl;
  }
  /* for(int i=0; i<loopnums.size(); i++){
    cout << i << ", " << loopnums[i] << endl;}
  */
}   

void LOOPS::change_operators()
{
  antipar = init_antipar;
  whichbonds = init_whichbonds;
  int neighbs(0);
  if(dim2==1){neighbs=2};
  else{neighbs=6};

  for(int op=0; op<number_of_sites/2; op++){ //for the first N/2 operators (i.e. the edge)
    if(superbops(op,1)==1){                   //if operator is offdiagonal
      for(int i=0; i<neighbs; i++){           //change the antiparallelness of neighbouring bonds
	int loc = Nnnbonds(superbops(op,0),i);
	if(antipar(loc)>-1){                  //if bond is already antiparallel change to parallel
	  whichbonds.erase(antipar(loc));
	  antipar(loc) = -99;
	}
	else{
	  antipar(loc) = whichbonds.size();   //if bond is parallel change to antiparallel
	  whichbonds.push_back(whichbonds.size());
	}                                        
      }
    }                                          //otherwise (if diagonal) do nothing
  }

  // Now look at the *real* operators
  for(int op=number_of_sites/2; op<number_of_bondops+number_of_sites/2; op++){
    if(superbops(op,1)==1){                   //if operator is offdiagonal
      for(int i=0; i<neighbs; i++){           //change the antiparallelness of neighbouring bonds
	int loc = Nnnbonds(superbops(op,0),i);
	if(antipar(loc)>-1){                  //if bond is already antiparallel change to parallel
	  whichbonds.erase(antipar(loc));
	  antipar(loc) = -99;
	}
	else{
	  antipar(loc) = whichbonds.size();   //if bond is parallel change to antiparallel
	  whichbonds.push_back(whichbonds.size());
	}                                       
      }
    }
    else{          //if the operator is diagonal we need to change it randomly using whichbond..
      superbops(op,0) = whichbond(irand()%whichbonds.size());
    }
  }
}


#endif
