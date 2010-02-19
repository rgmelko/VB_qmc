//Jan 18, 2010 --- starting loop code
//Feb 19, 2010 --- trying to swap-ify

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

  int dim1, dim2,  number_of_sites; //the dimensions and number of sites
  int cross; /*the number of loops crossing the boundary i.e. the number of 
	       loops created by overlapping the propagated |VL> and |VR> */
  long long energyint; /*keeps track of the energy: 
			    energy = energyint*0.75/iterations */
  long long number_of_bondops, number_of_nnbonds;
  long long vlegs;/*the number of vertex legs including legs from the vertices
		    used to simulate the edge states |VL> and |VR> */
  long long iterations; //number of iterations per loop. Used for energy calc.
  bool OBC; //0 for PBC, 1 for OBC
  double energy; //the energy in non-integer form
  string bopfile; //the name of the file in which the bondops are stored
  
  vector <long long> Vlinks, Hlinks;//the vert and horizontal links for the LL
  vector <int> antipar, init_antipar; //keeps track of antiparallelness
  vector <int> isgood, init_isgood; //keeps track of 'good'ness
  vector <int> sides; /*store which side of the boundary a leg is on
			removing the need to check if it's higher or lower
			that the "middle" of the number of legs.... I'm not
			sure if this is more efficient or not */
  vector <int> VL, VR; //the propagated left and right VB states
  vector <int> whichloop; /* stores which loop number each site is in. used
			     in the energy measurement */

  iMatrix nnbonds;//list of all possible nnbonds.  Index is the bond number
  iMatrix nn_mat; /*matrix of the nnbonds. indices are sites
		    and contents are bond numbers.*/
  iMatrix Nnnbonds;//list of all nearest nnbonds
  iMatrix bops; // list of bond operators
  iMatrix superbops; //list of bond operators plus edges simulated via bops

  //CONSTRUCTOR
  LOOPS::LOOPS(int xsites, int ysites, int bondops, bool ob, long long its,
	       long long rseed, string bondopfile);

  void nnbondlist(); //creates list of nnbonds
  void Nnnbondlist(); //creates list of Nnnbonds
  void generate_ops(); //generates initial operators
  void create_Vlinks(); //creates vertical links
  void create__Hlinks(); //creates horizontal linkts (harder)
  void make_flip_loops(); //generates and flips loops (w/ prob 1/2)
  void take_measurement(); //measures energy at the moment
  void change__operators(); //changes the diagonal operators randomly

  void calculate_stuff(); //calculates energy at the moment
  void print_bops(); //prints the bond operators after every loop of iterations
  void read_bops(); //reads in bondops if theyre there, or runs generate_ops()

};

//*************** CONSTRUCTOR ******************************************
LOOPS::LOOPS(int xsites, int ysites, int bondops, bool ob, long long its, 
	     long long rseed, string bondopfile)
{
  irand.seed(rseed); //uses the random seed from the parameter file
  drand.seed(rseed); //enabling us to run fakely parallelize simulations
  
  dim1 = xsites; 
  dim2 = ysites;
  OBC = ob;
  number_of_sites = dim1*dim2; //calculates total number of sites
  //****changed**** multiplied by 2
  number_of_bondops = 2*2*bondops; /*the *real* number of bondops is multiplied
				   by 2, one set for |VL> and one for |VR> */
    //****changed**** multiplied first term by 2
  vlegs = 2*4*number_of_sites + 4*number_of_bondops; //number of vertex legs
  energyint = 0; energy = 0; //initialize the energy counters
  iterations = its; 
  bopfile = bondopfile; //name of the bond operator file

  Vlinks.assign(vlegs, -99); //set size and initialize
  Hlinks.assign(vlegs, -99); 
  //Initialize Hlinks so the "edge" sites are linked to themselves
  for(long long i=0; i<vlegs; i++){ Hlinks[i]=i; } 

  bops.resize(number_of_bondops,2); //set size of bops
  //create "sides"
  sides.assign(vlegs,0); 
  for(long long i=vlegs/2; i<vlegs; i++){sides[i]=1;}
  VL.assign(number_of_sites, -99); //set size of VL and VR
  VR=VL;
}

/********** nnbondlist() ***********************************************
   Uses:
     Global:
         dim2                   //
         number_of_nnbonds______sets value based on 1D/2D and PBC/OBC
         number_of_sites        //
         OBC                    //
         nnbonds[#nnbonds][2]___sized and filled
         nn_mat[#sites][#sites]_sized and filled
         dim1                   //
         init_antipar[#nnbonds]_sized and filled
         init_isgood[#nnbodds]__sized and filled
         antipar[#nnbonds]______sized
         isgood[#nnbonds]_______sized
      
      Local:
         counter //represents the number of the nnbond we're on

   Create the list of nnbonds and a matrix called nn_mat, which gives 
   the bond number for a pair of sites, and garbage if they're not nn 
   sites... maybe I should make it -99...
**********************************************************************/
void LOOPS::nnbondlist()
{
  if(dim2==1){
    number_of_nnbonds = number_of_sites;
    if(OBC){number_of_nnbonds--;}//if we have open BCs there's one less bond
    
    //****changed**** multiplied first dimension by 2
    nnbonds.resize (number_of_nnbonds*2,2);//make nnbonds the proper size
    //****changed**** multiplied dimensions by 2
    nn_mat.resize(number_of_sites*2, number_of_sites*2);
   
    for(int i=0; i<number_of_nnbonds; i++){
      nnbonds(i,0) = i;
      nnbonds(i,1) = i+1;
      nn_mat(i,(i+1)%number_of_sites) = i;
      nn_mat((i+1)%number_of_sites,i) = i;
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
      }
    }

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
  
  //****changed**** added this part in to double nnbonds and nnmat
      for(int b=number_of_nnbonds; b<2*number_of_nnbonds; b++){
	int a = number_of_nnbonds;
	nnbonds(b,0) = nnbonds(b-a,0)+number_of_sites;
	nnbonds(b,1) = nnbonds(b-a,1)+number_of_sites;
	nn_mat(nnbonds(b,0),nnbonds(b,1))=b;
	nn_mat(nnbonds(b,1),nnbonds(b,0))=b;
      }
  /**** end of this change *****/
	
  //****changed**** now multiplying #nnbonds by 2
      number_of_nnbonds *=2; 
  
  // resizing antiparallelness vectors
  init_antipar.assign(number_of_nnbonds, 0);
   
  init_isgood.resize(number_of_nnbonds);
  for(int i=0; i<number_of_nnbonds; i++){init_antipar[i]=1; init_isgood[i]=i;}
  antipar.resize(number_of_nnbonds);
  isgood.resize(number_of_nnbonds);

  //****changed**** changing #nnbonds back now
      number_of_nnbonds /=2;
}

/************* Nnnbondlist() *************************************************
 Uses:
  Global:
   dim2                   //
   Nnnbonds[#Nnnbonds][2] //sizes and fills
   number_of_nnbonds      //
   nn_mat[#sites][#sites] //
   OBC                    //
   number_of_sites        //

  Local:
   counter__goes from 0-3 (i think) and represents the number of neighbours
            we've found so far.
   bnum_____the number of the bond we're looking at

   Creates a list of the neighbouring bonds for a given bond.
   Works for OBC, PBC, 1D, and 2D.
   For 2D, uses the matrix nn_mat, which is created in nnbondlist().
*****************************************************************************/
void LOOPS::Nnnbondlist()
{
  int num_neighbs = 0;
  if(dim2==1){  // 1D case
    num_neighbs = 2;
    //****changed**** multiplied first dimension by 2
    Nnnbonds.resize(2*number_of_nnbonds,num_neighbs);

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
    num_neighbs = 6;
    // Resize and initialize Nnnbonds
    //****changed**** multiplied first dimension by 2
    Nnnbonds.resize(2*number_of_nnbonds,num_neighbs);
    for(int i=0; i<number_of_nnbonds; i++){
      for(int j=0; j<num_neighbs; j++){
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
  } 
  //****changed**** added this in to copy Nnnbonds.. to double it.. sorta
  for(int q=number_of_nnbonds; q<2*number_of_nnbonds; q++){
    for(int s=0; s<num_neighbs; s++){
      int t=number_of_nnbonds;
      Nnnbonds(q,s)=Nnnbonds(q-t,s)+number_of_nnbonds;
    }
  }


  //****changed**** for realisies doubling #nnbonds and #sites
				    number_of_nnbonds *= 2;
  number_of_sites *= 2;
}

/***** generate_ops() *****************************************************
 Uses:
   number_of_bondops      //
   bops[#bondops][2]______fills randomly
   number_of_nnbonds      //
   superbops[#bondops][2]_sizes and fills with edges and ops from bops
   number_of_sites        //
   nn_mat[#sites][#sites] //

***************************************************************************/
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
      superbops(i01/2, 0) = nn_mat(i01,i01+1); 
      superbops(i01/2, 1) = 0;
      superbops(number_of_sites/2+number_of_bondops+i01/2,0)=nn_mat(i01,i01+1);
      superbops(number_of_sites/2+number_of_bondops+i01/2,1) = 0;
    }

    for(int i=0; i<number_of_bondops; i++){
      superbops(number_of_sites/2+i,0)=bops(i,0);
      superbops(number_of_sites/2+i,1)=bops(i,1);
    }
}
/************ create_Vlinks() ************************************************
 Uses:
   vlegs         //
   Vlinks[vlegs] //filled with the vertical links

******************************************************************************/
void LOOPS::create_Vlinks()
{
  for(long long i=0; i<vlegs; i+=2){
    Vlinks[i] = i+1;
    Vlinks[i+1] = i;
  }

  
}
/************ create__Hlinks() ************************************************
 Uses:
  Global:
   number_of_sites        //
   number_of_bondops      //
   Hlinks[vlegs]__________filled with horizontal links
   nnbonds[#nnbonds][2]   //
   superbops[#bondops][2] //

  Local:
   last[#sites]_stores the last vertex leg corresponding to a site
   legnum_______the vertex leg number of the current bondop
   bopnum_______counts through the bondops
*******************************************************************************/
void LOOPS::create__Hlinks()
{
  vector <long long> last (number_of_sites,-99);
  for(int i=0; i<number_of_sites; i+=2){ 
    last[i]=i+2*(i/2+1); 
    last[i+1]=i+1+2*(i/2+1); 
  }
  
  long long legnum = 0;
  // iterate through bond operators and create horizontal links
  // ******definitely fix the *order* problem if the probabilities change******
  long long bopnum = number_of_sites/2;
  for(bopnum; bopnum<number_of_bondops+number_of_sites; bopnum++){
    legnum = 4*bopnum;

    Hlinks[legnum] = last[nnbonds(superbops(bopnum,0),0)];//does it matter if
    Hlinks[last[nnbonds(superbops(bopnum,0),0)]] = legnum;//I screw up the order
    last[nnbonds(superbops(bopnum,0),0)] = legnum + 2; //? because I am

    Hlinks[legnum+1] = last[nnbonds(superbops(bopnum,0),1)]; 
    Hlinks[last[nnbonds(superbops(bopnum,0),1)]] = legnum+1;
    last[nnbonds(superbops(bopnum,0),1)] = legnum + 3;
  }
}

/************ make_flip_loops() **********************************************
 Creates loops and flips them with probability 1/2

 Global:
   VL[#sites]___________get filled by looking at links crossing the boundary
   VR[#sites]___________same
   whichloop[#sites]____get filled. stores the loop number for each site
   cross________________counts the number of loops crossing the boundary
   vlegs                //
   Hlinks[vlegs]        //
   sides[vlegs]         //
   nnbonds[#nnbonds][2] //
   superbops[#bops][2]  //
   Vlinks[vlegs]        //

 Local:
   loopnums[vlegs].stores the loop number for each vertex leg
   loopnum_________the loop number we're currently looking at
   startsite_______the start site for the current loop
   counter_________starts at the beginning of the vertex legs, goes to the end
   site____________the current vertex leg we're looking at
   which___________0 for looking at vertical links, 1 for horizontal
   flip____________0 if we're not flipping this loop, 1 if we are
   firstcross______the first site involved in a loop crossing the boundary. 
                   Used to get the last bond at the boundary.
   lastcross_______the last site involved in a loop crossing the boundary
   current_________the current site we're looking at (corresponds to some 
                   vertex leg.
   right___________0 if the current loop is going left over the boundary, 1
                   if it's going right (the operators go from left to right)
   boolcross_______1 if the loop that was just completed crosses the boundary
                   so we can increase the cross counter. O otherwise.
   
*****************************************************************************/
void LOOPS::make_flip_loops()
{
  vector <int> loopnums(vlegs,-99);
  int loopnum(1), startsite(0); 
  long long counter(0), site(0);
  bool which(0), flip=0;
  int firstcross(-99),lastcross(-99), current(0); 
  int right=-99;
  bool boolcross=0;

  VL.assign(number_of_sites, -99);
  VR = VL;
  whichloop = VL;
  cross = 0;

  while(counter < vlegs-2){

    if(drand()<0.5){flip=1;}
    else{flip=0;}

    startsite = counter;        // set the initial site (startsite)
    site = Hlinks[counter];     // the site connected to startsite horizontally

    //making sure site isn't bonded to itself (edges are bonded to themselves)
    //and checking that the startsite isn't already in a loop
    while(((site == Hlinks[site])|(loopnums[counter]>0))&&(counter<vlegs-2)){ 
      counter++;  //changing startsite
      site = Hlinks[counter]; //site connected to new startsite
      startsite = counter;    //setting new startsite
    }
    //breaks if we get to the last or 2nd last site
    if(counter > vlegs-2){break;} 
    loopnums[counter] = loopnum; //including startsite in new loop
    loopnums[site] = loopnum;    //adding the next site to the loop
    which = 0;      //"which" checks if we're on horiz(1) or vert(0)     

    if(sides[startsite]!=sides[site]){
   
      firstcross =  nnbonds(superbops(site/4,0),site%2);
      lastcross = firstcross;
      if(sides[startsite]<sides[site]){right=1;}
      else{right=0;}
      boolcross++;
      whichloop[firstcross]=loopnum;
    }
    
    //while loops ends when we get back to the startsite
    while(site!=startsite){
      
      //VERTICAL LINKS
      if(!which){                
	site = Vlinks[site];
	if(flip){superbops(site/4,1) = (superbops(site/4,1)+1)%2;}
      }
      //HORIZONTAL LINKS
      else{
    
	if(sides[site]!=sides[Hlinks[site]]){  //if it crosses the boundary
	  current = nnbonds(superbops(site/4,0),site%2);
	  boolcross++;
	  //if this is the first crossing for this loop set firstcross
	  if(lastcross<0){//and set right to show if it's crossing left or right
	    firstcross =  current;
	    if (sides[site]<sides[Hlinks[site]]){right=1;}
	    else {right=0;}
	  }
	  //if it's not the first crossing
	  else{
	    if(right){
	      VR[current] = lastcross;
	      VR[lastcross] = current;
	    }
	    else{
	      VL[current] = lastcross;
	      VL[lastcross] = current;
	    }
	    right = (right+1%2);  //change cross direction
	  }
	  lastcross = current;
	  whichloop[current]=loopnum;
	}
	
	site = Hlinks[site];
      }
      which = !which; //changes from horiz(1) to vert(0) or vice versa
      loopnums[site] = loopnum; //adds next site to loop
    }
    if(lastcross>-1){
      if(right){ VR[firstcross]=lastcross; VR[lastcross]=firstcross;}
      else{ VL[firstcross]=lastcross; VL[lastcross]=firstcross;}
    }
    if(boolcross){cross++;}
    boolcross = 0;
    loopnum++; //loop is finished, go to next loop
    lastcross = -99; right = -99;
    //if counter is already in a loop increase counter
    while(loopnums[counter]>0){counter++; if(counter>vlegs){break;}}
  }
}   
/************ take_measurement() *********************************************
 Global:
   number_of_nnbonds
   nnbonds[#nnbonds][2]
   whichloop[#sites]
   energyint

 Local:
   mdiff
   a
   b

*****************************************************************************/
void LOOPS::take_measurement()
{
  //using VL, VR, cross....
  int mdiff=0;

  for(int i=0; i<number_of_nnbonds; i++){

    int a(0),b(0);
    a = nnbonds(i,0);
    b = nnbonds(i,1);
    
    if(whichloop[a]!=whichloop[b]){mdiff++;}
  }
  energyint += mdiff - number_of_nnbonds;
}
/************ change__operators() ********************************************
*****************************************************************************/
void LOOPS::change__operators()
{
  antipar = init_antipar;
  isgood = init_isgood;
  int neighbs(0);
  if(dim2==1){neighbs=2;}
  else{neighbs=6;}

  //for the first N/2 operators (i.e. the edge)
  for(int op=0; op<number_of_sites/2; op++){ 
    if(superbops(op,1)==1){                   //if operator is offdiagonal

      for(int i=0;i<neighbs;i++){//change antiparallelness of neighboring bonds

	int loc = Nnnbonds(superbops(op,0),i);
       	if(loc<0){continue;}//goto start of loop if its OBC&that nn doesnt exist
	if(antipar[loc]==1){//if bond is already antiparallel change to parallel

	  int i=0;
	  do{ i++; }
	  while(isgood[i]!=loc);
	  isgood.erase(isgood.begin() + i);
	  antipar[loc]--;
	}
	else{
	  antipar[loc]++;   //if bond is parallel change to antiparallel
	  isgood.push_back(loc);
	}                                        
      }
    }                       //otherwise (if diagonal) do nothing
  }
  // Now look at the *real* operators
  long long op = number_of_sites/2;
  for(op; op<number_of_bondops+number_of_sites/2; op++){
    if(superbops(op,1)==1){                         //if operator is offdiagonal
      for(int i=0;i<neighbs; i++){//change antiparallelness of neighboring bonds

	int loc = Nnnbonds(superbops(op,0),i);

	if(antipar[loc]==1){//if bond's already antiparallel change to parallel

	  int j=-1;
	  do{ j++; }
	  while(isgood[j]!=loc);
	  isgood.erase(isgood.begin()+j);
	  antipar[loc]--;
	}
	else{
	  antipar[loc]++;   //if bond is parallel change to antiparallel
	  isgood.push_back(loc);
	}                                       
      }
    }      //if the operator is diagonal we need to change it randomly
    else{       //using whichbond..              
      superbops(op,0) = isgood[irand()%isgood.size()];
    }
  }
}

/************ calculate_stuff() ***********************************************
******************************************************************************/
void LOOPS::calculate_stuff()
{
  energy = 0.5*(energyint*0.75)/iterations;
  energyint = 0;
}

/************ print_bops() ****************************************************
******************************************************************************/
void LOOPS::print_bops()
{
  ofstream bout(bopfile.c_str());
  for(int i=0; i<number_of_bondops+number_of_sites; i++){
    bout << superbops(i,0) << endl << superbops(i,1) << endl;
  }
}

/************ read_bops() *****************************************************
| If the bond operator file (filename given in parameters file)
******************************************************************************/
void LOOPS::read_bops()
{
  ifstream bin(bopfile.c_str());
  
  if(bin.fail()){ generate_ops(); }
  
  else{ 
    superbops.resize(number_of_bondops+number_of_sites,2);
    for(int i=0; i<number_of_bondops+number_of_sites; i++){
      bin >> superbops(i,0) >> superbops(i,1);
    }
  }
}
#endif
