//March 21, 2012 --- Trying to modify this code to do the J-Q model
//March 15, 2012 --- Modifying latest cylinder code so operator choice
//                   works differently.  Instead of keeping a list of
//                   "good" operators, follow the steps on today's date 
//                   in my "maze" notebook:
// 1. Choose Operator type w/ prob: P_p = Q N_p / (Q N_P + J N_b) and 
//    P_b = 1 - P_p
// 2. Randomly choose location for given type of operator
// 3. If operator type is allowed by spins, move on. If not go back to 1.
//Feb 2012 --- try to adapt for cylinder geometry
//Nov 2010 --- 2D RATIO LOOP HEADER
//Feb 19, 2010 --- trying to swap-ify
//Jan 18, 2010 --- starting loop code


#ifndef jq_header
#define jq_header

#include "header.h"
#include "matrix.h"
#include <vector>

class LOOPS
{
 public:

  MTRand drand; //drand() gives you a random double precision number
  MTRand_int32 irand; // irand() gives you a random integer

  int flip; // -1 for bare swap, 0 for ratio 1, etc

  double J,Q;
  int Lx, Ly,  number_of_sites; //the dimensions and number of sites
  //Note: things are periodic in the Ly direction, open in Lx

  int cross; /*the number of loops crossing the boundary i.e. the number of 
	       loops created by overlapping the propagated |VL> and |VR> */
  long long energyint; /*keeps track of the energy: 
			    energy = energyint*0.75/iterations */
  long long number_of_bondops, number_of_nnbonds, numPlaquettes;
  int num2site;
  double middle, vlegMiddle;
  double plaqProb;
  long long vlegs;/*the number of vertex legs including legs from the vertices
		    used to simulate the edge states |VL> and |VR> */
  long long iterations; //number of iterations per loop. Used for energy calc.
  double energy; //the energy in non-integer form
  string bopfile; //the name of the file in which the bondops are stored
  
  vector <long long> Vlinks, Hlinks;//the vert and horizontal links for the LL
  vector <double> vleg2op; //list of which leg/4 corresponds to what operator
  vector <int> spins, init_spins; //keeps track of spins for swaperation

  vector <int> VL, VR; //the propagated left and right VB states
  vector <int> whichloop; /* stores which loop number each site is in. used
			     in the energy measurement */
  vector <double> entropy, entropy_final;

  iMatrix nnbonds;//list of all possible nnbonds.  Index is the bond number
  iMatrix plaquettes;/*list of all possible plaquette operators.  Index
		       is plaquette operator number. 2 per plaquette */
  iMatrix nn_mat; /*matrix of the nnbonds. indices are sites
		    and contents are bond numbers.*/
  iMatrix bops; // list of bond operators
  iMatrix superbops; //list of bond operators plus edges simulated via bops
                     //superbops(:,0) is the bond the operator acts on 
                     //         (:,1) is 0 for Diag, 1 for OffDiag
                     //         (:,2) is 1 for Bond, 2 for Plaquette
                     //               and 3 for Edge state.

  //CONSTRUCTOR
  LOOPS(double jay, double que, int xsites, int ysites, int flips, 
	int bondops, long long its, long long rseed, string bondopfile);

  void operatorLists(); //creates list of nnbonds & plaquettes
  void generate_ops(); //generates initial operators
  void create__Hlinks(); //creates horizontal links (harder)
  void make_flip_loops(); //generates and flips loops (w/ prob 1/2)
  void take_measurement(); //measures energy at the moment
  void change__operators(); //changes the diagonal operators randomly

  void swaperator(); //the swap operator

  void calculate_stuff(); //calculates energy at the moment
  void print_bops(); //prints the bond operators after every loop of iterations
  void read_bops(); //reads in bondops if theyre there, or runs generate_ops()

};

//*************** CONSTRUCTOR ******************************************
LOOPS::LOOPS(double jay, double que, int xsites, int ysites, int flips, 
	     int bondops, long long its, long long rseed, string bondopfile)
{
  irand.seed(rseed); //uses the random seed from the parameter file
  drand.seed(rseed); //enabling us to run fakely parallelize simulations
  
  J = jay;
  Q = que;

  Lx = xsites; 
  Ly = ysites;
  flip = flips;
  number_of_sites = Lx*Ly; //calculates total number of sites
  //****changed**** multiplied by 2
  number_of_bondops = 2*2*bondops; /*the *real* number of bondops is multiplied
				   by 2, one set for |VL> and one for |VR> */

  //initialize the energy counters
  energyint = 0; energy = 0; 
  iterations = its; 
  //name of the bond operator file
  bopfile = bondopfile; 

  entropy.assign(xsites,0);
  entropy_final = entropy;

  int maxVlegs =  2*4*number_of_sites + 8*number_of_bondops;
  //Initialize Vlinks
  Vlinks.assign(maxVlegs, -99); //set size and initialize
  for(long long i=0; i<maxVlegs; i+=2){
    Vlinks[i] = i+1;
    Vlinks[i+1] = i;
  }

  bops.resize(number_of_bondops,3); //set size of bops

  VL.assign(number_of_sites*2, -99); //set size of VL and VR
  VR=VL;

  //set initial spin values
  init_spins.assign(number_of_sites*2,0); 
  int k;
  //for one replica
  for(int i=0; i<Ly; i++){
    for(int j=0; j<Lx; j++){
      k = i+Ly*j;
      if((i+j)%2==0){init_spins[k]=1;}
    }
  }
  //fill the other replica
  for(int i=0; i<number_of_sites; i++){
    init_spins[i+number_of_sites]=init_spins[i];
  }
  spins = init_spins; 

}

/********** operatorLists() *******************************************
   Uses:
     Global:
         Ly                     //
         number_of_nnbonds______sets value based on 1D/2D and PBC/OBC
         number_of_sites        //
	 numPlaquettes
         nnbonds[#nnbonds][2]___sized and filled
         nn_mat[#sites][#sites]_sized and filled
	 plaquettes[#plaquettes][4]___sized and filled
         Lx                     //
      
      Local:
         counter //represents the number of the nnbond we're on

   Create the list of nnbonds and a matrix called nn_mat, which gives 
   the bond number for a pair of sites, and garbage if they're not nn 
   sites... maybe I should make it -99...
**********************************************************************/
void LOOPS::operatorLists()
{
  // should put in a special case for Ly==2 ??

  //Case of 1D open chain
  //This case doesn't exist for J-Q.. maybe delete it.
  if(Ly==1){
    number_of_nnbonds = (Lx-1);
    //changed multiplied first dimension by 2
    nnbonds.resize (number_of_nnbonds*2,2);//make nnbonds the proper size
    for(int i=0;i<number_of_nnbonds*2;i++){nnbonds(i,0)=-99;nnbonds(i,1)=-99;}
    //changed multiplied dimensions by 2
    nn_mat.resize(number_of_sites*2, number_of_sites*2);
    for(int i=0;i<number_of_sites*2;i++){
      for(int j=0;j<number_of_sites*2;j++){
	nn_mat(i,j)=-99;
      }
    }
    for(int i=0; i<number_of_nnbonds; i++){
      nnbonds(i,0) = i;
      nnbonds(i,1) = i+1;
      nn_mat(i,(i+1)%number_of_sites) = i;
      nn_mat((i+1)%number_of_sites,i) = i;
    }  
  }

  //Generic 2D case and 1D ring.
  else{
    number_of_nnbonds = (Lx-1)*Ly + Lx*Ly; // = Ly*(2*Lx-1)
    numPlaquettes = 2*(Lx-1)*Ly; // = Ly*(2*Lx-2)
    
    //changed multiplied first dimension by 2
    nnbonds.resize (2*number_of_nnbonds,2);
    plaquettes.resize (2*numPlaquettes,4);
    for(int i=0;i<number_of_nnbonds*2;i++){
      nnbonds(i,0)=-99;nnbonds(i,1)=-99;
    }
    for(int i=0;i<numPlaquettes*2;i++){
      plaquettes(i,0)=-99;
      plaquettes(i,1)=-99;
      plaquettes(i,2)=-99;
      plaquettes(i,3)=-99;
    }

    //resize and initialize the matrix of nnbonds
    //changed multiplied dimensions by 2
    nn_mat.resize(2*number_of_sites, 2*number_of_sites);
    for(int i=0; i<2*number_of_sites; i++){
      for(int j=0; j<2*number_of_sites; j++){
	nn_mat(i,j) = -99;
      }
    }

    //fill nnbonds and nn_mat
    //ybonds!! (vertical, periodic)
    int bCounter = 0;
    int pCounter = 0;
    for(int x=0; x<Lx; x++){
      for(int y=0; y<Ly; y++){
	nnbonds(bCounter,0) = Ly*x+y;
	nnbonds(bCounter,1) = Ly*x+(y+1)%Ly;
	nn_mat(nnbonds(bCounter,0),nnbonds(bCounter,1))=bCounter;
	nn_mat(nnbonds(bCounter,1),nnbonds(bCounter,0))=bCounter;
	bCounter++;
	if((x+1)<Lx){
	  //first two sites are the same as the bond
	  plaquettes(pCounter,0) = Ly*x+y;
	  plaquettes(pCounter,1) = Ly*x+(y+1)%Ly;
	  //next two sites are shifted right (in the Lx direction)
	  plaquettes(pCounter,2) = Ly*(x+1)+y;
	  plaquettes(pCounter,3) = Ly*(x+1)+(y+1)%Ly;
	  pCounter++;
	}
      }
    }
    
    //xbonds... (horizontal, not periodic)
    for(int i=0; i<(Lx-1)*Ly; i++){
      nnbonds(bCounter,0) = i;
      nnbonds(bCounter,1) = i+Ly;
      nn_mat(nnbonds(bCounter,0),nnbonds(bCounter,1))=bCounter;
      nn_mat(nnbonds(bCounter,1),nnbonds(bCounter,0))=bCounter;
      bCounter++;
      plaquettes(pCounter,0)=i;
      plaquettes(pCounter,1)=i+Ly;
      plaquettes(pCounter,2)=(i+1)%Ly;
      plaquettes(pCounter,3)=(i+1)%Ly+Ly;
      pCounter++;
    }

    // make sure we have the proper number of nnbonds
    if(bCounter!=number_of_nnbonds){cout << "supererror" << endl;}
    if(pCounter!=numPlaquettes){cout << "super Plaquette counting error" << endl;}
    
  }
  
  int a = number_of_nnbonds;
  //changed added this part in to double nnbonds and nnmat
  for(int b=number_of_nnbonds; b<2*number_of_nnbonds; b++){
    nnbonds(b,0) = nnbonds(b-a,0)+number_of_sites;
    nnbonds(b,1) = nnbonds(b-a,1)+number_of_sites;
    nn_mat(nnbonds(b,0),nnbonds(b,1))=b;
    nn_mat(nnbonds(b,1),nnbonds(b,0))=b;
  }
  a = numPlaquettes;
  for(int b=numPlaquettes;b<2*numPlaquettes;b++){
    plaquettes(b,0) = plaquettes(b-a,0)+number_of_sites;
    plaquettes(b,1) = plaquettes(b-a,1)+number_of_sites;
    plaquettes(b,2) = plaquettes(b-a,2)+number_of_sites;
    plaquettes(b,3) = plaquettes(b-a,3)+number_of_sites;
  }
  // end of this change //
  
  //print out nn_mat !!!!! 
  //    for(int i=0;i<nn_mat.length();i++){
  //      for(int j=0; j<nn_mat.width();j++){
  //	if(nn_mat(i,j)==-99){cout<<". ";}
  //	else{cout << nn_mat(i,j) <<" " ;}
  //     }
  //      cout << endl;
  //   }
  
  number_of_nnbonds *= 2;
  number_of_sites *= 2;
  numPlaquettes *=2;

  plaqProb = Q*numPlaquettes/(J*number_of_nnbonds + Q*numPlaquettes);
  
}


/***** generate_ops() *****************************************************
 Uses:
   number_of_bondops      //
   bops[#bondops][3]______fills randomly
   number_of_nnbonds      //
   superbops[#bondops][3]_sizes and fills with edges and ops from bops
   number_of_sites        //
   nn_mat[#sites][#sites] //

***************************************************************************/
void LOOPS::generate_ops()
{
    
  //picks initial operators acting on antiparallel spins
  int temp(-99);
  
  //num2site is the effective number of 2 site operators (counting a plaquette
  // operator as 1 bond operators)
  num2site=0;
  //middle gives the middle bond operator if you equate a 4site operator 
  //to 2 bondops
  middle=0;
  
  //Just put this here in case I want to experiment with different Qs 
  //without changing the param file.
  //Can delete it later.  The real one is at the end of operatorLists()
  plaqProb = Q*numPlaquettes/(J*number_of_nnbonds + Q*numPlaquettes);
  
  for(int i=0; i<number_of_bondops;){
    if(drand()<plaqProb){
      //choose a plaquette!!
      temp = irand()%numPlaquettes;
      //if the spins are antiparallel, assign it
      if(spins[plaquettes(temp,0)]+spins[plaquettes(temp,1)]==1 && 
	 spins[plaquettes(temp,2)]+spins[plaquettes(temp,3)]==1){
	bops(i,0) = temp;
	bops(i,1) = 0;
	bops(i,2) = 2; //2 for plaquette, right?
	i++;
	num2site+=2;
      }
    }
    else{
      //choose a bondop
      temp = irand()%number_of_nnbonds;
      //if the spins are antiparallel, assign it
      if(spins[nnbonds(temp,0)]+spins[nnbonds(temp,1)]==1){
	bops(i,0) = temp;
	bops(i,1) = 0;
	bops(i,2) = 1; //1 for bond operator
	i++;
	num2site++;
      }
    }
    if(i==number_of_bondops/2){middle = num2site-0.5;}
  }

  vlegs = 4*number_of_sites + 4*num2site;
  vlegMiddle = 2*number_of_sites + 4*(middle+0.5) - 0.5;

  superbops.resize(number_of_bondops+number_of_sites,3);
  //  for(int i=0;i<number_of_sites+number_of_bondops;i++){
  //    superbops(i,0)=-99;
  //    superbops(i,1)=-99;
  //    }

  /*==============================================================
    Setting the initial dimerized state in terms of bond operators
    ==============================================================*/
  //initial state is dimerized, but either Lx or Ly can be odd.
  //Dimerize in x direction if Lx is even
  //Otherise dimerize in Ly direction
  int bondd=0;
  if(Lx%2==0){
    for(int ix=0; ix<Lx; ix+=2){
      for(int iy=0; iy<Ly; iy++){
	superbops(bondd,0) = nn_mat(iy+ix*Ly,iy+(ix+1)*Ly);
	superbops(bondd,1) = 0;
	superbops(bondd,2) = 1;
	superbops(bondd+Lx*Ly/2,0) = nn_mat(Lx*Ly+iy+ix*Ly,Lx*Ly+iy+(ix+1)*Ly);
	superbops(bondd+Lx*Ly/2,1) = 0;
	superbops(bondd+Lx*Ly/2,2) = 1; //Edges are bond operators
	//copy into the right side initial state
	bondd++;
      }
    }
  }
  else{//Lx is odd, so Ly *must* be even
    if(Ly%2==1){ cout<<"Error, both dimensions are odd!!"<<endl; exit(1);}

    for(int ix=0; ix<Lx; ix++){
      for(int iy=0; iy<Ly; iy+=2){
	//cout << iy+ix*Ly<< ", "<<iy+1+ix*Ly<<"   
	// "<<Lx*Ly+iy+ix*Ly<<", "<<Lx*Ly+iy+1+ix*Ly<<endl;
	superbops(bondd,0) = nn_mat(iy+ix*Ly,iy+1+ix*Ly);
	superbops(bondd,1) = 0;
	superbops(bondd,2) = 1;
	superbops(bondd+Lx*Ly/2,0) = nn_mat(Lx*Ly+iy+ix*Ly,Lx*Ly+iy+1+ix*Ly);
	superbops(bondd+Lx*Ly/2,1) = 0;
	superbops(bondd+Lx*Ly/2,2) = 1; //These are bond operators
	//copy into the right side initial state
	bondd++;
      }
    }
  }
  //check to make sure there's the right number of bonds in the initial state
  if(bondd!=Lx*Ly/2){
    cout<<"ERROR: Number of bonds didn't work. Check generate_ops()\n"; 
    cout<<"bondd = "<< bondd << endl;
    exit(1);
  }

  //Copy the starting states for all 4 systems
  for(int iall=0; iall<Lx*Ly; iall++){
    superbops(Lx*Ly + number_of_bondops + iall,0) = superbops(iall,0);
    superbops(Lx*Ly + number_of_bondops + iall,1) = 0;
    superbops(Lx*Ly + number_of_bondops + iall,2) = superbops(iall,2);
  }

  for(int i=0; i<number_of_bondops; i++){
    superbops(number_of_sites/2+i,0)=bops(i,0);
    superbops(number_of_sites/2+i,1)=bops(i,1);
    superbops(number_of_sites/2+i,2)=bops(i,2);
  }

  //  for(int i=number_of_bondops;i<number_of_bondops+number_of_sites; i++){
  //    cout << "superbops("<<i<<",0) = "<<superbops(i,0)<<endl;
  //  }
  
}

/************ create__Hlinks() ************************************************
 Uses:
  Global:
   number_of_sites        //
   number_of_bondops      //
   Hlinks[vlegs]__________filled with horizontal links
   nnbonds[#nnbonds][2]   //
   superbops[#bondops][3] //

  Local:
   last[#sites]_stores the last vertex leg corresponding to a site
   legnum_______the vertex leg number of the current bondop
   bopnum_______counts through the bondops


i assume vertex legs are labelled like so (Feb 2012):

                0 ----------- 2
                       |
                       |
                       |
                1 ----------- 3

*******************************************************************************/
void LOOPS::create__Hlinks()
{ 
  //The last vertex leg number for a given site
  vector <long long> last (number_of_sites,-99);

  //put initial VB state as initial vertices in linked list.
  int vertex=0; 
  for(int i=0; i<number_of_sites/2; i++){   
    last[nnbonds(superbops(i,0),0)]= 4*vertex+2;
    last[nnbonds(superbops(i,0),1)]= 4*vertex+3;
    vertex++;
  }
  
  //check that everything got filled 
  //  for(int i=0;i<last.size();i++){
  //    if(last[i]<0){cout <<"omg!!! "<< i <<endl;}// exit(1);}
  //  }

  long long legnum = 0;
  // iterate through bond operators and create horizontal links

  // Size Hlinks and vleg2op
  Hlinks.assign(vlegs, -99); 
  vleg2op.assign(vlegs, -99);
  //Initialize Hlinks so the "edge" sites are linked to themselves
  for(int i=0; i<number_of_sites/2; i++){
    for(int j=0; j<4; j++){Hlinks[i*4+j]=i*4+j;}
    vleg2op[i] = i;
  }
  for(int i=vlegs/4-number_of_sites/2; i<vlegs/4; i++){
    for(int j=0; j<4; j++){Hlinks[i*4+j]=i*4+j;}
    vleg2op[i]= i-vlegs/4+number_of_sites+number_of_bondops;
  }

  // the first N/2 bondops are the initial VB config. Pick up from there
  long long bopnum = number_of_sites/2;

  // create Hlinks for the operators after the init VB config.. 
  // but also include the bondops for the other initial VB config 
  // (on the right side) |VR> 

  long long templegnum;
  //initialize the leg we start on
  legnum = 4*bopnum;

  //for the LHS of the bondops
  for(bopnum; legnum<vlegMiddle; bopnum++){

    //figure out operator type
    //bond operator
    if(superbops(bopnum,2)==1){
      //join the first leg to the... 
      templegnum = nnbonds(superbops(bopnum,0),0);
      Hlinks[legnum] = last[templegnum];//does it matter if
      Hlinks[last[templegnum]] = legnum;//I screw up the order
      last[templegnum] = legnum + 2; //? because I am
      
      templegnum = nnbonds(superbops(bopnum,0),1);
      Hlinks[legnum+1] = last[templegnum]; 
      Hlinks[last[templegnum]] = legnum+1;
      last[templegnum] = legnum + 3;
      
      vleg2op[legnum/4] = bopnum;
      legnum+=4;
    }
    else{ //it's a plaquette operator!!!
      templegnum = plaquettes(superbops(bopnum,0),0);
      Hlinks[legnum] = last[templegnum];
      Hlinks[last[templegnum]] = legnum;
      last[templegnum] = legnum + 2;

      templegnum = plaquettes(superbops(bopnum,0),1);
      Hlinks[legnum+1] = last[templegnum]; 
      Hlinks[last[templegnum]] = legnum+1;
      last[templegnum] = legnum + 3;
      
      vleg2op[legnum/4] = bopnum;
      legnum +=4;
      
      //now the other side of the plaquette
      templegnum = plaquettes(superbops(bopnum,0),2);
      Hlinks[legnum] = last[templegnum];
      Hlinks[last[templegnum]] = legnum;
      last[templegnum] = legnum + 2;

      templegnum = plaquettes(superbops(bopnum,0),3);
      Hlinks[legnum+1] = last[templegnum]; 
      Hlinks[last[templegnum]] = legnum+1;
      last[templegnum] = legnum + 3;

      vleg2op[legnum/4] = bopnum+0.5;
      legnum +=4;  
    }
  }
  
  /*************************************************************
    NOW stop to switch the connections within region A
   *************************************************************/
  long long a,b,c;
  
  //flip 1 column of Ly at a time
  for(int iz=0; iz<=flip; iz++){
    //flip each site in the column
    for(int jz=0; jz<Ly; jz++){
      c = iz*Ly+jz;
      a = last[c];
      b = last[c + number_of_sites/2];
      last[c] = b;
      last[c+number_of_sites/2]=a;
    }
  }
  /*************************************************************
    End of switching the connections within region A
   *************************************************************/

  //now the RHS bondops
  for(bopnum; bopnum<number_of_bondops+number_of_sites; bopnum++){

    //figure out operator type
    //bond operator
    if(superbops(bopnum,2)==1){
      //join the first leg to the... 
      templegnum = nnbonds(superbops(bopnum,0),0);
      Hlinks[legnum] = last[templegnum];//does it matter if
      Hlinks[last[templegnum]] = legnum;//I screw up the order
      last[templegnum] = legnum + 2; //? because I am
      
      templegnum = nnbonds(superbops(bopnum,0),1);
      Hlinks[legnum+1] = last[templegnum]; 
      Hlinks[last[templegnum]] = legnum+1;
      last[templegnum] = legnum + 3;

      vleg2op[legnum/4] = bopnum;
      legnum+=4;
    }
    else{ //it's a plaquette operator!!!
      templegnum = plaquettes(superbops(bopnum,0),0);
      Hlinks[legnum] = last[templegnum];
      Hlinks[last[templegnum]] = legnum;
      last[templegnum] = legnum + 2;
      
      templegnum = plaquettes(superbops(bopnum,0),1);
      Hlinks[legnum+1] = last[templegnum]; 
      Hlinks[last[templegnum]] = legnum+1;
      last[templegnum] = legnum + 3;
      
      vleg2op[legnum/4] = bopnum;
      legnum +=4;
      
      //now the other side of the plaquette
      templegnum = plaquettes(superbops(bopnum,0),2);
      Hlinks[legnum] = last[templegnum];
      Hlinks[last[templegnum]] = legnum;
      last[templegnum] = legnum + 2;
      
      templegnum = plaquettes(superbops(bopnum,0),3);
      Hlinks[legnum+1] = last[templegnum]; 
      Hlinks[last[templegnum]] = legnum+1;
      last[templegnum] = legnum + 3;
      
      vleg2op[legnum/4] = bopnum+0.5;
      legnum +=4;  
    }
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
   nnbonds[#nnbonds][2] //
   superbops[#bops][3]  //
   Vlinks[vlegs]        //

 Local:
   loopnums[vlegs]_stores the loop number for each vertex leg
   loopnum_________the loop number we're currently looking at
   startleg________the start vertex for the current loop
   counter_________starts at the beginning of the vertex legs, goes to the end
   leg_____________the current vertex leg we're looking at
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
  int loopnum(1), startleg(0); 
  long long counter(0), leg(0);
  bool which(0), flip=0;
  int rfirstcross(-99),rlastcross(-99), rcurrent(0); 
  int right=-99;
  bool boolcross=0;

  VL.assign(number_of_sites, -99);
  VR = VL;
  whichloop = VL;
  cross = 0;

  while(counter < vlegs-2){

    if(drand()<0.5){flip=1;}
    else{flip=0;}

    startleg = counter;        // set the initial site (startsite)
    leg = Hlinks[counter];     // the site connected to startsite horizontally

    //making sure site isn't bonded to itself (edges are bonded to themselves)
    //and checking that the startsite isn't already in a loop
    while(((leg == Hlinks[leg])|(loopnums[counter]>0))&&(counter<vlegs-2)){ 
      counter++;  //changing startsite
      leg = Hlinks[counter]; //site connected to new startsite
      startleg = counter;    //setting new startsite
    }
    //breaks if we get to the last or 2nd last site
    if(counter > vlegs-2){break;} 
    loopnums[startleg] = loopnum; //including startsite in new loop
    loopnums[leg] = loopnum;    //adding the next site to the loop
    which = 0;      //"which" checks if we're on horiz(1) or vert(0) 
    

    if(leg>vlegMiddle&&startleg<vlegMiddle){
      //if it's a bond operator
      int tempOp = floor (vleg2op[leg/4]);
      if(superbops(tempOp,2)==1){
	rfirstcross = nnbonds(superbops(tempOp,0),leg%2);
      }
      //if it's a plaquette operator
      else{
	rfirstcross = plaquettes(superbops(tempOp,0),
				 (floor ((vleg2op[leg/4]-tempOp)+0.6))*4+leg%2);
      }
      right = 1;
      boolcross++;
      whichloop[rfirstcross]=loopnum;
      rlastcross = rfirstcross;
    }
    if(leg<vlegMiddle&&startleg>vlegMiddle){
      //if it's a bond operator
      int tempOp = floor (vleg2op[startleg/4]);
      if(superbops(tempOp,2)==1){
	rfirstcross = nnbonds(superbops(tempOp,0),leg%2);
      }
      //if it's a plaquette operator
      else{
	rfirstcross = plaquettes(superbops(tempOp,0),
				 (floor ((vleg2op[startleg/4]-tempOp)+0.6))*4+startleg%2);
      }

      right = 0;
      boolcross++;
      whichloop[rfirstcross]=loopnum;
      rlastcross = rfirstcross;
    }
  
    //while loop ends when we get back to the startsite
    while(leg!=startleg){
      
      //VERTICAL LINKS
      if(!which){                
	leg = Vlinks[leg];

	if(flip){
	  // superbops(leg/4,1) = (superbops(leg/4,1)+1)%2;
	  
	  //if it's a bond operator
	  int tempOp = floor (vleg2op[leg/4]);
	  if(superbops(tempOp,2)==1){ 
	  superbops(tempOp,1)=(superbops(tempOp,1)+1)%2;
	  }
	  //if it's a plaquette operator
	  else{
	    //figure out which of the two ops we're on
	    //if it's the second one
	    if(vleg2op[leg/4]-tempOp){
	      //add 2 mod 4
	      // 0->2,1->3,3->1,2->0
	      superbops(tempOp,1)=(superbops(tempOp,1)+2)%4;
	    }
	    //otherwise 0->1, 1->0, 2->3, 3->2
	    else{superbops(tempOp,1)=Vlinks[superbops(tempOp,1)];}
	  } 
	}
      }
      //HORIZONTAL LINKS
      //****************************
      else{
	//if it crosses the boundary
	if((leg<vlegMiddle)^(Hlinks[leg]<vlegMiddle)){
	  if(leg<vlegMiddle){

	    //if it's a bond operator
	    int tempOp = floor (vleg2op[Hlinks[leg]/4]);
	    if(superbops(tempOp,2)==1){
	      rcurrent = nnbonds(superbops(tempOp,0),Hlinks[leg]%2);
	    }
	    //if it's a plaquette operator
	    else{
	      rcurrent = plaquettes(superbops(tempOp,0),
				    (floor ((vleg2op[Hlinks[leg]/4]-tempOp)+0.6))*4
				    +Hlinks[leg]%2);
	    } 
	    //    rcurrent = nnbonds(superbops(Hlinks[leg]/4,0),Hlinks[leg]%2);
	  }
	  else{
	    
	    //if it's a bond operator
	    int tempOp = floor (vleg2op[leg/4]);
	    if(superbops(tempOp,2)==1){
	      rcurrent = nnbonds(superbops(tempOp,0),leg%2);
	    }
	    //if it's a plaquette operator
	    else{
	      rcurrent = plaquettes(superbops(tempOp,0),
				    (floor ((vleg2op[leg/4]-tempOp)+0.6))*4
				    +leg%2);
	    } 
	    //	    rcurrent = nnbonds(superbops(leg/4,0),leg%2);
	  }
	  boolcross++;

	  //if this is the first crossing for this loop set firstcross
	  if(rlastcross<0){//and set right to show if it's crossing left or right
	    rfirstcross = rcurrent;
	    if (leg<vlegMiddle){right=1;}
	    else {right=0;}
	  }

	  //if it's not the first crossing
	  else{
	    if(right){
	      VR[rcurrent] = rlastcross;
	      VR[rlastcross] = rcurrent;
	    }
	    else{
	      VL[rcurrent] = rlastcross;
	      VL[rlastcross] = rcurrent;
	    }
	    right = (right+1)%2;  //change cross direction
	  }
	  rlastcross = rcurrent;
	  whichloop[rcurrent]=loopnum;
	}
	
	leg = Hlinks[leg];
	
      }
      which = !which; //changes from horiz(1) to vert(0) or vice versa
      loopnums[leg] = loopnum; //adds next site to loop
    }
    if(rlastcross>-1){
      if(right){ VR[rfirstcross]=rlastcross; VR[rlastcross]=rfirstcross;}
      else{ VL[rfirstcross]=rlastcross; VL[rlastcross]=rfirstcross;}
    }
    if(boolcross){cross++;}
    boolcross = 0;
    loopnum++; //loop is finished, go to next loop
    rlastcross = -99;  right = -99;
    //if counter is already in a loop increase counter
    while(loopnums[counter]>0){counter++; if(counter>vlegs){break;}}
  }

  //check of states
  //  for(int i=0;i<number_of_sites;i++){
  //    cout << "VL["<<i<<"] = "<<VL[i]<<endl;
  //  }

  //Hlinks isn't used in any other functions.  Clear it.
  Hlinks.clear();
  vleg2op.clear();
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
  spins = init_spins;

  int neighbs(0);
  if(Ly==1){neighbs=2;}
  else{neighbs=6;}

  //for the first N/2 operators (i.e. the edge)
  for(int op=0; op<number_of_sites/2; op++){ 
    if(superbops(op,1)==1){                   //if operator is offdiagonal
      //update spins
      spins[nnbonds(superbops(op,0),0)] = (spins[nnbonds(superbops(op,0),0)]+1)%2;
      spins[nnbonds(superbops(op,0),1)] = (spins[nnbonds(superbops(op,0),1)]+1)%2;
    }                       //otherwise (if diagonal) do nothing
  }

  // Now look at the first half of the *real* operators
  long long op = number_of_sites/2;
  int otemp(-99);
  int oldtype=-1;
  bool flag = false;
  for(op; op<number_of_bondops/2+number_of_sites/2; op++){
    //if operator is offdiagonal
    if(superbops(op,1)==1){                       
      //update spins
      spins[nnbonds(superbops(op,0),0)] = (spins[nnbonds(superbops(op,0),0)]+1)%2;
      spins[nnbonds(superbops(op,0),1)] = (spins[nnbonds(superbops(op,0),1)]+1)%2;
    }     
    //if the operator is diagonal we need to change it randomly 
    else{                
    
      // int temp(-99);
      //do{temp = irand()%number_of_nnbonds;}
      //while(spins[nnbonds(temp,0)]+spins[nnbonds(temp,1)]!=1);
      //superbops(op,0) = temp; 
      flag = false;
      
      while(!flag){
	if(drand()<plaqProb){
	  //choose a plaquette!!
	  otemp = irand()%numPlaquettes;
	  //if the spins are antiparallel, assign it
	  if(spins[plaquettes(otemp,0)]+spins[plaquettes(otemp,1)]==1 && 
	     spins[plaquettes(otemp,2)]+spins[plaquettes(otemp,3)]==1){
	    superbops(op,0) = otemp;
	    superbops(op,1) = 0;
	    oldtype = superbops(op,2);
	    superbops(op,2) = 2; //2 for plaquette, right?

	    flag = true;
	    num2site += 2-oldtype;
	    middle += 2-oldtype;
	  }
	}
	else{
	  //choose a bondop
	  otemp = irand()%number_of_nnbonds;
	  //if the spins are antiparallel, assign it
	  if(spins[nnbonds(otemp,0)]+spins[nnbonds(otemp,1)]==1){
	    superbops(op,0) = otemp;
	    superbops(op,1) = 0;
	    oldtype = superbops(op,2);
	    superbops(op,2) = 1; //1 for bond operator
	    
	    flag = true;
	    num2site += 1-oldtype;
	    middle += 1-oldtype;
	  }
	}
      }
    }
  }

  vlegMiddle = 2*number_of_sites + 4*(middle+0.5) - 0.5;
  /*******************************************************************
               Swap some of the spins and stuff, y'know?
  ********************************************************************/
  long long a,b,c;
  
  //flip 1 column of Ly at a time
  for(int iz=0; iz<=flip; iz++){
    //flip each site in the column
    for(int jz=0; jz<Ly; jz++){
      c = iz*Ly+jz;
      a = spins[c];
      b = spins[c + number_of_sites/2];
      spins[c] = b;
      spins[c+number_of_sites/2]=a;
    }
  }
  /*******************************************************************
                               end of that
  ********************************************************************/
  
  for(op; op<number_of_bondops+number_of_sites/2; op++){
    if(superbops(op,1)==1){                         //if operator is offdiagonal
      //update spins
      spins[nnbonds(superbops(op,0),0)] = (spins[nnbonds(superbops(op,0),0)]+1)%2;
      spins[nnbonds(superbops(op,0),1)] = (spins[nnbonds(superbops(op,0),1)]+1)%2;
    }      //if the operator is diagonal we need to change it randomly

    else{                  
      // int temp(-99);
      //   do{temp = irand()%number_of_nnbonds;}
      //   while(spins[nnbonds(temp,0)]+spins[nnbonds(temp,1)]!=1);
      //  superbops(op,0) = temp;     

      flag = false;
      
      while(!flag){
	if(drand()<plaqProb){
	  //choose a plaquette!!
	  otemp = irand()%numPlaquettes;
	  //if the spins are antiparallel, assign it
	  if(spins[plaquettes(otemp,0)]+spins[plaquettes(otemp,1)]==1 && 
	     spins[plaquettes(otemp,2)]+spins[plaquettes(otemp,3)]==1){
	    superbops(op,0) = otemp;
	    superbops(op,1) = 0;
	    oldtype = superbops(op,2);
	    superbops(op,2) = 2; //2 for plaquette, right?

	    flag = true;
	    num2site += 2-oldtype;
	  }
	}
	else{
	  //choose a bondop
	  otemp = irand()%number_of_nnbonds;
	  //if the spins are antiparallel, assign it
	  if(spins[nnbonds(otemp,0)]+spins[nnbonds(otemp,1)]==1){
	    superbops(op,0) = otemp;
	    superbops(op,1) = 0;
	    oldtype = superbops(op,2);
	    superbops(op,2) = 1; //1 for bond operator
	    // cout << "oldtype = " << oldtype << endl;

	    flag = true;
	    num2site += 1-oldtype;
	  }
	}
      }   
    }
  }
  //cout << "vlegs1" << vlegs << endl;
  //cout << num2site << endl;
  vlegs = 4*number_of_sites + 4*num2site;
  // cout << vlegs << endl;
}
/************ swaperator() ****************************************************
what geometry does this even use?  squares?

i should change it to ladder geometry... but of course... bus error.
******************************************************************************/
void LOOPS::swaperator()
{
  vector <int> tempbonds;
  tempbonds = VR;
  int a,b,c,d;

  
  //flip 1 column of Ly at a time
  for(int iz=flip+1; iz<Lx; iz++){
    
    //flip each site in the column
    for(int jz=0; jz<Ly; jz++){
      a = iz*Ly+jz;
      d = a+number_of_sites/2;
      b = tempbonds[d];
      c = tempbonds[a];

      tempbonds[a] = b;
      tempbonds[b] = a;
      tempbonds[d] = c;
      tempbonds[c] = d;
      
    }

    //take measurement for each flipped column
    int counter(0), temploopnum(0), startsite(0), mite(-99), which(0);
    vector <int> site(number_of_sites+2,0);
    
    while(counter < number_of_sites){
      
      site[counter]=1;
      startsite = counter;
      
      mite = VL[counter];
      which=0;
      
      while(mite!=startsite){
	
	if(mite==-99){
	  cout << "SUPER ERROR in Swaterator: VL["<<counter<<"] = "<<mite << endl; 
	  exit(1); 
	}

	site[mite]=1;
	
	if(which==0){
	  mite = tempbonds[mite];
	  which++;
	}
	else{
	  mite = VL[mite];
	  which--;
	}
      }
      temploopnum++;
      while(site[counter]==1){counter++;}
    }
    int loopdiff = temploopnum - cross;
    entropy[iz] += pow(2,loopdiff);
    
  }
}
/************ calculate_stuff() ***********************************************
******************************************************************************/
void LOOPS::calculate_stuff()
{
  energy = 0.5*(energyint*0.75)/iterations;
  energyint = 0;

  for(int i=0; i<entropy.size(); i++){
    //  entropy_final[i] = -log(entropy[i]/(1.0*iterations));
    entropy_final[i] = (entropy[i]/(1.0*iterations));
    entropy[i]=0;
  }
}

/************ print_bops() ****************************************************
******************************************************************************/
void LOOPS::print_bops()
{
  ofstream bout(bopfile.c_str());
  bout << Lx << endl;
  bout << Ly << endl;
  for(int i=0; i<number_of_bondops+number_of_sites; i++){
    bout << superbops(i,0) << endl << superbops(i,1) << endl << superbops(i,2) << endl;
  }
  
  bout << -99 << endl;
}

/************ read_bops() *****************************************************
| If the bond operator file (filename given in parameters file)
******************************************************************************/
void LOOPS::read_bops()
{
  ifstream bin(bopfile.c_str());
  
  if(bin.fail()){ generate_ops(); }
  else{ 
    int test=0;
    int tally(0);
    int midtally(0);
    
    //Check Lx value
    bin >> test;
    if(test!=Lx){cout<<"wrong Lx dimension!!\n"; exit(1);}

    //Check Ly value
    bin >> test;
    if(test!=Ly){cout<<"wrong Ly dimension!!\n"; exit(1);}

    //Read in bondOps
    superbops.resize(number_of_bondops+number_of_sites,3);

    for(int i=0; i<number_of_bondops+number_of_sites; i++){
      bin >> superbops(i,0) >> superbops(i,1) >> superbops(i,2);

      if(i>=number_of_sites/2 && i<number_of_sites/2+number_of_bondops){
	tally+=superbops(i,2);
	if(i<number_of_sites/2+number_of_bondops/2){
	  midtally+=superbops(i,2);
	}
      }
    }
    
    middle = midtally-0.5;
    //  cout << middle << endl;
    num2site = tally;
    //  cout << num2site << endl;
    
    vlegs = 4*number_of_sites + 4*num2site;
    //  cout << vlegs << endl;
    vlegMiddle = 2*number_of_sites + 4*(middle+0.5) - 0.5;
    //   cout << vlegMiddle << endl;

    //Read in -99 at EOF
    bin >> test;
    if(test!=-99){
      cout << "Problem with bondop file! Wrong system size???" << endl; 
      exit(1);
    }  
  }
}
#endif
