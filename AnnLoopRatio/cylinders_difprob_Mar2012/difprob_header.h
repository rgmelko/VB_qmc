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


#ifndef difprob_header
#define difprob_header

#include "header.h"
#include "matrix.h"
#include <vector>

class LOOPS
{
 public:

  MTRand drand; //drand() gives you a random double precision number
  MTRand_int32 irand; // irand() gives you a random integer

  int flip; // -1 for bare swap, 0 for ratio 1, etc

  int Lx, Ly,  number_of_sites; //the dimensions and number of sites
  //Note: things are periodic in the Ly direction, open in Lx

  int cross; /*the number of loops crossing the boundary i.e. the number of 
	       loops created by overlapping the propagated |VL> and |VR> */
  long long energyint; /*keeps track of the energy: 
			    energy = energyint*0.75/iterations */
  long long number_of_bondops, number_of_nnbonds;
  long long vlegs;/*the number of vertex legs including legs from the vertices
		    used to simulate the edge states |VL> and |VR> */
  long long iterations; //number of iterations per loop. Used for energy calc.
  double energy; //the energy in non-integer form
  string bopfile; //the name of the file in which the bondops are stored
  
  vector <long long> Vlinks, Hlinks;//the vert and horizontal links for the LL
  vector <int> spins, init_spins; //keeps track of spins for swaperation
  vector <int> antipar, init_antipar; //keeps track of antiparallelness
  vector <int> isgood, init_isgood; //keeps track of 'good'ness
  vector <int> sides; /*store which side of the boundary a leg is on
			removing the need to check if it's higher or lower
			that the "middle" of the number of legs.... I'm not
			sure if this is more efficient or not */
  vector <int> VL, VR; //the propagated left and right VB states
  vector <int> whichloop; /* stores which loop number each site is in. used
			     in the energy measurement */
  vector <double> entropy, entropy_final;

  iMatrix nnbonds;//list of all possible nnbonds.  Index is the bond number
  iMatrix nn_mat; /*matrix of the nnbonds. indices are sites
		    and contents are bond numbers.*/
  iMatrix Nnnbonds;//list of all nearest nnbonds
  iMatrix bops; // list of bond operators
  iMatrix superbops; //list of bond operators plus edges simulated via bops
                     //superbops(:,0) is the bond the operator acts on 
                     //         (:,1) is 0 for diag, 1 for offdiag

  //CONSTRUCTOR
  LOOPS(int xsites, int ysites, int flips, int bondops, long long its,
	       long long rseed, string bondopfile);

  void nnbondlist(); //creates list of nnbonds
  void Nnnbondlist(); //creates list of Nnnbonds
  void generate_ops(); //generates initial operators
  void create_Vlinks(); //creates vertical links
  void create__Hlinks(); //creates horizontal linkts (harder)
  void make_flip_loops(); //generates and flips loops (w/ prob 1/2)
  void take_measurement(); //measures energy at the moment
  void change__operators(); //changes the diagonal operators randomly

  void swaperator(); //the swap operator

  void calculate_stuff(); //calculates energy at the moment
  void print_bops(); //prints the bond operators after every loop of iterations
  void read_bops(); //reads in bondops if theyre there, or runs generate_ops()

};

//*************** CONSTRUCTOR ******************************************
LOOPS::LOOPS(int xsites, int ysites, int flips, int bondops, long long its, 
	     long long rseed, string bondopfile)
{
  irand.seed(rseed); //uses the random seed from the parameter file
  drand.seed(rseed); //enabling us to run fakely parallelize simulations
  
  Lx = xsites; 
  Ly = ysites;
  flip = flips;
  number_of_sites = Lx*Ly; //calculates total number of sites
  //****changed**** multiplied by 2
  number_of_bondops = 2*2*bondops; /*the *real* number of bondops is multiplied
				   by 2, one set for |VL> and one for |VR> */
    //****changed**** multiplied first term by 2
  vlegs = 2*4*number_of_sites + 4*number_of_bondops; //number of vertex legs
  energyint = 0; energy = 0; //initialize the energy counters
  iterations = its; 
  bopfile = bondopfile; //name of the bond operator file

  entropy.assign(xsites,0);
  entropy_final = entropy;

  Vlinks.assign(vlegs, -99); //set size and initialize
  Hlinks.assign(vlegs, -99); 
  //Initialize Hlinks so the "edge" sites are linked to themselves
  for(long long i=0; i<vlegs; i++){ Hlinks[i]=i; } 

  bops.resize(number_of_bondops,2); //set size of bops
  //create "sides"
  sides.assign(vlegs,0); 
  for(long long i=vlegs/2; i<vlegs; i++){sides[i]=1;}
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

/********** nnbondlist() ***********************************************
   Uses:
     Global:
         Ly                     //
         number_of_nnbonds______sets value based on 1D/2D and PBC/OBC
         number_of_sites        //
         nnbonds[#nnbonds][2]___sized and filled
         nn_mat[#sites][#sites]_sized and filled
         Lx                     //
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
  // should put in a special case for Ly==2.

  //Case of 1D open chain
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
    
    //changed multiplied first dimension by 2
    nnbonds.resize (2*number_of_nnbonds,2);
    for(int i=0;i<number_of_nnbonds*2;i++){nnbonds(i,0)=-99;nnbonds(i,1)=-99;}
    //resize and initialize the matrix of nnbonds
    //changed multiplied dimensions by 2
    nn_mat.resize(2*number_of_sites, 2*number_of_sites);
    for(int i=0; i<2*number_of_sites; i++){
      for(int j=0; j<2*number_of_sites; j++){
	nn_mat(i,j) = -99;
      }
    }

    //fill nnbonds and nn_mat
    //ybonds!!
    int counter = 0;
    for(int x=0; x<Lx; x++){
      for(int y=0; y<Ly; y++){
	nnbonds(counter,0) = Ly*x+y;
	nnbonds(counter,1) = Ly*x+(y+1)%Ly;
	nn_mat(nnbonds(counter,0),nnbonds(counter,1))=counter;
	nn_mat(nnbonds(counter,1),nnbonds(counter,0))=counter;
	counter++;
      }
    }

    //xbonds...
    for(int i=0; i<(Lx-1)*Ly; i++){
      nnbonds(counter,0) = i;
      nnbonds(counter,1) = i+Ly;
      nn_mat(nnbonds(counter,0),nnbonds(counter,1))=counter;
      nn_mat(nnbonds(counter,1),nnbonds(counter,0))=counter;
      counter++;
    }

    // make sure we have the proper number of nnbonds
    if(counter!=number_of_nnbonds){cout << "supererror" << endl;}

    //OTHERCHECKS
    cout<<"first write them... test... then comment out OTHERCHECKS:\n";
   
  }

  //changed added this part in to double nnbonds and nnmat
      for(int b=number_of_nnbonds; b<2*number_of_nnbonds; b++){
	int a = number_of_nnbonds;
	nnbonds(b,0) = nnbonds(b-a,0)+number_of_sites;
	nnbonds(b,1) = nnbonds(b-a,1)+number_of_sites;
	nn_mat(nnbonds(b,0),nnbonds(b,1))=b;
	nn_mat(nnbonds(b,1),nnbonds(b,0))=b;
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


      //changed now multiplying #nnbonds by 2
      number_of_nnbonds *=2; 
      
      
      // resizing antiparallelness vectors
      init_antipar.assign(number_of_nnbonds, 0);
      //  init_isgood.resize(number_of_nnbonds);
      //changed it to *check* antiparallelness (since there can be odd dims now)
      for(int i=0; i<number_of_nnbonds; i++){
	if(spins[nnbonds(i,0)]+spins[nnbonds(i,1)]==1){
	  init_antipar[i]=1; 
	  init_isgood.push_back(i);
	}
	else{ init_antipar[i]=0;}
      }
      cout << "number of antiparallel bonds: " << init_isgood.size()<< endl;
  
      antipar.resize(number_of_nnbonds);
      isgood.resize(init_isgood.size());
      antipar = init_antipar;
      isgood = init_isgood;
      
      //changed**** changing #nnbonds back now
      number_of_nnbonds /=2;


      
}


/************* Nnnbondlist() *************************************************
 Uses:
  Global:
   Ly                   //
   Nnnbonds[#Nnnbonds][2] //sizes and fills
   number_of_nnbonds      //
   nn_mat[#sites][#sites] //
   number_of_sites        //

  Local:
   counter__goes from 0-3 (i think) and represents the number of neighbours
            we've found so far.
   bnum_____the number of the bond we're looking at

   Creates a list of the neighbouring bonds for a given bond.
   Works for  PBC, 1D, and 2D.
   For 2D, uses the matrix nn_mat, which is created in nnbondlist().
*****************************************************************************/
void LOOPS::Nnnbondlist()
{
  // There's no option for Lx==1 FYI.  So no open chain.

 
  int num_neighbs = 0;
  if(Ly==1){  // 1D open chain case
    num_neighbs = 2;
    //changed**** multiplied first dimension by 2
    Nnnbonds.resize(2*number_of_nnbonds,num_neighbs);

    // generate nearest nnbonds
    for(int i=0; i<number_of_nnbonds; i++){
      Nnnbonds(i,0)=(i+1)%number_of_nnbonds; // forward bond
      Nnnbonds(i,1)=(i+number_of_nnbonds-1)%number_of_nnbonds; //backward bond
    }
    // It's always OBC for Lx now!!!
    // if(OBC){  // For open BCs first and last bonds only have 1 neighbour
    Nnnbonds(0,1) = -99; //no backward bond for 1st bond
    Nnnbonds(number_of_nnbonds-1,0) = -99; //no forward bond for last bond
    // }
  }
  
  else{  // the 2D case... more complicated...
    num_neighbs = 6;
    // Resize and initialize Nnnbonds
    //changed**** multiplied first dimension by 2;
    Nnnbonds.resize(2*number_of_nnbonds,num_neighbs);
    for(int i=0; i<2*number_of_nnbonds; i++){
      for(int j=0; j<num_neighbs; j++){
	Nnnbonds(i,j)=-99;
      }
    }
    
    //This looks at the full row and column of the nn_mat corresponding
    //to bond number i.
    //The matrix is symmetric, so you could equivalently go through the
    //nnbonds(i,0) and nnbonds(i,1) rows (or columns even)
    for(int i=0; i<number_of_nnbonds; i++){
      int counter=0;
      for(int j=1; j<number_of_sites; j++){
	
	//going through the nn matrix
	int bnum = nn_mat(nnbonds(i,0), (nnbonds(i,1) + j)%number_of_sites);
	//if the bond is a nn of the initial bond;
	if(bnum != -99){Nnnbonds(i,counter) = bnum; counter++;}
	

	//doing the same but moving down instead of across
	bnum = nn_mat((nnbonds(i,0) + j)%number_of_sites, nnbonds(i,1));
	if(bnum != -99){Nnnbonds(i,counter) = bnum; counter++;}
      }
    } 
  } 


  //changed**** added this in to copy Nnnbonds.. to double it.. sorta
  for(int q=number_of_nnbonds; q<2*number_of_nnbonds; q++){
    for(int s=0; s<num_neighbs; s++){
      int t=number_of_nnbonds;
      if(Nnnbonds(q-t,s)==-99){Nnnbonds(q,s)=-99;}
      else{Nnnbonds(q,s)=Nnnbonds(q-t,s)+number_of_nnbonds;}
    }
  }



  //  for(int i=0;i<Nnnbonds.length();i++){
  //    for(int j=0;j<Nnnbonds.width();j++){
  //      cout << Nnnbonds(i,j) << " ";
  //    }
  //    cout << endl;
  //  }

  //changed**** for realisies doubling #nnbonds and #sites
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
    
  //changed so it doesn't just pick random bonds... the spins must be antiparallel
  int temp(-99);
  
  for(int i=0; i<number_of_bondops; i++)
    {
      //  bops(i,0) = init_isgood[irand() % init_isgood.size()]; //the bond being operated on
      //  bops(i,1) = 0; //0 = diagonal, 1 = off-diagonal
      //  cout << "bops("<<i<<")="<<bops(i,0)<< endl;

      do{temp = irand()%number_of_nnbonds;}
      while(spins[nnbonds(temp,0)]+spins[nnbonds(temp,1)]!=1);
      bops(i,0) = temp;
      bops(i,1) = 0;
    }

  superbops.resize(number_of_bondops+number_of_sites,2);
  //  for(int i=0;i<number_of_sites+number_of_bondops;i++){
  //    superbops(i,0)=-99;
  //    superbops(i,1)=-99;
  //    }
     
  //initial state is dimerized, but either Lx or Ly can be odd.
  //Dimerize in x direction if Lx is even
  //Otherise dimerize in Ly direction
  int bondd=0;
  if(Lx%2==0){
    for(int ix=0; ix<Lx; ix+=2){
      for(int iy=0; iy<Ly; iy++){
	
	//	cout << iy+ix*Ly<< ", "<<iy+(ix+1)*Ly<<"    "<<Lx*Ly+iy+ix*Ly<<", "<<Lx*Ly+iy+(ix+1)*Ly<<endl;

	superbops(bondd,0) = nn_mat(iy+ix*Ly,iy+(ix+1)*Ly);
	superbops(bondd,1) = 0;
	superbops(bondd+Lx*Ly/2,0) = nn_mat(Lx*Ly+iy+ix*Ly,Lx*Ly+iy+(ix+1)*Ly);
	superbops(bondd+Lx*Ly/2,1) = 0;
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
	superbops(bondd+Lx*Ly/2,0) = nn_mat(Lx*Ly+iy+ix*Ly,Lx*Ly+iy+1+ix*Ly);
	superbops(bondd+Lx*Ly/2,1) = 0;
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
  }

  for(int i=0; i<number_of_bondops; i++){
    superbops(number_of_sites/2+i,0)=bops(i,0);
    superbops(number_of_sites/2+i,1)=bops(i,1);
  }

  //  for(int i=number_of_bondops;i<number_of_bondops+number_of_sites; i++){
  //    cout << "superbops("<<i<<",0) = "<<superbops(i,0)<<endl;
  //  }
  
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


i assume vertex legs are labelled like so (Feb 2012):

                0 ----------- 2
                       |
                       |
                       |
                1 ----------- 3

*******************************************************************************/
void LOOPS::create__Hlinks()
{
  //what does this part do???
  //might be implying an initial dimerization!!!
  //yeah totally.  fix it up!
  //I think it's the same problem for the end too!!
  vector <long long> last (number_of_sites,-99);

  //put initial state as initial vertices in linked list.
  //THIS IS CONFUSING
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

  // the first N/2 bondops are actually the initial VB configuration
  long long bopnum = number_of_sites/2;

  // create Hlinks for the operators after the init VB config.. 
  // but also include the bondops for the other initial VB config 
  // (on the right side) ... 

  //for the LHS of the bondops
  for(bopnum; bopnum<(number_of_bondops+number_of_sites)/2; bopnum++){

    //the top left (first) leg of the bond op
    legnum = 4*bopnum;


    //join the first leg to the 
    Hlinks[legnum] = last[nnbonds(superbops(bopnum,0),0)];//does it matter if
    Hlinks[last[nnbonds(superbops(bopnum,0),0)]] = legnum;//I screw up the order
    last[nnbonds(superbops(bopnum,0),0)] = legnum + 2; //? because I am

    Hlinks[legnum+1] = last[nnbonds(superbops(bopnum,0),1)]; 
    Hlinks[last[nnbonds(superbops(bopnum,0),1)]] = legnum+1;
    last[nnbonds(superbops(bopnum,0),1)] = legnum + 3;
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

    //the top left (first) leg of the bond op
    legnum = 4*bopnum;

    //join the first leg to the 
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
    loopnums[startsite] = loopnum; //including startsite in new loop
    loopnums[site] = loopnum;    //adding the next site to the loop
    which = 0;      //"which" checks if we're on horiz(1) or vert(0) 
    

    /*    if(sides[startsite]!=sides[site]){ 
	  cout << "cross!!!\n";
	  firstcross =  nnbonds(superbops(site/4,0),site%2);//************
	  lastcross = firstcross;
	  cout << "firstcross = " << firstcross << endl;
	  right = sides[site];
	  boolcross++;
	  whichloop[firstcross]=loopnum;  
      }*/

    if(sides[startsite]<sides[site]){
      rfirstcross = nnbonds(superbops(site/4,0),site%2);
      right = 1;
      boolcross++;
      whichloop[rfirstcross]=loopnum;
      rlastcross = rfirstcross;
    }
    if(sides[startsite]>sides[site]){
      rfirstcross = nnbonds(superbops(startsite/4,0),startsite%2);
      right = 0;
      boolcross++;
      whichloop[rfirstcross]=loopnum;
      rlastcross = rfirstcross;
    }
  
    //while loops ends when we get back to the startsite
    while(site!=startsite){
      
      //VERTICAL LINKS
      if(!which){                
	site = Vlinks[site];

	if(flip){superbops(site/4,1) = (superbops(site/4,1)+1)%2;}
      }
      //HORIZONTAL LINKS
      //****************************
      else{
    
	if(sides[site]!=sides[Hlinks[site]]){  //if it crosses the boundary
	  if(sides[site]<sides[Hlinks[site]]){
	    rcurrent = nnbonds(superbops(Hlinks[site]/4,0),Hlinks[site]%2);
	  }
	  else{
	    rcurrent = nnbonds(superbops(site/4,0),site%2);
	  }
	  boolcross++;

	  //if this is the first crossing for this loop set firstcross
	  if(rlastcross<0){//and set right to show if it's crossing left or right
	    rfirstcross = rcurrent;
	    if (sides[site]<sides[Hlinks[site]]){right=1;}
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
	
	site = Hlinks[site];
	
      }
      which = !which; //changes from horiz(1) to vert(0) or vice versa
      loopnums[site] = loopnum; //adds next site to loop
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
  //  antipar = init_antipar;
  // isgood = init_isgood;
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



      /*    for(int i=0;i<neighbs;i++){//change antiparallelness of neighboring bonds
	    
	    int loc = Nnnbonds(superbops(op,0),i);
	    
	    if(loc<0){continue;}//goto start of loop if its OBC&that nn doesnt exist
	    if(antipar[loc]==1){//if bond is already antiparallel change to parallel
	    
	    int i=0;
	    
	    while(isgood[i]!=loc){ i++;};
	    isgood.erase(isgood.begin() + i);
	    antipar[loc]--;
	    }
	    else{
	    antipar[loc]++;   //if bond is parallel change to antiparallel
	    isgood.push_back(loc);
	    }                                        
	    }
      */


    }                       //otherwise (if diagonal) do nothing
  }

  // Now look at the first half of the *real* operators
  long long op = number_of_sites/2;
  for(op; op<number_of_bondops/2+number_of_sites/2; op++){
    if(superbops(op,1)==1){                         //if operator is offdiagonal
      //update spins
      spins[nnbonds(superbops(op,0),0)] = (spins[nnbonds(superbops(op,0),0)]+1)%2;
      spins[nnbonds(superbops(op,0),1)] = (spins[nnbonds(superbops(op,0),1)]+1)%2;



      /* for(int i=0;i<neighbs; i++){//change antiparallelness of neighboring bonds
	 
	 int loc = Nnnbonds(superbops(op,0),i);
	 
	 if(loc<0){continue;}//goto start of loop ifOBC&that nn doesnt exist
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
      */



    }      //if the operator is diagonal we need to change it randomly
    else{       //using whichbond..              
      //   superbops(op,0) = isgood[irand()%isgood.size()];
      int temp(-99);
      do{temp = irand()%number_of_nnbonds;}
      while(spins[nnbonds(temp,0)]+spins[nnbonds(temp,1)]!=1);
      superbops(op,0) = temp;
      
    }
  }

  /*******************************************************************
               Swap some of the spins and stuff, y'know?
  ********************************************************************/
  //swap the spins

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

  /*************************************************************
   *************************************************************/
  
  //find the new antiparallelness
  /*
    isgood.clear();
    for(int i=0; i<number_of_nnbonds; i++){
    //if they're different it's 1, otherwise 0)
    antipar[i]=(spins[nnbonds(i,0)]+spins[nnbonds(i,1)])%2;
    if(antipar[i]==1){isgood.push_back(i);} //if it's antiparallel add it to the list
    }
  */
  
  /*******************************************************************
                               end of that
  ********************************************************************/
  
  for(op; op<number_of_bondops+number_of_sites/2; op++){
    if(superbops(op,1)==1){                         //if operator is offdiagonal
      //update spins
      spins[nnbonds(superbops(op,0),0)] = (spins[nnbonds(superbops(op,0),0)]+1)%2;
      spins[nnbonds(superbops(op,0),1)] = (spins[nnbonds(superbops(op,0),1)]+1)%2;


      /*      for(int i=0;i<neighbs; i++){//change antiparallelness of neighboring bonds
	      
	      int loc = Nnnbonds(superbops(op,0),i);
	      
	      if(loc<0){continue;}//goto start of loop ifOBC&that nn doesnt exist
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
      */


    }      //if the operator is diagonal we need to change it randomly
    else{       //using whichbond..              
      //superbops(op,0) = isgood[irand()%isgood.size()];
       int temp(-99);
       do{temp = irand()%number_of_nnbonds;}
       while(spins[nnbonds(temp,0)]+spins[nnbonds(temp,1)]!=1);
       superbops(op,0) = temp;     
    }
  }
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
