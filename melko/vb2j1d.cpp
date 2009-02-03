// vb2j1d.cpp  Last updated Jan 13, 2009
// Trying adapt vb2j.cpp back to one dimension

#include<iostream>
#include<math.h>
#include"mtrand.h" // random number generator


using namespace std;

void shuffle(int [][2]); //randomizes initial bonds (currently not working)
void print_chain(int chain [][2]); //prints the bonds
void generate_operator(int operater[2], int neighbours[][2], int Js[], int indexx);
//generates 1 bond operator
double apply_operator(int op0, int op1, int chain[][2], int Jayrock);
//applies 1 bond operator
void change_operators(int operaters[][2], int Js[], int a, int neighbours[][2]);/*randomly 
changes a number of bond operators... where that number is "a"*/

MTRand drand; //drand() gives you a random double precision number
MTRand_int32 irand; // irand() gives you a random integer


const int lattice_type = 0; // 0 for columnar, 1 for staggered
const long int superseed = 583409361; // ********You************************
const int L = 16; // 1-D length of the lattice *******Can********************
const int zone = 1; // the size of "the zone" *********Change***************
const double jprime =1.0; // ****************************These*Values*********
double J = 1.0;
const int L2 = L*L; // total number of sites
const int half_L = L/2; // total number of sites divided by 2
const int n = L*5; // number of bond operators
const int start = 10000000; /* number of iterations until the programs takes ***
			     measurements  */
const int iterations = 10*start; // total number of iterations
int chain [half_L][2] = {0}; // the bonds are stored in here
int operater[2] = {0}; //it's an operator
int initial_state[half_L][2] ={0,2,1,4,3,5,10,7,8,9,6,11,13,15,14,12}; 
//stores the initial bond configuration

int Js[n] = {0};   // Stores the interaction strength for each operator 0=J,1=J'
double w_old[2]={0};  // the new and old weights
double w_new[2]={0};  // the new and old weights
double x_old[2]={0};  // the new and old weights
double x_new[2]={0};  // the new and old weights
double probb = 0; // prob of keeping new operators

int neighbours[L][2]; //lists the 2 nearest neighbours for each site
int box[L] = {0}; // the zone! <-- exclamation point

main() // the main program..
{

  irand();
  // print_chain(initial_state);

  cout << "one dimensional spin chain" << endl;

  cout << "L = " << L << "    " << "zone = " << zone << "   " << 
    iterations << " iterations" << "    n = " << n;
  cout << "     J = " << J << "   J' = " << jprime << endl;
  cout.precision(10); // ten digits of precision..  or ten decimal places?

  //******Finding Nearest Neighbours********************************
      for(int iii=0; iii<L; iii++)
      {
	if(iii==0){neighbours[iii][0]=9999;}  //L-1;}   //North
	else neighbours[iii][0]=iii-1;
	if(iii==L-1){neighbours[iii][1]=9999;}  //=0;} //South
	else neighbours[iii][1]=iii+1;
	
	//    cout << neighbours[iii][0] << ", ";
	//    cout << neighbours[iii][1] << ", ";
	//    cout << endl;
      }

  //****Define the zone-box-type-thing*******************************
      for(int abox=0; abox<zone; abox+=1)
      {
	box[abox]=1;
      }
  //****Print out the zone box*****************************************
//       for (int zzz = 0; zzz<L2; zzz+=L)   
//       {
// 	for(int xxx = 0; xxx<L; xxx++)
// 	  {
// 	    cout << box[xxx+zzz] << ", ";
// 	  }
// 	cout << endl;
//       }
  //*************************************************************************
  //  cout << endl;

   //  shuffle(initial_state);//"randomize" the initial state (but not really)
 

  //  print_chain(initial_state);

      for(int supacount = 0; supacount < 100; supacount++)
      {

  int operaters[n][2], new_operaters[n][2];   // old and new operators     
  int bond[2] = {0};   //number of NN J bonds
  int bondprime[2] = {0};  // number of NN J' bonds
  int  acc = 0, rej = 0;         //number of changes accepted and rejected
  int cross[2] = {0};   //the number of bonds crossing the zone boundary
  double energy = 0;          // the energy
  double energyprime = 0;
  double entropy = 0;
  double jaybo = 0;

  //-------Generate Operators----------------------------------------
  for(int i0=0; i0<n; i0++)
    {  
      generate_operator(operater, neighbours, Js, i0);
      operaters[i0][0] = operater[0];
      operaters[i0][1] = operater[1];
    }
  //-----------------------------------------------------------------

  //--------Initialize Bonds------------------------------------------
  for(int jai=0; jai<half_L; jai++)
    {
      chain[jai][0] = initial_state[jai][0];
      chain[jai][1] = initial_state[jai][1];
    }
  //------------------------------------------------------------------

  w_old[0] = 0;
  w_old[1] = 0;  // this is really unnecessary because I already set
  x_old[0] = 0;  // these to zero and I haven't changed them
  x_old[1] = 0;
  w_new[0] = 0;
  w_new[1] = 0;
  x_new[0] = 0;
  x_new[1] = 0;
  

  //-------Apply Operators-------(also get the weight)------------------
  for(int k=0; k<n; k++)
    {
      apply_operator(operaters[k][0],operaters[k][1],chain,Js[k]); 
    }
  //-------------------------------------------------------------------

  w_old[0] = w_new[0];
  w_old[1] = w_new[1];
  x_old[0] = x_new[0];
  x_old[1] = x_new[1];
  w_new[0] = 0;
  w_new[1] = 0;          // now it is necessary
  x_new[0] = 0;
  x_new[1] = 0;

  int q = 0; // records the number of steps
  int q2 = 0; // used to show how far along the program is

  for (int i=0; i<iterations; i = i++)
    {  
      // if(i%start==0){cout<< q2 << "% "<< endl ;  q2+=10;}//can include this because i
      //don't have a % operator for bignum yet

      for(int i7=0; i7<n; i7++)
	{ 
	  new_operaters[i7][0] = operaters[i7][0];
	  new_operaters[i7][1] = operaters[i7][1];
	}

      change_operators(new_operaters,Js, n, neighbours);

       for(int jay=0; jay<half_L; jay++)
	{
	  chain[jay][0] = initial_state[jay][0];
	  chain[jay][1] = initial_state[jay][1];
	}

      for(int k=0; k<n; k++)
	{
	  apply_operator(new_operaters[k][0],new_operaters[k][1],chain,Js[k]);
	  // cout << Js[k] << "   " <<  new_operaters[k][0] << ", " << new_operaters[k][1] << endl;
	}


//       cout << endl;
//       print_chain(chain);

      probb =  pow(0.5,(w_new[0]+w_new[1]-w_old[0]-w_old[1]))*
	pow(J,(w_new[0]+x_new[0]-w_old[0]-x_old[0]))*
	pow(jprime,(w_new[1]+x_new[1]-w_old[1]-x_old[1]));

  
      //      cout << " prob = " <<  probb << endl;
//   	"wold0 = " << w_old[0] << "   wnew0 = " << w_new[0] << endl
//   	   <<"wold1 = " << w_old[1] << "   wnew1 = " << w_new[1] << endl
//   	   <<"xold0 = " << x_old[0] << "   xnew0 = " << x_new[0] << endl
//   	   <<"xold1 = " << x_old[1] << "   xnew1 = " << x_new[1] << endl;
	  

	if(drand() < probb)
	{
	  if(i >= start)
	    {
	      for(int id=0; id<half_L; id++)
		{
		  if(
		     (chain[id][1] == neighbours[chain[id][0]][0])|    // if bond operator is diagonal
		     (chain[id][1] == neighbours[chain[id][0]][1])    // spin chain is not changed
		        // but bond is counted
		     )
		    {
		      
		      
		      if( ((chain[id][0]+chain[id][1]-1)/2)%2 == 0 )  // condition for J bond
			{bond[1]+=1;}                                 // (instead of J')
		      else{bondprime[1]+=1;}
		      
		    }
		
		  if(box[chain[id][0]]+box[chain[id][1]]==1){cross[1]+=1;}
		}
	      energy += bond[1]*1.0;
	      energyprime += bondprime[1]*1.0;
	      entropy += cross[1]*1.0;
	      q=q+1;
	      acc=acc+1;
	      
	    }

	  for(int i8=0; i8<n; i8++)
	    {
	      operaters[i8][0] = new_operaters[i8][0];
	      operaters[i8][1] = new_operaters[i8][1];
	    }
	  w_old[0] = w_new[0];
	  w_old[1] = w_new[1];
	  x_old[0] = x_new[0];
	  x_old[1] = x_new[1];
	  bond[0] = bond[1];
	  bondprime[0] = bondprime[1];
	  cross[0] = cross[1];
	}

      else{if(i>=start)
	  {
	    if(q==0){} 
	    else
	      {
		energy += bond[0]*1.0;
		energyprime += bondprime[0]*1.0;
		entropy += cross[0]*1.0;
		q=q+1;
		rej=rej+1;
	      }
	  }
      }
      
      w_new[0] = 0;
      x_new[0] = 0;
      w_new[1] = 0;
      x_new[1] = 0;
      cross[1] = 0;
      bond[1] = 0;
      bondprime[1] = 0;
    }

  //cout << "bonds " <<  energy << "     bonds' "<< energyprime << "      q = " << q << endl;

  energy /= ((L-1)*q*1.0);   //all 'L's changed to L-1 due to OBCs
  energy += 0.25;             //is that right?
  energy *= -(L-1)*J*0.5; 
  energyprime /= ((L-1)*q*1.0);
  energyprime +=  0.25;
  energyprime *= -(L-1)*jprime*0.5; 
  energy = (energy + energyprime);

  //energy /= L*1.0;

  entropy /= q*1.;
  entropy *= log(2);
 
  //cout << endl;
  //cout << "energy = " << energy << endl;
  cout << energy << endl;
      } //******end of supacount loop*********

  //  cout << "entropy = " << entropy << endl;
  // cout << "entropy/x = " << entropy/zone << endl;
  // cout << q << " energies" << endl;
  //cout << (100.*acc)/(q*1.) << "% accepted"<< endl;


  return 0;
}

void shuffle(int chain[][2])
{

  for(int i=0; i<half_L; i++){chain[i][0] = -1;}
  for(int i=0; i<half_L; i++){chain[i][1] = -1;}

  int row=0;
  int col=0;
  int rowcol=0;
  int check = 0;

//   for (int site = 0; site < L2;)
//     {
//       rowcol =  (irand() + 2) %2;
//       if(rowcol==0)
// 	{
// 	  row =  (irand()+half_L) %half_L;
// 	  while((chain[row][0] != -1)&(rowcol==0))
// 	    { 
// 	      row = (irand()+half_L) %half_L; 
// 	      check += 1;
// 	      if(check==half_L*half_L*half_L){rowcol=1;}
// 	    }
// 	  check = 0;
// 	  if(rowcol==0){chain[row][0]=site; site++;}
// 	}
//       if(rowcol==1)
// 	{
// 	  col =  (irand()+half_L) %half_L;
      
// 	  while((chain[col][1] != -1)&(rowcol==1))
// 	    { 
// 	      col = (irand()+half_L) %half_L; 
// 	      check += 1;
// 	      if(check==half_L*half_L*half_L){rowcol=0;}
// 	    }
// 	  check = 0;
// 	  if (rowcol==1){chain[col][1]=site;site++;}
// 	}

  for(int site=0; site<(L2-1); site+=2)
    {
      chain[row][0] = site;
      chain[row][1] = site+1;
      row += 1; 
    }
}

void print_chain(int chain[][2])
{
  for (int i2 = 0; i2 < half_L; i2++)
    {
      cout << chain[i2][0] << ',' << chain[i2][1] << endl;
    }

  cout << endl;
}

void generate_operator(int operater[2], int neighbours[][2], int Js[], int k)
{
  Js[k] = 0;
  int initb = (irand()+L) %L;
  int neighb = (irand()+2)%2;

  //OPEN BCs -------------------------------------------------------------------------------------
  while((initb == 0 && neighb == 0) || (initb == L-1 && neighb == 1))
    {
      initb = (irand()+L) %L;
      neighb = (irand()+2)%2;
    }
  //END OF OPEN BCs --------(except for the part in change_operators)------------------------------*/
  
  operater[0] = initb;
  operater[1] = neighbours[initb][neighb];


  if((operater[0]+operater[1])%2 == 1){Js[k]=1;}  
  //  cout << "(" << operater[0] << "," << operater[1] << ")    " << Js[k] <<  endl;

}

double apply_operator(int op0, int op1,int chain[][2], int JJ1)
{
  int index1[2] = {0};
  int index2[2] = {0};

  //find indices of the sites in question
  for (int i4 = 0; i4 < half_L; i4++)
    {
      for (int jb2 = 0; jb2 < 2; jb2++)
	{
	  if(op0 == chain[i4][jb2])
	    {
	      index2[0] = i4;
	      index2[1] = jb2;
	    }

	  if(op1 == chain[i4][jb2])
	    {
	      index1[0] = i4;
	      index1[1] = jb2;
	    }
	}
    }

  if(index1[0] == index2[0])
    { 
      if(JJ1 == 0){x_new[0] +=1;}
      else{x_new[1] += 1;}
    }

  //applying operator and changing chain
  else
    {
      chain[index1[0]][index1[1]] = chain[index2[0]][(index2[1]+1)%2];   
      chain[index2[0]][(index2[1]+1)%2] = op1;
      
      //  print_chain(chain, half_L);
      if(JJ1 == 0){w_new[0] += 1;}
      else{w_new[1] += 1;}
    }
  
  return 0;

}
void change_operators(int operaters[][2], int Js[],int a, int neighbours[][2])
{

  int first = irand() %a;
  int second = irand() %a;
  // int third = irand() %a;
  //  int fourth = irand() %a; 

  while(first == second){second = irand() %a;}
  //  while(first == third) {third = irand() %a;} while(third == second){third = irand() %a;}
  //  while(fourth == first){fourth = irand() %a;} while(fourth == second){fourth = irand() %a;} while(fourth == third){fourth = irand() %a;}
  
  int changings[2] = {first,second};//,third};//,fourth};

  for(int m=0; m<2; m++)
    {
      int news[2] = {0};
      
      Js[changings[m]] = 0;
      int initb = (irand()+L) %L;
      int neighb = (irand()+2)%2;

      //OPEN BCs -------------------------------------------------------------------------------------
      while((initb == 0 && neighb == 0) || (initb == L-1 && neighb == 1))
	{
	  initb = (irand()+L) %L;
	  neighb = (irand()+2)%2;
	}
      //END OF OPEN BCs -------------------------------------------------------------------------------*/
      
      news[0] = initb;
      news[1] = neighbours[initb][neighb];
      
      operaters[changings[m]][0] = news[0];
      operaters[changings[m]][1] = news[1];


      if( news[1] > L)
	{
	  cout<<"ERROR"<<endl;
	  cout << initb << "," << neighb << endl;
	  exit(1);
	}
    
    
      
      if((news[0]+news[1])%2 == 1){Js[changings[m]] = 1;}
    } 
}
