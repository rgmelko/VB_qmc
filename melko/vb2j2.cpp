// vb2j.cpp  Last updated Oct 23, 2008
// Trying to add J''

#include<iostream>
#include<math.h>
#include"mtrand.h" // ramdom number generator
using namespace std;

void shuffle(int [][2]); //randomizes initial bonds (currently not working)
void print_chain(int chain [][2]); //prints the bonds
void generate_operator(int operater[2], int neighbours[][4], int Js[], int indexx);
//generates 1 bond operator
double apply_operator(int op0, int op1, int chain[][2], int Jayrock);
//applies 1 bond operator
void change_operators(int operaters[][2], int Js[], int a, int neighbours[][4]);/*randomly 
changes a number of bond operators... where that number is "a"*/

MTRand drand; //drand() gives you a random double precision number
MTRand_int32 irand; // irand() gives you a random integer

const long int superseed = 827193545; // ********You************************
const int L = 4; // 1-D length of the lattice *******Can********************
const int zone = 1; // the size of "the zone" *********Change***************
const double jprime =2; // ****************************These*Values*******
double J = 1;
const int L2 = L*L; // total number of sites
const int half_L = L2/2; // total number of sites divided by 2
const int n = L2*4; // number of bond operators
const int start = 100000; /* number of iterations until the programs takes ***
			     measurements  */
const int iterations = start*10; // total number of iterations
int chain [half_L][2] = {0}; // the bonds are stored in here
int operater[2] = {0}; //it's an operator
int initial_state[half_L][2] ={0}; //stores the initial bond configuration

int Js[n] = {0};   // Stores the interaction strength for each operator 0=J,1=J'
double w_old[2]={0};  // the new and old weights
double w_new[2]={0};  // the new and old weights
double x_old[2]={0};  // the new and old weights
double x_new[2]={0};  // the new and old weights
double probb = 0; // prob of keeping new operators

int neighbours[L2][4]; //lists the 4 nearest neighbours for each site
int box[L2] = {0}; // the zone! <-- exclamation point

main() // the main program..
{
  irand.seed(superseed);

  cout << "L = " << L << "    " << "zone = " << zone << "   " << 
    iterations << " iterations" << "    n = " << n;
  cout << "     J = " << J << "   J' = " << jprime << endl;
  cout.precision(10); // ten digits of precision..  or ten decimal places?

  //******Finding Nearest Neighbours********************************
  for(int iii=0; iii<L2; iii++)
    {
      if(iii%L==0){neighbours[iii][0]=iii+L-1;}   //North
      else neighbours[iii][0]=iii-1;
      neighbours[iii][1]=(iii+L)%L2;                //East
      if(iii%L==L-1){neighbours[iii][2]=iii-L+1;} //South
      else neighbours[iii][2]=iii+1;
      neighbours[iii][3]=(iii-L+L2)%L2;            // West
      
      //    cout << neighbours[iii][0] << ", ";
      //    cout << neighbours[iii][1] << ", ";
      //    cout << neighbours[iii][2] << ", ";
      //    cout << neighbours[iii][3] << ", ";
      //    cout << endl;
    }

  //****Define the zone-box-type-thing*******************************
      for(int abox=0; abox<L*zone; abox+=L)
      {
	for(int bbox=0; bbox<zone; bbox++)
	  {
	    box[abox+bbox]=1;
	    //  cout << box[abox+bbox] << ", " ;
	  }
	//  cout << endl;
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

  shuffle(initial_state);//"randomize" the initial state (but not really)
 
  //  print_chain(initial_state);

  int operaters[n][2], new_operaters[n][2];   // old and new operators     
  int bond[2] = {0};   //number of NN J bonds
  int bondprime[2] = {0};  // number of NN J' bonds
  int acc = 0, rej =0;         //number of changes accepted and rejected
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

  for (int i=0; i<iterations; i++)
    {  
      if(i%start==0){cout<< q2 << "% "<< endl ;  q2+=10;}

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
	}


//       cout << endl;
//       print_chain(chain);

      probb =  pow(0.5,(w_new[0]+w_new[1]-w_old[0]-w_old[1]))*
	pow(J,(w_new[0]+x_new[0]-w_old[0]-x_old[0]))*
	pow(jprime,(w_new[1]+x_new[1]-w_old[1]-x_old[1]));

  
//     cout << probb << endl << 
//  	"wold0 = " << w_old[0] << "   wnew0 = " << w_new[0] << endl
//  	   <<"wold1 = " << w_old[1] << "   wnew1 = " << w_new[1] << endl
//  	   <<"xold0 = " << x_old[0] << "   xnew0 = " << x_new[0] << endl
//  	   <<"xold1 = " << x_old[1] << "   xnew1 = " << x_new[1] << endl;
	  

	if(drand() < probb)
	{
	  if(i >= start)
	    {
	      for(int id=0; id<half_L; id++)
		{
		  if(
		     (chain[id][1] == neighbours[chain[id][0]][0])|
		     (chain[id][1] == neighbours[chain[id][0]][1])|
		     (chain[id][1] == neighbours[chain[id][0]][2])|
		     (chain[id][1] == neighbours[chain[id][0]][3])
		     )
		    {
		      if(
			 ((chain[id][1] == neighbours[chain[id][0]][0])&
			  (((chain[id][0]/L)%2 + chain[id][0]%2)%2 == 1))|
			 ((chain[id][1] == neighbours[chain[id][0]][2])&
			  (((chain[id][0]/L)%2 + chain[id][0]%2)%2 == 0))
			 )
			{bondprime[1]+=1;}
		      else{bond[1]+=1;}
		    }
		
		  if(box[chain[id][0]]+box[chain[id][1]]==1){cross[1]+=1;}
		}
	      energy += bond[1];
	      energyprime += bondprime[1];
	      entropy += cross[1];
	      q+=1;
	      acc+=1;
	      
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
		energy += bond[0];
		energyprime += bondprime[0];
		entropy += cross[0];
		q+=1;
		rej+=1;
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

  cout << "bonds " <<  energy << "     bonds' "<< energyprime << endl;



  energy /= (2*L2*0.75*q);
  energy += 1 - 0.5;
  energy *= -0.5*2*L2*0.75*J; 
  energyprime /= (2*L2*0.25*q);
  energyprime += 1 - 0.5;
  energyprime *= -0.5*2*L2*0.25*jprime; 
  energy = (energy + energyprime);
  
  //  energy *= L2;

  entropy /= q;
  entropy *= log(2);
 
  cout << endl;
  cout << "energy = " << energy << endl;
  cout << "entropy = " << entropy << endl;
  cout << "entropy/x = " << entropy/zone << endl;
 	   
  cout << q << " energies" << endl;

  cout << 100*double(acc)/q << "% accepted"<< endl;

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

void generate_operator(int operater[3], int neighbours[][4], int Js[], int k)
{
  Js[k] = 0;
  int initb = (irand()+L2) %L2;
  operater[0] = initb;
  int neighb = (irand()+4)%4;
  operater[1] = neighbours[operater[0]][neighb];
  if(
     ((((initb/L)%2+initb%2)%2==1)&(neighb==0))
     |((((initb/L)%2+initb%2)%2==0)&(neighb==2)))
    {Js[k]=1;}
      
}

double apply_operator(int op0, int op1,int chain[][2], int JJ1)
{
  int index1[2] = {0};
  int index2[2] = {0};

  //find indices of the sites in question
  for (int i4 = 0; i4 < half_L; i4++)
    {
      for (int j2 = 0; j2 < 2; j2++)
	{
	  if(op0 == chain[i4][j2])
	    {
	      index2[0] = i4;
	      index2[1] = j2;
	    }

	  if(op1 == chain[i4][j2])
	    {
	      index1[0] = i4;
	      index1[1] = j2;
	    }
	}
    }

  if(index1[0] == index2[0])
    { if(JJ1 == 0){x_new[0] +=1; return 0;}
      x_new[1] += 1;
      return 0;
    }

  //applying operator and changing chain
  chain[index1[0]][index1[1]] = chain[index2[0]][(index2[1]+1)%2];   
  chain[index2[0]][(index2[1]+1)%2] = op1;

  //  print_chain(chain, half_L);
  if(JJ1 == 0){w_new[0] += 1; return 0;}
  w_new[1] += 1; 

  return 0;

}
void change_operators(int operaters[][2], int Js[],int a, int neighbours[][4])
{

  int first = irand() %a;
  int second = irand() %a;
  int third = irand() %a;
  //  int fourth = irand() %a; 

  while(first == second){second = irand() %a;}
  while(first == third) {third = irand() %a;} while(third == second){third = irand() %a;}
  //  while(fourth == first){fourth = irand() %a;} while(fourth == second){fourth = irand() %a;} while(fourth == third){fourth = irand() %a;}
  
  int changings[3] = {first,second,third};//,fourth};

  for(int m=0; m<3; m++)
    {
      int news[2] = {0};

      generate_operator(news, neighbours, Js, changings[m]);


      while((operaters[changings[m]][0] == news[0]) && (operaters[changings[m]][1] == news[1]))
	{
	  generate_operator(news, neighbours, Js, changings[m]);
	}

      operaters[changings[m]][0] = news[0];
      operaters[changings[m]][1] = news[1];
    } 
}
