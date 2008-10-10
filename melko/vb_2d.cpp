// vb_2d.cpp trying to add a dimension. Last updated May 22, 2008
// Comment by R Melko

#include<iostream>
#include<math.h>
#include"mtrand.h" // ramdom number generator
using namespace std;

void shuffle(int [][2]); //randomizes initial bonds (currently not working)
void print_chain(int chain [][2]); //prints the bonds
void generate_operator(int operater[2], int neighbours[][4]);//generates 1 bond 
                                                                   //operator
int apply_operator(int op0, int op1, int chain[][2]);//applies 1 bond operator
void change_operators(int operaters[][2], int a, int neighbours[][4]);/*randomly 
changes a number of bond operators... where that number is "a"*/

MTRand drand; //drand() gives you a random double precision number
MTRand_int32 irand; // irand() gives you a random integer

const int L = 6; // 1-D length of the lattice
const int zone = 2; // the size of "the zone"
const int L2 = L*L; // total number of sites
const int half_L = L2/2; // total number of sites divided by 2
const int n = L2*6; // number of bond operators
const int start = 310000; /* number of iterations until the programs takes 
			     measurements  */
const int iterations = start*2; // total number of iterations
const int array_sz = iterations;
int chain [half_L][2] = {0}; // the bonds are stored in here
int operater[2] = {0}; //it's an operator
int initial_state[half_L][2] ={0}; /*stores the initial bonds which I'm using
				     instead of randomizing them*/

int neighbours[L2][4]; //lists the 4 nearest neighbours for each site
int box[L2] = {0}; // the zone! <-- exclamation point

main() // the main program..
{

  irand();

  cout << "L = " << L << "    " << "zone = " << zone << "     " << iterations << " iterations"<< endl;
  cout.precision(10); // ten digits of precision..  or ten decimal places?

  //******Finding Nearest Neighbours********************************
  for(int iii=0; iii<L2; iii++)
    {
      if(iii%L==0){neighbours[iii][0]=iii+L-1;}
      else neighbours[iii][0]=iii-1;
      neighbours[iii][1]=(iii+L)%L2;
      if(iii%L==L-1){neighbours[iii][2]=iii-L+1;}
      else neighbours[iii][2]=iii+1;
      neighbours[iii][3]=(iii-L+L2)%L2;
      
      //    cout << neighbours[iii][0] << ", ";
      //    cout << neighbours[iii][1] << ", ";
      //    cout << neighbours[iii][2] << ", ";
      //    cout << neighbours[iii][3] << ", ";
      //    cout << endl;
    }



  //****Define the zone-box-type-thing*******************************
  for(int abox=0; abox<L*zone; abox+=L)
    {
      box[abox]=1;
      // cout << box[abox] << ", ";
      for(int bbox=1; bbox<zone-1; bbox++)
	{
	  box[abox+bbox]=1;
	  box[neighbours[abox+bbox][2]]=1;
	  //  cout << box[abox+bbox] << ", " ;
	}
      //  cout << endl;
    }

  /*****Print out the zone box*****************************************
        for (int zzz = 0; zzz<L2; zzz+=L)   
           {
             for(int xxx = 0; xxx<L; xxx++)
       	       {
       	          cout << box[xxx+zzz] << ", ";
       	       }
             cout << endl;
           }
  *************************************************************************/

  //  cout << endl;

  shuffle(initial_state);//"randomize" the initial state (but not really)
 
  print_chain(initial_state);

  FILE * bondss;                       // open a file to print the number of
  bondss = fopen("bonds.txt", "w");  // nearest neighbour bonds in

  int operaters[n][2], new_operaters[n][2];   // old and new operators
  int bonds[array_sz] = {0};                 /* number of nn bonds 
						  for each iteration*/
  int acc = 0, rej =0;           //number of changes accepted and rejected
  long int w_old=0, w_new=0;   // the new and old weights
  long int superweights[array_sz] = {0}; // storing alllllll the weights
  int cross[array_sz] = {0}; //the number of bonds crossing the zone boundary
  double energy = 0;          // the energy
  double entropy = 0;

  //-------Generate Operators----------------------------------------
  for(int i0=0; i0<n; i0++)
    {  
      generate_operator(operater, neighbours);
      operaters[i0][0] = operater[0];
      operaters[i0][1] = operater[1];
    }
  //-----------------------------------------------------------------

  //--------Initialize Bonds------------------------------------------
  for(int j=0; j<half_L; j++)
    {
      chain[j][0] = initial_state[j][0];
      chain[j][1] = initial_state[j][1];
    }
  //------------------------------------------------------------------

  //-------Apply Operators-------(also get the weight)------------------
  for(int k=0; k<n; k++)
    {
      w_old += apply_operator(operaters[k][0],operaters[k][1],chain); 
    }
  //-------------------------------------------------------------------
  int q = 0; //the current index for bonds which records the number of nn bonds

  for (int i=0; i<iterations; i++)
    {  
      if(i%100000==0){cout<<i<<" steps"<<endl;}

      for(int i7=0; i7<n; i7++)
	{ 
	  new_operaters[i7][0] = operaters[i7][0];
	  new_operaters[i7][1] = operaters[i7][1];
	}

      change_operators(new_operaters, n, neighbours);

       for(int j=0; j<half_L; j++)
	{
	  chain[j][0] = initial_state[j][0];
	  chain[j][1] = initial_state[j][1];
	}

      for(int k=0; k<n; k++)
	{
	  w_new += apply_operator(new_operaters[k][0],new_operaters[k][1],chain);
	}

      cout << "prob = " << pow(2,w_old-w_new) << endl;

      if(drand() < pow(2,w_old-w_new))
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
		    {bonds[q]+=1;}
		  if(box[chain[id][0]]+box[chain[id][1]]==1){cross[q]+=1;}
		}
		  superweights[q]=w_new;
	      q+=1;
	      acc+=1;
	      
	    }

	  for(int i8=0; i8<n; i8++)
	    {
	      operaters[i8][0] = new_operaters[i8][0];
	      operaters[i8][1] = new_operaters[i8][1];
	    }
	  w_old = w_new; 
	}

      else{if(i>=start)
	  {
	    if(q==0){} 
	    else
	      {
		bonds[q]=bonds[q-1];
		superweights[q]=superweights[q-1];
		cross[q]=cross[q-1]; 
		q+=1;
		rej+=1;
	      }
	  }
      }
      
      w_new = 0;
    }

  double superduperweights = 0;
  
//   cout<< superweights[0] << ","<< superweights[1]<< "," << superweights[2]<<"," << superweights[3] <<"," ;
//   cout<< superweights[4] << ","<< superweights[5]<< "," << superweights[6]<<"," << superweights[7]<<"," ;
//   cout<< superweights[8] << ","<< superweights[9]<< "," << superweights[10]<<"," << superweights[11]<< endl;

  for(int i9=0; i9<q; i9++)
    {
      //   cout << "bonds[" << i9 << "] = " << bonds[i9] << endl;
      energy += bonds[i9];
      entropy += cross[i9]*pow(0.5,superweights[i9]);
      superduperweights += pow(0.5,superweights[i9]);
    }

  //  cout << endl << energy << " bonds" << endl;

  energy /= -L2*(q);
  energy += -1;
  energy *= 0.5;
  //  energy -= L2*(0.25);
  entropy /= superduperweights;
  entropy *= log(2);

 
  cout << "energy = " << energy << endl;
  cout << "entropy = " << entropy << endl;
  cout << "entropy/x = " << entropy/zone << endl;
 	   
  //  cout << q << " energies" << endl;

  cout << 100*double(acc)/q << "% accepted"<< endl;

  for(int j=0; j<half_L; j++)
    {
      chain[j][0] = initial_state[j][0];
      chain[j][1] = initial_state[j][1];
    }

  for(int k=0; k<n; k++)
    {
      w_new += apply_operator(operaters[k][0],operaters[k][1],chain);
    }
  
  //  print_chain(chain, half_L);

  int z = 0;
  for(int i=0; i<q; i+=16)
    {
      bonds[z] = bonds[i]+bonds[i+1]+bonds[i+2]+bonds[i+3]+bonds[i+4]+bonds[i+5]+bonds[i+6]+bonds[i+7]
  	+bonds[i+8]+bonds[i+9]+bonds[i+10]+bonds[i+11]+bonds[i+12]+bonds[i+13]+bonds[i+14]+bonds[i+15];
      z +=1;
    }
  
  for(int i=0; i<z; i++)
    {
      fprintf(bondss, "%i", bonds[i]);
      fprintf(bondss, ",\n ");
    }  

  fclose(bondss);

  //  cout << endl;

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

  for (int site = 0; site < L2;)
    {
      rowcol =  (irand() + 2) %2;
      if(rowcol==0)
	{
	  row =  (irand()+half_L) %half_L;
	  while((chain[row][0] != -1)&(rowcol==0))
	    { 
	      row = (irand()+half_L) %half_L; 
	      check += 1;
	      if(check==half_L*half_L*half_L){rowcol=1;}
	    }
	  check = 0;
	  if(rowcol==0){chain[row][0]=site; site++;}
	}
      if(rowcol==1)
	{
	  col =  (irand()+half_L) %half_L;
      
	  while((chain[col][1] != -1)&(rowcol==1))
	    { 
	      col = (irand()+half_L) %half_L; 
	      check += 1;
	      if(check==half_L*half_L*half_L){rowcol=0;}
	    }
	  check = 0;
	  if (rowcol==1){chain[col][1]=site;site++;}
	}


//         chain[row][0] = site;
//         chain[row][1] = site+1;
//         row += 1;
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

void generate_operator(int operater[2], int neighbours[][4])
{
  operater[0] = irand() %L2;
  operater[1] = neighbours[operater[0]][irand()%4];
  
  //cout << '(' << operater[0] << ',' << operater[1] << ") ";       
}

int apply_operator(int op0, int op1, int chain[][2])
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

  if(index1[0] == index2[0]) {// print_chain(chain);
  return 0;}

  //applying operator and changing chain
  chain[index1[0]][index1[1]] = chain[index2[0]][(index2[1]+1)%2];   
  chain[index2[0]][(index2[1]+1)%2] = op1;

  //  print_chain(chain, half_L);

  return 1;

}
void change_operators(int operaters[][2], int a, int neighbours[][4])
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

      generate_operator(news, neighbours);

      while((operaters[changings[m]][0] == news[0]) && (operaters[changings[m]][1] == news[1]))
	{
	  generate_operator(news, neighbours);
	}

      operaters[changings[m]][0] = news[0];
      operaters[changings[m]][1] = news[1];
    } 
}
