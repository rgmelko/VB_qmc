// vb.cpp attempt number 2.  Last updated Wed Mar 26, 2008

#include<iostream>
#include<math.h>
#include"mtrand.h"
using namespace std;

void shuffle(int [][2], int a);
void print_chain(int chain [][2], int a);
void generate_operator(int operater[2]);
int apply_operator(int op0, int op1, int chain[][2], int n);
void change_operators(int operaters[][2], int a);

MTRand drand;
MTRand_int32 irand;

const int L = 10;
const int half_L = L/2;
const int n = L*6;
const int start = 100000;
const int iterations[5] = {start+20000, start+40000, start+80000, start+160000, start+320000};
int chain [half_L][2] = {0};
int operaters[n][2], new_operaters[n][2] = {0};
int operater[2] = {0};
int initial_state[half_L][2] ={0};

main()
{
  cout.precision(10);
  cout << endl;

  shuffle(initial_state, half_L);
 
  print_chain(initial_state, half_L);

  for(int i0=0; i0<n; i0++) // generate operators
    {  
      generate_operator(operater);
      operaters[i0][0] = operater[0];
      operaters[i0][1] = operater[1];
    }


 FILE * bondss;
  bondss = fopen("bonds.txt", "w");

  for(int super_its=0; super_its<5; super_its++)
    {

      drand();irand();
      if(super_its == 0)bondss = fopen("bonds_1.txt", "w");
      if(super_its == 1)bondss = fopen("bonds_2.txt", "w");
      if(super_its == 2)bondss = fopen("bonds_3.txt", "w");
      if(super_its == 3)bondss = fopen("bonds_4.txt", "w");
      if(super_its == 4)bondss = fopen("bonds_5.txt", "w");
      if(super_its == 5)bondss = fopen("bonds_6.txt", "w");
      if(super_its == 6)bondss = fopen("bonds_7.txt", "w");
      if(super_its == 7)bondss = fopen("bonds_8.txt", "w");
      if(super_its == 8)bondss = fopen("bonds_9.txt", "w");
      if(super_its == 9)bondss = fopen("bonds_10.txt", "w");
      if(super_its == 10)bondss = fopen("bonds_11.txt", "w");
      if(super_its == 11)bondss = fopen("bonds_12.txt", "w");
      if(super_its == 12)bondss = fopen("bonds_13.txt", "w");

      int bonds[90000000] = {0};
      int acc = 0, rej =0; 
      long int w_old=0, w_new=0;
      double energy = {0};

    
  for(int j=0; j<half_L; j++)
    {
      chain[j][0] = initial_state[j][0];
      chain[j][1] = initial_state[j][1];
    }
  
  for(int k=0; k<n; k++)
    {
      w_old  += apply_operator(operaters[k][0],operaters[k][1],chain, half_L); 
    }
  
  int q = 0;

  for (int i=0; i<iterations[super_its]; i++)
    {  
      if(i%100000==0){cout<<i<<endl;}

      for(int i7=0; i7<n; i7++)
	{ 
	  new_operaters[i7][0] = operaters[i7][0];
	  new_operaters[i7][1] = operaters[i7][1];
	}

      change_operators(new_operaters, n);

      for(int j=0; j<half_L; j++)
	{
	  chain[j][0] = initial_state[j][0];
	  chain[j][1] = initial_state[j][1];
	}

      for(int k=0; k<n; k++)
	{
	  w_new += apply_operator(new_operaters[k][0],new_operaters[k][1],chain, half_L);
	}

      if(drand() < pow(2.0,(w_old-w_new)))
	{
	  if(i >= start)
   	    {
		  for(int id=0; id<half_L; id++)
		    {
		      if((chain[id][0] == (chain[id][1]+L+1)%L) | (chain[id][0] == (chain[id][1]+L-1)%L)){bonds[q]+=1;}
		    }
	      q+=1;
	      acc+=1;
	    }

	  //  cout <<"hep" <<endl;
	  for(int i8=0; i8<n; i8++)
	    {
	      operaters[i8][0] = new_operaters[i8][0];
	      operaters[i8][1] = new_operaters[i8][1];
	    }
	  w_old = w_new; 
	}
      else{if(i>=start){if(q==0){} else{bonds[q]=bonds[q-1]; q+=1; rej+=1;}}}
      
      w_new = 0;
    }

  for(int i9=0; i9<q-1; i9++){energy += bonds[i9];}
  
  energy /= -2*(q-1);
  energy -= L*(.25);
  energy /= L;
 
  cout << "energy = " << energy << endl;
 	   
  cout << q << " energies" << endl;

  cout << 100*double(acc)/q << "% accepted"<< endl;

  for(int j=0; j<half_L; j++)
    {
      chain[j][0] = initial_state[j][0];
      chain[j][1] = initial_state[j][1];
    }

  for(int k=0; k<n; k++)
    {
      w_new += apply_operator(operaters[k][0],operaters[k][1],chain, half_L);
    }
  
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

    }

  cout << endl;

  return 0;
}

void shuffle(int chain[][2], int a)
{

  for(int i=0; i<a; i++){chain[i][0] = -1;}

  int row, column;

  for (int site = 0; site < L; site+=2)
    {
      row = irand() %a;
      
      while(chain[row][0] != -1)
	{
	  row = irand() %a;
    	}

      chain[row][0] = site;
      chain[row][1] = site+1;
    }
}

void print_chain(int chain[][2], int a)
{
  for (int i2 = 0; i2 < a; i2++)
    {
      cout << chain[i2][0] << ',' << chain[i2][1] << endl;
    }

  cout << endl;
}

void generate_operator(int operater[2])
{
  operater[0] = irand() %L;
  operater[1] = (operater[0] + ((irand() %2)*2 - 1) + L)%L;
  
  // cout << '(' << operater[0] << ',' << operater[1] << ") ";         
}

int apply_operator(int op0, int op1, int chain[][2], int a)
{
  int index1[2] = {0};
  int index2[2] = {0};

  //find indices of the sites in question
  for (int i4 = 0; i4 < a; i4++)
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

  if(index1[0] == index2[0]) {// print_chain(chain, half_L);
  return 0;}

  //applying operator and changing chain
  chain[index1[0]][index1[1]] = chain[index2[0]][(index2[1]+1)%2];   
  chain[index2[0]][(index2[1]+1)%2] = op1;

  //  print_chain(chain, half_L);

  return 1;
}

void change_operators(int operaters[][2], int a)
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

      generate_operator(news);

      while((operaters[changings[m]][0] == news[0]) && (operaters[changings[m]][1] == news[1]))
	{
	  generate_operator(news);
	}

      operaters[changings[m]][0] = news[0];
      operaters[changings[m]][1] = news[1];
    } 
}
