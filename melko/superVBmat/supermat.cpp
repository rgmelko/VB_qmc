// May 7, 2010 // Making supermatrix to show overcompleteness

#include"header.h"
#include"lapack.h"
#include<vector>

int loop_counter(vector <int> leftside, vector <int> rightside);

int main(){  

  int supermat[24][24] = {0};
  Array<double , 2> doublemat(24,24);

  int states[8][24] = {0};

  states[0][0] = 1;
  states[2][0] = 3;
  states[5][0] = 4;
  states[7][0] = 6;

  states[0][1] = 1;
  states[2][1] = 3;
  states[5][1] = 6;
  states[7][1] = 4;

  states[0][2] = 1;
  states[2][2] = 4;
  states[5][2] = 6;
  states[7][2] = 3;

  states[0][3] = 1;
  states[2][3] = 4;
  states[5][3] = 3;
  states[7][3] = 6;

  states[0][4] = 1;
  states[2][4] = 6;
  states[5][4] = 3;
  states[7][4] = 4;

  states[0][5] = 1;
  states[2][5] = 6;
  states[5][5] = 4;
  states[7][5] = 3;

  states[0][6] = 3;
  states[2][6] = 1;
  states[5][6] = 4;
  states[7][6] = 6;

  states[0][7] = 3;
  states[2][7] = 1;
  states[5][7] = 6;
  states[7][7] = 4;

  states[0][8] = 3;
  states[2][8] = 4;
  states[5][8] = 1;
  states[7][8] = 6;

  states[0][9] = 3;
  states[2][9] = 4;
  states[5][9] = 6;
  states[7][9] = 1;

  states[0][10] = 3;
  states[2][10] = 6;
  states[5][10] = 1;
  states[7][10] = 4;

  states[0][11] = 3;
  states[2][11] = 6;
  states[5][11] = 4;
  states[7][11] = 1;

  states[0][12] = 4;
  states[2][12] = 1;
  states[5][12] = 3;
  states[7][12] = 6;

  states[0][13] = 4;
  states[2][13] = 1;
  states[5][13] = 6;
  states[7][13] = 3;

  states[0][14] = 4;
  states[2][14] = 3;
  states[5][14] = 1;
  states[7][14] = 6;

  states[0][15] = 4;
  states[2][15] = 3;
  states[5][15] = 6;
  states[7][15] = 1;

  states[0][16] = 4;
  states[2][16] = 6;
  states[5][16] = 1;
  states[7][16] = 3;

  states[0][17] = 4;
  states[2][17] = 6;
  states[5][17] = 3;
  states[7][17] = 1;

  states[0][18] = 6;
  states[2][18] = 1;
  states[5][18] = 3;
  states[7][18] = 4;

  states[0][19] = 6;
  states[2][19] = 1;
  states[5][19] = 4;
  states[7][19] = 3;

  states[0][20] = 6;
  states[2][20] = 3;
  states[5][20] = 1;
  states[7][20] = 4;

  states[0][21] = 6;
  states[2][21] = 3;
  states[5][21] = 4;
  states[7][21] = 1;

  states[0][22] = 6;
  states[2][22] = 4;
  states[5][22] = 1;
  states[7][22] = 3;

  states[0][23] = 6;
  states[2][23] = 4;
  states[5][23] = 3;
  states[7][23] = 1;

  for(int i=0; i<24; i++){

      states[states[0][i]][i] = 0;
      states[states[2][i]][i] = 2;
      states[states[5][i]][i] = 5;
      states[states[7][i]][i] = 7;

  }

  vector <int> leftvect(8), rightvect(8);

  for(int i=0; i<24; i++){

    for(int k=0; k<8; k++){
      leftvect[k] = states[k][i];
    }
    
    for(int j=0; j<24; j++){
      
      for(int l=0; l<8; l++){
	rightvect[l] = states[l][j];
      }

      supermat[i][j] = loop_counter(leftvect, rightvect);

    }
  }


    for(int i=0; i<24; i++){
      for(int j=0; j<24; j++){
	doublemat(i,j) = pow(2.0,supermat[i][j]-4);
      }
    }



  cout.flush();
  vector<double> dd;
  cout << "starting diagonalization\n";
  diagWithLapack_R(doublemat,dd);  //*** LAPACK DIAG
			
  cout<<"Done diagonalizting the Hamiltonian \n";
  cout.flush();
  
  for (int ii=0; ii<dd.size(); ii++) {
   cout<<setprecision(12)<<dd.at(ii)<<"\n";
  }
  cout<<endl;
  
  return 0;
}


int loop_counter(vector <int> leftside, vector <int> rightside)
{
  int counter(0), loopnumber(0), startsite(0), AA(-99), which(0);
  vector <int> site(10,0);

  while(counter < 8){

    site[counter]=1;
    startsite = counter;
    

    AA = leftside[counter];
    which = 0;
   
    while(AA!=startsite){
      site[AA] = 1;
      
      if(which==0){
	AA = rightside[AA];
	which++;
      }
      else{
	AA = leftside[AA];
	which--;
      }
      site[AA] = 1;
    }
    loopnumber++;
    while(site[counter]==1){counter++;}
  }
  return loopnumber;
}
