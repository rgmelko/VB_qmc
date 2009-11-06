// File for averaging VB EE data
//
// February 9, 2009
// Roger Melko
//--------------------------------------------------------
// Changin' it up Feb 19, 2009
//
// I think it's actually working. It calculates the average
// entropies and their errors (sdom)
//--------------------------------------------------------
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
using namespace std;

int main(int argc, char *argv[]){
  
  if (argc < 4) {
    cout<<"I expect a command line argument: filename, number of columns, number to exclude.\n";
    exit(1);
  }
  
  double temp;
    
  //Number of spins
  string filename;
  int nSpin;
  int nExclude;
  filename = argv[1];
  nSpin = atoi(argv[2]); //first command line parameter
  nExclude = atoi(argv[3]);
  
  ifstream cfin(filename.c_str());
  // cout << filename << endl << nSpin << endl << nExclude << endl;
  
  vector<double> EA;
  vector<double> nontemp;
  EA.assign(nSpin,0);
  nontemp.assign(nSpin,0);
  

  int index = 0;
  int numRow = 0;

  double dummy;
  for(int i=0; i<nExclude; i++)
    {
      for(int j=0; j<nSpin; j++){cfin>>dummy;}
    }
  
  cfin>>nontemp[index];
 
  while( !cfin.eof() ){
    numRow++;
    nontemp.resize(nSpin*numRow+1000);
    
    EA.at(0) += nontemp[index];
    index++;
    
    for (int i=1; i<EA.size(); i++){
      cfin>>nontemp[index];
      EA.at(i) += nontemp[index];
      index++;
    }
    
    cfin>>nontemp[index];
  };
  // cout<<numRow<<endl;

  
  vector<double> summ(nSpin);
  summ.assign(nSpin,0);


  for(int i=0; i<nSpin; i++)
    {
      for(int j=0; j<numRow; j++)
	{
	  index = j*nSpin+i;
	  nontemp[index] -= EA.at(i)/numRow;
	  summ[i] += nontemp[index]*nontemp[index];
	}
      summ[i] = sqrt(summ[i]/(numRow*(numRow-1)));
      
    }

  
  //average
  for (int i=0; i<EA.size(); i++){
    EA.at(i) /= numRow;
    cout<<i+1<<" "<<EA.at(i)<<" "<<summ[i]<<endl;
  }
  
  cfin.close();

}
