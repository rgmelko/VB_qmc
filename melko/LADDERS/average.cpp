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
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h>
using namespace std;

int main(int argc, char *argv[]){

  cout.precision(10);
  
  if (argc < 4) {
    cout<<"I expect a command line argument: filename, number of columns, number to exclude.\n";
    exit(1);
  }
    
  //Number of spins
  string filename;
  int nSpin;
  int nExclude;
  filename = argv[1];
  nSpin = atoi(argv[2]); //first command line parameter
  nExclude = atoi(argv[3]);//second command line parameter
  
  ifstream cfin(filename.c_str());
  // cout << filename << endl << nSpin << endl << nExclude << endl;
  
  vector<double> Sum, Sum_Sq;
  double temp(0);
  Sum.assign(nSpin,0);
  Sum_Sq = Sum;

  int numRow = 0;

  double dummy;
	//throw out the excluded data
	for(int i=0; i<nExclude; i++) //number of rows to exclude
    {
      for(int j=0; j<nSpin; j++){cfin>>dummy;} //number of columns per row
    }
  
  cfin>>temp;//read in the first entry

		
  while( !cfin.eof() ){
    numRow++; //increase the number of rows each time
    
    Sum.at(0) += temp; //add the first value to sum
	Sum_Sq.at(0) += temp*temp;
    
    for (int i=1; i<Sum.size(); i++){ //for the rest of the row
      cfin>>temp;          //read in the rest of the values
      Sum.at(i) += temp;    //add them to sum
	  Sum_Sq.at(i) += temp*temp;
	  
    }
    
    cfin>>temp;            //read in the first value of the next row
  };
  // cout<<numRow<<endl;

  
  //average
  for (int i=0; i<Sum.size(); i++){
    Sum.at(i) /= numRow;
	Sum_Sq.at(i)/= numRow;
    cout<<i+1<<" "<<Sum.at(i)<<" "<<sqrt((Sum_Sq[i]-Sum[i]*Sum[i])/((numRow-1))) << endl;
  }
  
  cfin.close();

}
