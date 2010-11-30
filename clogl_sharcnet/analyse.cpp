#include<iostream>
#include<fstream>
#include<iomanip>
#include<sstream>
#include<vector>
#include<cmath>
using namespace std;

void average(string fname, int rows, string syze, int exclude);

int main(){

  int flipnum = 1;
  const int num_of_sizes = 2;
  int sys_size[num_of_sizes] = {20, 30};
  int size;
  string flipstring, sizestring, filename;
  stringstream out;
  int nExclude = 0; //number of rows of data to exclude;

  //loop through all the system sizes
  for(int j=0; j<num_of_sizes; j++)
    {
      size = sys_size[j];
      out.str(""); //clear out
      out << size; //put size in out
      sizestring = out.str(); //make out a string
      
      //loop through the swap ratios.. the denominators..
      for(flipnum=1; flipnum<size-1; flipnum++)
	{
	  //make the swap number into a string (for the filename)
	  out.str("");
	  if(flipnum < 10){ out << "0"; } //make sure it's two digits
	  out << flipnum;
	  flipstring = out.str();

	  //these are the filenames.. now i can read them in
	  filename = "entrpy"+sizestring+"swap"+flipstring+".dat";
	  
	  average(filename, size, sizestring, nExclude);
	}
  
  
    }



  return 0;
}


void average(string fname, int rows, string syze, int exclude){
  
  ifstream cfin(fname.c_str());
  
  vector<double> Sum, Sum_Sq;
  double temp(0);
  Sum.assign(rows,0);
  Sum_Sq = Sum;

  int numRow = 0;

  double dummy;
  //throw out the excluded data
  for(int i=0; i<exclude; i++) //number of rows to exclude
    {
      for(int j=0; j<rows; j++){cfin>>dummy;} //number of columns per row
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
 

  string outfile = "./"+syze+"/av_"+fname;
  ofstream fout(outfile.c_str());

  //average
  for (int i=0; i<Sum.size(); i++){
    Sum.at(i) /= numRow;
    Sum_Sq.at(i)/= numRow;
    fout<<i+1<<" "<<-log(Sum.at(i))<<" "<<(sqrt((Sum_Sq[i]-Sum[i]*Sum[i])/((numRow-1))))/Sum.at(i) << endl;
  }
  
  cfin.close(); 
  fout.close();
}
