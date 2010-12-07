#ifndef SIMPARAM_H
#define SIMPARAM_H

#include <fstream>
#include <vector>

class PARAMS
{
 public:
  int Dim_;           //lattice dimension (1,2,or3)
  vector <int> Size_; //vector of length Dim containing x,y,z dimensions
  string lattice_type;//Options are: square, kondo, J1J2, ...
  long SEED_;         //the random seed
  
  class RATIO{
  public:
    int Ncorners_;
    vector <int> Corners_;
    RATIO(); //constructor
  };//class RATIO

  PARAMS();//constructor

};//class PARAMS

PARAMS::RATIO::RATIO(){
  //specify the denominator region for the ratio
  //using the number of corners
  //NOTE: should the dimension of the region be specified too??
  ifstream rfin;
  rfin.open("ratio.txt");

  rfin >> Ncorners_;
  Corners_.resize(Ncorners_,-99);
  for(int i=0; i<Ncorners_; i++){rfin>>Corners_[i];}

  rfin.close();

}//RATIO constructor

PARAMS::PARAMS(){
  //reads in paramters from a file
  ifstream pfin;
  pfin.open("param.txt");

  pfin >> Dim_; //number of dimensions

  Size_.resize(Dim_,-99); //read in x,y,z dimensions
  for(int i=0; i<Dim_; i++){pfin>>Size_[i];}

  pfin >> lattice_type;
  pfin >> SEED_;

  
  pfin.close();

}//PARAMS constructor

#endif
