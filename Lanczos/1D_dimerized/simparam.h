#ifndef PARAMSIM_H
#define PARAMSIM_H


//Class to read in the simulation parameters from a file
// Adapted to the JQ ED code

class PARAMS
{
  public:
    double JJ_; //the heisenberg exchange
    double J2_; //the dimerized exchange
    int Sz_; // z-component of total spin
    int L_;  //1D system size
    // FOR LANCZOS
    int Neigen_;    //: # of eigenvalues to converge
    int valvec_; //  1 for -values only, 2 for vals AND vectors
    // FULL_DIAG?

    PARAMS(){
      //initializes commonly used parameters from a file
      ifstream pfin;
      pfin.open("param.dat");
    
      pfin >> L_;
      pfin >> JJ_;
      pfin >> J2_;
      pfin >> Sz_;
      pfin >> Neigen_;
      pfin >> valvec_;
    
      pfin.close();
    }//constructor

}; //PARAMS

#endif
