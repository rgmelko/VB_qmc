//hopefully i'll be able to measure the vbEE of a ladder with OBC and 
//an arbitrary number of legs

#include<iostream>
#include<fstream>
#include<math.h>
#include"ladder_header.h" //my ladder class
#include"mtrand.h" // random number generator

using namespace std;

int main()
{
  cout << endl;
  
  int legs, length; // system dimensions
  int y; // number of bond operators per site  
  int n; // total number of bond operators
  int r; // number of bond operators changed per MC step

  ifstream fin("param.txt"); // read in paramaters from file
  fin >> legs >> length >> y >> r;
  fin.close();
  
  n = legs*length*y;

  cout << legs << " legs, " << length << " sites long each" << endl;
  cout << "r = " << r << "     " << "y = " << y << "     "<< "n = " << n <<"\n\n";

  LADDER system (legs, length, n, r);

  system.nnbondlist();


  // PRINTS OUT THE POSSIBLE NN BONDS ------------------- 
    for(int i=0; i< system.nnbonds0.size(); i++)
      {
        cout << system.nnbonds0[i] << "," << system.nnbonds1[i] << endl;
      }

  // initial state stuff..

  //read in bond operators from file.. if file is empty generate n 
  // bond operators
  // also read in weights from the last step (number of non-diag ops
  // from the last step)
 
  // if program hasn't been run yet (# steps = 0) warm up
  // I possibly don't need to record numbers during the warm up

  // first MC step (if # steps was zero)
    // apply n bond operators
    // count nnbonds, entropy crossings, weight(# non-diag ops)
    // change r bond operators

  // non-first MC  steps (first step if # steps was not zero)
    // apply n bond operators
    // count nnbonds, entropy crossing, weight
    // use weights to get a probability of accepting changes
    // generate random number to decide whether to keep changes
    // if YES (keep new entropy, energy & add to totals, keep weight)
    // if NO (discard new measurements & weight, record #s from the previous step)
    // change r bond operators from the new(if YES)/old(if NO) list of bond ops

  //At the end of some number of steps:
    // Calculate energy & entropy & output them to data file
    // output: # steps, bond operators, weights

  

  return 0;
}
