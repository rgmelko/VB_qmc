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
  int its;
 
  ifstream fin("param.txt"); // read in paramaters from file
  fin >> legs >> length >> y >> r >> its;
  fin.close();
  
  n = legs*length*y;

  cout << legs << " x " << length << " system" << endl;
  cout << "r = " << r << "     " << "y = " << y << "     "<< "n = " << n;
  cout << "   " << its << " iterations" << endl;

  LADDER system (legs, length, n, r, its);

  system.nnbondlist();

  /* PRINTS OUT THE POSSIBLE NN BONDS ------------------- 
     for(int i=0; i< system.bonds.size(); i++)
     {
     cout << i << "," << system.init[i] << endl;
     }
  */
  system.first_step();
  
  for(int i=0; i<its; i++)
    {
      system.change_ops();
      system.apply_ops();
      system.decide();
      system.measure();
      system.reinitialize();
    }
  
  system.calculate_stuff();

  cout << endl << system.accept/its*100 << "% accepted" << endl;

  cout << "energy = " << system.energy << endl;
  cout << "for zone(2) entropy = " << system.entropies[1] << endl;

  // PRINTS THE SUPER AWESOME BOND CHECKING MATRIX ------------
  /*
    for(int i=0; i< system.nncheck.length(); i++)
    {
    for(int j=0; j< system.nncheck.width(); j++)
    {
    cout << system.nncheck(i,j) << "," ;
    }
    cout << endl;
    }
  */
  

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
  cout << endl;
  return 0;
}
