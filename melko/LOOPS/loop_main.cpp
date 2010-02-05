//Jan 18, 2010 --- starting loop code

#include "header.h"
#include "loop_header.h"

int main(){

  //  read in parameters: system dimensions, number of bond operators, 
  //   filenames, iterations per loop, number of loops, a random seed

  int dim1=6, dim2=1, its_per_loop=1, loops=1;
  double bops_per_site=0.5;
  bool OBC=1;
  long long ranseed=43289841;
  string filenames; //need to figure out what should be saved in files

  //read in params from file*******************

  int total_bops = dim1*dim2*bops_per_site;

  cout << dim1 << " x " << dim2 << " system,  N = " 
       << dim1*dim2 << " sites \n" << bops_per_site << " bops/site,  " 
       << total_bops << " bops total" << endl;
  cout << "------------------------------------------------ \n"; 

  if(dim1==2|dim2==2){cout<<"warning! nnbonds get screwed up for a x 2 \n";}
  LOOPS system (dim1, dim2, total_bops, OBC, ranseed);
 
  // create initial VB config and initial spin config
  system.nnbondlist();
  system.Nnnbondlist();
  system.generate_ops();   // generate 2m bond operators
  system.create_Vlinks();  // build vertical LL from init VBs and operators
  system.create__Hlinks(); // build horizontal linked list from operators
  system.make_flip_loops();// generate loops and flip w/ prob 0.5
  system.change_operators();// Change the diagonal operators


 

  
  

  //7  relabel operators and vertex types

  //8  measure things

  //9  change all diag bondops (go to 4)

  return 0;
}
