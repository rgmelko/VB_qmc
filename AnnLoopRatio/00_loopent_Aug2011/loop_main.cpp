//Jan 18, 2010 --- starting loop code

#include "header.h"
#include "loop_header.h"

int main(){

  cout.precision(10);
  //  read in parameters: system dimensions, number of bond operators, 
  //   filenames, iterations per loop, number of loops, a random seed

  int dim1, dim2;
  long long its_per_loop=10000, loops=100;
  long long initialization;
  double bops_per_site=10;
  bool OBC=0;
  long long ranseed=43289841;
  string enerfilename, entrofilename, bondopfilename;

  ifstream fin("param.in");
  fin >> enerfilename >> entrofilename >> bondopfilename
      >> dim1 >> dim2
      >> its_per_loop >> loops
      >> bops_per_site >> OBC
      >> ranseed
      >> initialization;
  fin.close();
  
  int total_bops = dim1*dim2*bops_per_site;

  cout << dim1 << " x " << dim2 << " system,  N = " 
       << dim1*dim2 << " sites \n" << bops_per_site << " bops/site,  " 
       << total_bops << " bops total" << endl;
  cout << "------------------------------------------------ \n"; 

  if(dim1==2|dim2==2){cout<<"warning! nnbonds get screwed up for a x 2 \n";}
  LOOPS system (dim1, dim2, total_bops, OBC, its_per_loop, ranseed, 
		bondopfilename);
 
  // create initial VB config and initial spin config
  system.nnbondlist();
  system.Nnnbondlist();
  system.read_bops(); //checks if file has bops, otherwise generates new ones

  for(int jj=0; jj<initialization; jj++){
    system.create_Vlinks();
    system.create__Hlinks();
    system.make_flip_loops();
    system.change__operators();
  }
  for(int kk=0; kk<loops; kk++){
    for(int jk=0; jk<its_per_loop; jk++){
      system.create_Vlinks();    //build vertical LL from init VBs and operators
      system.create__Hlinks();   //build horizontal linked list from operators
      system.make_flip_loops();  //generate loops and flip w/ prob 0.5
      system.take_measurement();
      system.swaperator();
      system.change__operators(); //Change the diagonal operators
    }

    system.calculate_stuff();

    ofstream energy_out(enerfilename.c_str(),ios::app);
    ofstream entrpy_out(entrofilename.c_str(),ios::app);
    energy_out.precision(10);
    entrpy_out.precision(10);

    cout << left << setw(12) << system.energy << "    ";
    energy_out << system.energy << endl;

    cout << system.entropy_final[dim1-2] << endl;
    energy_out.close();
    
    
    for(int i=0; i<system.entropy_final.size(); i++){
      entrpy_out << setw(18) << system.entropy_final[i];
    }
    entrpy_out << endl;
    
    entrpy_out.close();
    system.print_bops();
  }

  return 0;
}
