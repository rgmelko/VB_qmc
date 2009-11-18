// Nov 4, 2009 // New double-projector vb qmc program

#include"header.h"
#include"header_dubproj.h"

int main()
{
  cout.precision(10);

  /*** READ IN PARAMETERS **********************************************/
  string enerfilename, bondlengthfilename, corrfilename;
  string bondopfile1, bondopfile2;
  int xsites, ysites, zsites, corrsite, change_number;
  int iterations_per_loop, loops;
  double bondops_per_site;
  long long ranseed;
    
  ifstream fin("param.txt"); // read in paramaters from file
  fin >> enerfilename >> bondopfile1 >> bondopfile2 >> bondlengthfilename
      >> corrfilename >> xsites >> ysites >> zsites >> corrsite
      >> bondops_per_site >> change_number
      >> iterations_per_loop >> loops >> ranseed;
  
  fin.close();

  int bondops = int(floor(xsites*ysites*zsites*bondops_per_site+.1));
  
  cout << endl;
  cout << xsites << " x " << ysites << " x " << zsites << " system" << endl;
  cout << "r   = " << change_number << "  bondops changed per step \n";
  cout << "n/N = " << bondops_per_site << "  bond operators / site" << endl;
  cout << "its = " << iterations_per_loop << "  iterations/loop for " 
       << loops << " loops" <<  endl << endl;

  /*** CREATE LATTICE ***(object that contains all the information)*****/

  LATTICE system (xsites, ysites, zsites, corrsite, bondops, change_number, 
		  iterations_per_loop, bondopfile1, bondopfile2, ranseed);

  /*** CREATE LIST OF NEAREST NEIGHBOUR BONDS **************************/

  system.nnbondlist();

  for(int SUPERLOOP=0; SUPERLOOP<loops; SUPERLOOP++)
    {
      ifstream fin3(bondopfile1.c_str());
      if (fin3.fail() ) { system.first_step(); }
      else{ system.read_bonds(); }
    
      for(int i=0; i<iterations_per_loop; i++)
	{
	  system.change_ops(i);
	  system.apply_ops(i);
	  system.decide(i);
	  system.measure_energy();
	  system.reinitialize(i);
	}

      system.calculate_stuff();

      cout << "C(L/2,L/2) = " << system.corrfinal << endl;
      cout << "energy = " << system.energy << endl;
      cout << "energy/site = " << system.energy/(xsites*ysites*zsites)<< endl;
      cout << system.accept/(iterations_per_loop*1.0)*100 
	   << "% accepted" << endl << endl;
    
 
      system.print_quantities(enerfilename,system.energy);
      system.print_quantities(corrfilename,system.corrfinal);
      system.print_bondops(bondopfile1,system.bondops1,system.offdiag1);
      system.print_bondops(bondopfile2,system.bondops2,system.offdiag2);

      system.super_initialize();
    }

  return 0;
}
