// Nov 4, 2009 // New double-projector vb qmc program
// Nov 20,2009 // Trying to adapt this prog to measure renyi ent
// Apr 1, 2010 // adapting to use ratio trick

#include"header.h"
#include"ratio_header.h"

int main()
{
  cout.precision(10);

  /*** READ IN PARAMETERS **********************************************/
  int its_per_measurement = 100;
  string enerfilename, entrfilename;
  string bondopfile1, bondopfile2;
  int xsites, ysites, zsites, change_number;
  int iterations_per_loop, loops;
  double bondops_per_site;
  long long ranseed;
    
  ifstream fin("param.txt"); // read in paramaters from file
  fin >> enerfilename >> entrfilename >> bondopfile1 >> bondopfile2 
      >> xsites >> ysites >> zsites
      >> bondops_per_site >> change_number
      >> iterations_per_loop >> loops >> ranseed;
  
  fin.close();
 
  int bondops = int(floor(xsites*ysites*zsites*bondops_per_site+.1));
 
  /*** OUTPUT PARAMETERS *********************************************/
  cout << endl;
  cout << xsites << " x " << ysites << " x " << zsites << " system" << endl;
  cout << "r   = " << change_number << "  bondops changed per step \n";
  cout << "n/N = " << bondops_per_site << "  bond operators / site" << endl;
  cout << "its = " << iterations_per_loop << "  iterations/loop for " 
       << loops << " loops" <<  endl << endl;

  /*** CREATE LATTICE ***(object that contains all the information)*****/

  LATTICE system (xsites, ysites, zsites, bondops, change_number, 
		  iterations_per_loop, bondopfile1, bondopfile2, 
		  its_per_measurement, ranseed);

  /*** CREATE LIST OF NEAREST NEIGHBOUR BONDS **************************/

  system.nnbondlist();  //makes list of nn bonds

  for(int SUPERLOOP=0; SUPERLOOP<loops; SUPERLOOP++)
    {
      ifstream fin3(bondopfile1.c_str()); //If bondopfile1 doesn't exist yet
      if (fin3.fail() ) { system.first_step(); } //run the "first step"

      else{ system.read_bonds(); } //Otherwise, read in the last configuration
    
      for(int i=0; i<iterations_per_loop/its_per_measurement; i++)
	{
	  for(int j=0; j<its_per_measurement-1; j++){
	    system.change_ops(j);    //Change a few operators 
	    system.apply_ops(j);     //Apply them to the trial state
	    system.decide(j);        //Decide whether to keep the changes
	    system.reinitialize(j);  //Reinitialize some things
	  }
	  system.change_ops(its_per_measurement-1);//Change a few operators 
	  system.apply_ops(its_per_measurement-1);//Apply to the trial state
	  system.decide(its_per_measurement-1);//Decide if we keep the changes
	  system.measure_energy();//Measure the energy
	  system.measure_swap();//Measure the swap 
	  system.reinitialize(its_per_measurement-1);//Reinitialize some things
	}

   

      system.calculate_stuff();  // After "loops" iterations, calculate stuff

      cout << "energy = " << system.energy*.5 << endl;
      cout << "energy/site = " << system.energy/system.number_of_sites<< endl;
      cout << system.accept/(iterations_per_loop*1.0)*100 
	   << "% accepted" << endl << endl;
    
 
      system.print_quantities(enerfilename,system.energy);
      system.print_entropies(entrfilename,system.entropy);
      system.print_bondops(bondopfile1,system.bondops1,system.offdiag1);
      system.print_bondops(bondopfile2,system.bondops2,system.offdiag2);

      system.super_initialize();
    }

  return 0;
}
