// measure the vbEE of a ladder with PBC and 
//for only 1D

#include<iomanip>
#include"header.h"
#include"ladderPBC_header.h" //my ladder class


using namespace std;

int main()
{
    string entfilename, enerfilename, bondopfile;
    int legs, length; // system dimensions
    int y; // number of bond operators per site  
    int n; // total number of bond operators
    int r; // number of bond operators changed per MC step
    int its;
    int loops;
    long long ranseed;

    ifstream fin("param_PBC.txt"); // read in paramaters from file
    fin >> enerfilename  >> entfilename >> bondopfile 
    	>> legs >> length >> y >> r >> its >> loops 
	>> ranseed;
    
    fin.close();

    n = legs*length*y;

    cout << legs << " x " << length << " system" << endl;
    cout << "r = " << r << "     " << "y = " << y << "     "<< "n = " << n;
    cout << "   " << its << " iterations" << endl;

    LADDER system (legs, length, n, r, its, bondopfile, ranseed);
    system.nnbondlist();

    /* PRINTS OUT THE POSSIBLE NN BONDS ------------------- 
       for(int i=0; i< system.bonds.size(); i++)
       {
       cout << i << "," << system.init[i] << endl;
       }
     */

    for(int SUPERLOOP=0; SUPERLOOP<loops; SUPERLOOP++)
    {
	ifstream fin3(bondopfile.c_str());
	if (fin3.fail() ) {system.first_step();}
	else{ system.read_bonds();}

	for(int i=0; i<its; i++)
	{
	    system.change_ops();
	    system.apply_ops();
	    system.decide();
	    system.measure();
	    system.reinitialize();
	}

	system.calculate_stuff();

	ofstream fout(enerfilename.c_str(),ios::app);
	fout << setprecision(10) << system.energy << endl;
	fout.close();

	ofstream fout2(entfilename.c_str(),ios::app);
	for(int i=0; i<system.entropies.size(); i++)
	{
	    fout2 << setw(12) << setprecision(10) << system.entropies[i] << " ";
	}
	fout2 << endl;
	fout2.close();

	ofstream fout3(bondopfile.c_str());
	for(int i=0; i<system.bondopsA.size(); i++)
	{
	    fout3 << system.bondopsA[i] << endl;
	}
	fout3 << system.offdiagA << endl;
	fout3.close();

	cout << endl << system.accept/its*100 << "% accepted" << endl;

	cout << "energy = " << setprecision(10) <<system.energy << endl;
	cout << "for zone(2) entropy = " << system.entropies[1] << endl;

	cout << endl;

	system.super_initialize();
    }

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

    return 0;
}
