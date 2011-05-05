#ifndef RENYI_H
#define RENYI_H

#include "head_proj.h"
#include "basis.h"

class Renyi
{
    private: 
	  vector<double> entropy;
	  vector<double> TOTAL_H2;
      vector<vector<int> > inAreg; //definition of the A regions to measure
	  int nSwap;
	  int Nloop_den;   //number of loops in the denominator

    public: 
	  Renyi(const int &);  //direct <swap>
	  Renyi(const int &, const int &); //ratio 
      void zero();
      void measure_H2(const Basis &, const Basis &);
	  void measure_ratio(const Basis & , const Basis & );
      int calc_SWAP_2D(const Basis &, const Basis &, const int &);
      void record();
      void output(PARAMS &);
	  void printA(const int &);


};//Renyi

Renyi::Renyi(const int & nSpin){

    vector<int> Atemp;  //vector to be pushed back
    Atemp.assign(nSpin,0);

	ifstream fin;
	fin.open("regionA.dat");

	if (fin.fail() ) { //check for errors
		cout<<"Could not open a regionA.dat file"<<endl;
	}

    fin>>nSwap;
    if (nSwap < 1) cout<<"regionA.dat error 1 \n";

    int temp;

    for (int j=0; j<nSwap; j++){
        for (int i=0; i<nSpin/2; i++){
            fin>>temp;
            if (temp != 0 && temp != 1)  cout<<"regionA.dat error 2 \n";
            Atemp.at(i) = temp; //base layer
            Atemp.at(i+nSpin/2) = temp; //ancillary layer
        }
        fin>>temp;
        if (temp != -99) cout<<"regionA.dat error 3 \n";
        inAreg.push_back(Atemp); //the 1 spin region
    }//j

    fin.close();

	entropy.assign(nSwap,0);  //resize and initialize entropy
	TOTAL_H2.assign(nSwap,0);  //resize and initialize entropy total

	//--if there is a regionX.dat file, push this back as the last
	//--vector element in inAreg<>
	fin.open("regionX.dat");
	fin>>temp;
	if ( temp==0 )  //check for errors
		cout<<"not using Ratio! : renyi"<<endl;
	else{
		for (int i=0; i<nSpin/2; i++){
			fin>>temp;
			if (temp != 0 && temp != 1)  cout<<"regionX.dat error 3  \n";
			Atemp.at(i) = temp; //base layer
			Atemp.at(i+nSpin/2) = temp; //ancillary layer
		}
		fin>>temp;
		if (temp != -99) cout<<"regionX.dat error 4  BASIS\n";
		inAreg.push_back(Atemp); //the 1 spin region
	}

	fin.close();

};//constructor


//Renyi::Renyi(const int & Lsize, const int & x_num){
//
//	nSwap = Lsize-x_num; //exclude some number of points for the overlap
//
//	entropy.assign(nSwap-1,0);  //resize and initialize entropy
//	TOTAL_H2.assign(nSwap-1,0);  //resize and initialize entropy total
//
//};//constructor

void Renyi::zero(){

	TOTAL_H2.assign(nSwap,0);

}//zero


void Renyi::measure_H2(const Basis & A, const Basis & B){

	Basis Vl(A);   //copy constructors
	Basis Vr(B);
    Nloop_den = Vl|Vr ; 

    int Nloop_num;
	for (int i=0; i<entropy.size(); i++){
		Nloop_num = calc_SWAP_2D(A,B,i);
		entropy.at(i) = 1.0*pow(2,Nloop_num - Nloop_den);
	}

}//measure_H2


void Renyi::measure_ratio(const Basis & A, const Basis & B){

	Basis Vl(A);   //copy constructors
	Basis Vr(B);
	Nloop_den = calc_SWAP_2D(Vl,Vr,nSwap);  //the last element in inAreg

    int Nloop_num;
	for (int i=0; i<entropy.size(); i++){
        Nloop_num = calc_SWAP_2D(A,B,i);
		entropy.at(i) = 1.0*pow(2,Nloop_num - Nloop_den);
	}

}//measure_ratio


int Renyi::calc_SWAP_2D(const Basis & A, const Basis & B, const int & X){

	//for (int i=0; i < inAreg.size(); i++)
	//	cout<<i<<" "<<inAreg.at(i)<<endl;

	int Nloop_num; //number of loops in numerator and denominator
	Basis Vl(A);   //copy constructors
	Basis Vr(B);
    //Nloop_den = Vl|Vr ; 

	int a,b, bond1, bond2;
	int old1, old2, old3, old4;

    int sA; //spin in region "A"
    for (int i=0; i<A.numSpin/2; i++){ 

        sA = i;
        if (inAreg[X].at(sA) == 1) {//if the base layer spin is in region A

            a = sA;
            b = sA + A.numSpin/2; //b in the other layer

            bond1 = Vr.VBasis[a]; 
            bond2 = Vr.VBasis[b];

            if (inAreg[X].at(bond2) == 1 )
                Vr.VBasis[a] = bond2 - A.numSpin/2;
            else{
                Vr.VBasis[a] = bond2;
                Vr.VBasis[bond2] = a;
            }

            if (inAreg[X].at(bond1) == 1 )
                Vr.VBasis[b] = bond1 + A.numSpin/2;
            else{
                Vr.VBasis[b] = bond1;
                Vr.VBasis[bond1] = b;
            }

        }//if in A
    }//i

    Nloop_num = Vl|Vr; //new overlap

    return Nloop_num;

}//calc_SWAP_2D


void Renyi::record(){

    for (int i=0; i<TOTAL_H2.size(); i++)
	   TOTAL_H2.at(i) += entropy.at(i);

}//record


void Renyi::output(PARAMS & p){

	ofstream cfout;
    cfout<<setprecision(8);
	cfout.open("00.renyi",ios::app);

    for (int i=0; i<TOTAL_H2.size(); i++)
	  //cfout<<-log(TOTAL_H2.at(i)/(1.0*p.MCS_))<<" "; //careful with averaging logs
	  cfout<<TOTAL_H2.at(i)/(1.0*p.MCS_)<<" "; //careful with averaging logs
	cfout<<endl;

	cfout.close();

}//output


void Renyi::printA(const int & nX){

	for (int j=0; j<inAreg.size(); j++){
		for (int i=0; i<2*nX*nX; i++){
			cout<<inAreg[j].at(i)<<" ";
			if ((i+1)%nX == 0) cout<<endl;
		}
		cout<<endl;
	}//j

}//print


#endif
