#ifndef BASIS_H
#define BASIS_H

#include "head_proj.h"
#include "simparam.h"
#include "projector.h"


class Basis: public PARAMS
{

    public:
        vector<int>  VBasis;   //VB basis

        long int Weight;
        double Energy;

        Basis();
        Basis(const Basis &);  //copy constructor
        void print();
        void Propogate(const Projector&, Basis&);

		//Basis operator=(const Basis & );
		int operator|(const Basis & ); //returns number of loops in overlap

		void Basis::filewrite(const int & num);
		void Basis::fileread(const int & num);

};

Basis::Basis(){//Square lattice constructor

    int a, b;
    //int x,y;
    index2 temp;
    for (int i=0; i<numSpin; i+=2){
        a = i;
        b = i+1;
        if ((i+1)%nX_ == 0)
            b = a-nX_;
        VBasis.push_back(b);  //0 connected to 1
        VBasis.push_back(a);  //1 connected to 0
    }

};

Basis::Basis(const Basis & B){//Copy constructor

	int temp;
    for (int i=0; i<B.VBasis.size(); i++){
		temp = B.VBasis.at(i);
        (*this).VBasis.push_back(temp);  //0 connected to 1
    }

};

void Basis::print(){

    cout<<"VB basis: \n";
    for (int i=0;  i<VBasis.size(); i++)
        //VBasis[i].print();
        cout<<VBasis[i]<<endl;

    //cout<<"Is neighbor \n";
    //for(int i=0; i<numSpin; i++){
    //    for(int j=0; j<numSpin; j++)
    //        cout<<i<<" "<<j<<" "<<is_neighbor(i,j)<<endl;
    //}

};//print

void Basis::Propogate(const Projector& P, Basis& beta){

    for (int i=0; i<(*this).VBasis.size(); i++)
        beta.VBasis.at(i) = (*this).VBasis.at(i);

    beta.Weight = 0;
    beta.Energy = 0;

    int a,b;
    int bond1, bond2;
    for (int i=0; i<P.list_size; i++){
        a=Bst.at(P.O_list[i]).A;
        b=Bst.at(P.O_list[i]).B;

        bond1 = beta.VBasis[a]; 
        bond2 = beta.VBasis[b];

        //diagonal operation
        if (a == bond2) {
            if (b != bond1) cout<<"VB connection error \n";
            beta.Energy += 1;
        }
        //off diagonal operation
        else{
            beta.VBasis[a] = b;
            beta.VBasis[b] = a;
            beta.VBasis[bond1] = bond2;
            beta.VBasis[bond2] = bond1;
            beta.Weight += 1;
            beta.Energy += 0.5;
        }

    }//i
    beta.Energy /= 1.0*P.list_size;

}//Propogate


int Basis::operator|(const Basis & B){

    vector<int> is_in_loop;  //records whether a spin is counted in a loop 
    is_in_loop.assign(B.VBasis.size(),0);

    int next;
    int Nloop = 0;

    for (int i=0; i<B.VBasis.size(); i++){

        if (is_in_loop.at(i) == 0){
            is_in_loop.at(i) = 1;
            next = (*this).VBasis.at(i); //V_A basis
            while (is_in_loop.at(next) == 0){

                if  (is_in_loop.at(next) == 1) cout<<"loop error 1 \n";
                else is_in_loop.at(next) = 1;

                next = B.VBasis.at(next);      //V_B basis
                if  (is_in_loop.at(next) != 1) is_in_loop.at(next) = 1; 
                else break;

                next = (*this).VBasis.at(next); //V_A basis 
            }//while

            Nloop ++;

        }//if
    }//i

    return Nloop;

} //operator |
    

void Basis::filewrite(const int & num){

	char fname[8];

	if (num == 0) fname[1] = '0';
	else if (num%9 == 0) fname[1] = '9';
	else if (num%8 == 0) fname[1] = '8';
	else if (num%7 == 0) fname[1] = '7';
	else if (num%6 == 0) fname[1] = '6';
	else if (num%5 == 0) fname[1] = '5';
	else if (num%4 == 0) fname[1] = '4';
	else if (num%3 == 0) fname[1] = '3';
	else if (num%2 == 0) fname[1] = '2';
	else if (num%1 == 0) fname[1] = '1';

	fname[0] = '0';
	fname[2] = '.';
	fname[3] = 'b';
	fname[4] = 'a';
	fname[5] = 's';
	fname[6] = 'e';
	fname[7] = '\0';


	ofstream cfout;
	cfout.open(fname);

    for (int i=0;  i<VBasis.size(); i++)
        cfout<<VBasis[i]<<endl;

	cfout.close();

}//filewrite

void Basis::fileread(const int & num){

	char fname[8];

	if (num == 0) fname[1] = '0';
	else if (num%9 == 0) fname[1] = '9';
	else if (num%8 == 0) fname[1] = '8';
	else if (num%7 == 0) fname[1] = '7';
	else if (num%6 == 0) fname[1] = '6';
	else if (num%5 == 0) fname[1] = '5';
	else if (num%4 == 0) fname[1] = '4';
	else if (num%3 == 0) fname[1] = '3';
	else if (num%2 == 0) fname[1] = '2';
	else if (num%1 == 0) fname[1] = '1';

	fname[0] = '0';
	fname[2] = '.';
	fname[3] = 'b';
	fname[4] = 'a';
	fname[5] = 's';
	fname[6] = 'e';
	fname[7] = '\0';


	ifstream cfin;
	cfin.open(fname);

	int temp=0;
    for (int i=0;  i<VBasis.size(); i++){
        cfin>>VBasis[i];
		temp++;
	}
	if (temp != VBasis.size()) cout<<"Basis Read ERROR \n";

	cfin.close();

}//fileread

#endif 

