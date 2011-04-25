#ifndef BASIS_H
#define BASIS_H

#include "head_proj.h"
#include "simparam.h"


class Basis//: public PARAMS
{

	public:
		int LinX;
		int numSpin;
		int numVB;
		int numLattB;

		vector<int>  VBasis;   //VB basis

        Basis(const PARAMS &);
        Basis(const Basis &);  //copy constructor

        //The non-winding number fluctuating update
        void TwoBondUpdate(MTRand &, const PARAMS &, const vector<int> &);

		//The winding# fluctuating loop update
        void LoopUpdate(MTRand& , const PARAMS & , const vector<int> & );

        void print(); //print

        int TopoX(); //measures the X-topological sector of the VB wavefunction
        int TopoY(); //measures the Y-topological sector of the VB wavefunction

		//Basis operator=(const Basis & );
		int operator|(const Basis & ); //returns number of loops in overlap

		void filewrite(const int & num);
		void fileread(const int & num);

};

Basis::Basis(const PARAMS &p){//Square lattice constructor

	LinX = p.nX_;
	numLattB = p.numLattB;
	numSpin = p.numSpin;
	numVB = p.numVB;

	int a, b;
	//int x,y;
    index2 temp;
    for (int i=0; i<numSpin; i+=2){  //winding # (0.0)
        a = i;
        b = i+1;
        if ((i+1)%LinX == 0)
            b = a-LinX;
        VBasis.push_back(b);  //0 connected to 1
        VBasis.push_back(a);  //1 connected to 0
    }

	//winding # (1,0)
	//for (int i=0; i<LinX; i+=2){
    //    a = i;
    //    b = i-1;
	//	if (i == 0)
	//		b = LinX-1;
	//	VBasis.at(a) = b;
	//	VBasis.at(b) = a;
	//}

};

Basis::Basis(const Basis & B){//Copy constructor

	LinX = B.LinX;
	numLattB = B.numLattB;
	numSpin = B.numSpin;
	numVB = B.numVB;

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
        cout<<i<<"->"<<VBasis[i]<<endl;

    //cout<<"Is neighbor \n";
    //for(int i=0; i<numSpin; i++){
    //    for(int j=0; j<numSpin; j++)
    //        cout<<i<<" "<<j<<" "<<is_neighbor(i,j)<<endl;
    //}

};//print


//This measures the X-topolgical sector
int Basis::TopoX(){

    int topo=0;
    for (int i=0; i<LinX; i+=2){
        if (VBasis.at(i) == i+LinX)       //sublattice A->B
            topo += 1;
         if (VBasis.at(i+1) == i+1+LinX)  //sublattice B->A
            topo -= 1;
    }
    return topo;
}

//This measures the Y-topolgical sector
int Basis::TopoY(){

    int j;
    int topo=0;
    for (int i=0; i<LinX; i+=2){
        j = i*LinX;
        if (VBasis.at(j) == j+1)       //sublattice A->B
            topo += 1;
        if (VBasis.at(j+LinX) == j+1+LinX)  //sublattice B->A
            topo -= 1;
    }
    return topo;
}

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


//This function performs the two-bond update
//
//   D - C
//   |   |
//   A - B
//
void Basis::TwoBondUpdate(MTRand& ran, const PARAMS & p, const vector<int> & SS){

    int plaq;

    plaq = ran.randInt(p.Pst.size() - 1); //random spin state 0 or 1
    //cout<<"Plaq: "<<plaq<<endl;

    if (VBasis.at(p.Pst.at(plaq).A) == p.Pst.at(plaq).B && //bond connects A-B
            VBasis.at(p.Pst.at(plaq).C) == p.Pst.at(plaq).D){   //bond connects C-D

        //check to make sure spins are compatible on new bonds
        if ( (SS.at(p.Pst.at(plaq).A) != SS.at(p.Pst.at(plaq).D)) &&
                (SS.at(p.Pst.at(plaq).B) != SS.at(p.Pst.at(plaq).C)) ) { 

            VBasis.at(p.Pst.at(plaq).A) = p.Pst.at(plaq).D;
            VBasis.at(p.Pst.at(plaq).D) = p.Pst.at(plaq).A;
            VBasis.at(p.Pst.at(plaq).B) = p.Pst.at(plaq).C;
            VBasis.at(p.Pst.at(plaq).C) = p.Pst.at(plaq).B;
        }

    }
    else if (VBasis.at(p.Pst.at(plaq).A) == p.Pst.at(plaq).D && //bond connects A-D
            VBasis.at(p.Pst.at(plaq).C) == p.Pst.at(plaq).B){   //bond connects C-B

        //check to make sure spins are compatible on new bonds
        if ( (SS.at(p.Pst.at(plaq).A) != SS.at(p.Pst.at(plaq).B)) &&
                (SS.at(p.Pst.at(plaq).D) != SS.at(p.Pst.at(plaq).C)) ) { 

            VBasis.at(p.Pst.at(plaq).A) = p.Pst.at(plaq).B;
            VBasis.at(p.Pst.at(plaq).B) = p.Pst.at(plaq).A;
            VBasis.at(p.Pst.at(plaq).D) = p.Pst.at(plaq).C;
            VBasis.at(p.Pst.at(plaq).C) = p.Pst.at(plaq).D;
        }

    }

       //cout<<p.Pst.at(plaq).A<<" "<<p.Pst.at(plaq).B
        //    <<" "<<p.Pst.at(plaq).C<<" "<<p.Pst.at(plaq).D<<endl;

}
 
 
//This function performs the Loop Update
void Basis::LoopUpdate(MTRand& ran, const PARAMS & p, const vector<int> & SS){

	int origSite;
    origSite = ran.randInt(numSpin - 1); //random site to start

	cout<<"old : "<<origSite<<endl;
	int link = VBasis[origSite];
	cout<<"link: "<<link<<" "<<VBasis[link]<<endl;

    int oldsite;
    //int trial = ran.randInt(3); //choose one of 4 other nearest neighbors
	int nextSpin[3]; //pack an array with the 3 

    oldsite = origSite;
    int index=0;
	//pack the three new possible dimer positions onto an array
	if (p.Neighbor[link].A!= oldsite){ 
		nextSpin[index] = p.Neighbor[link].A; index++;}
	if (p.Neighbor[link].B!= oldsite){ 
		nextSpin[index] = p.Neighbor[link].B; index++;}
	if (p.Neighbor[link].C!= oldsite){ 
		nextSpin[index] = p.Neighbor[link].C; index++;}
	if (p.Neighbor[link].D!= oldsite){ 
		nextSpin[index] = p.Neighbor[link].D; index++;}

    int n0, n1, n2;
	n0 = nextSpin[0]; n1 = nextSpin[1]; n2 = nextSpin[2];

	for (int i=0; i<3; i++)
		cout<<nextSpin[i]<<" "; 
	cout<<endl;

	int temp = ran.randInt(5);
	//randomly reorder the array
	if (temp == 0) {nextSpin[0] = n0; nextSpin[1] = n1; nextSpin[2] = n2;}
	else if (temp == 1) {nextSpin[0] = n0; nextSpin[1] = n2; nextSpin[2] = n1;} 
	else if (temp == 2) {nextSpin[0] = n1; nextSpin[1] = n0; nextSpin[2] = n2;} 
	else if (temp == 3) {nextSpin[0] = n1; nextSpin[1] = n2; nextSpin[2] = n0;} 
	else if (temp == 4) {nextSpin[0] = n2; nextSpin[1] = n0; nextSpin[2] = n1;} 
	else if (temp == 5) {nextSpin[0] = n2; nextSpin[1] = n1; nextSpin[2] = n0;} 
	else cout<<"Reorder error \n";

	for (int i=0; i<3; i++)
		cout<<nextSpin[i]<<" "; 
	cout<<endl;

	int linkSpin = SS.at(link);
	cout<<"linkSpin: "<<linkSpin<<endl;

    int head;
	if (SS[nextSpin[0]] != linkSpin) head = nextSpin[0];
	else if (SS[nextSpin[1]] != linkSpin) head = nextSpin[1];
	else if (SS[nextSpin[2]] != linkSpin) head = nextSpin[2];
	else head = oldsite;

    cout<<head<<endl;



}//LoopUpdate
  

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

