#ifndef PROJETOR_H
#define PROJETOR_H

#include "simparam.h"
#include "head_proj.h"


class Projector: public PARAMS
{

    public:
        vector<int>  O_list;  //Operator list, the Bst number
        //vector<index2>  Basis;   //VB basis

        int list_size;

        Projector(MTRand& ran);           //constructor 1
        Projector(const int);  //constructor 2
        void print();
        void Sample_Ops(MTRand& ran); //swap a number of operators

        Projector operator=(const Projector & Pj);

        void filewrite(const int & num);
        void fileread(const int & num);


};

Projector::Projector(MTRand& ran){ //overloaded constructor 1

    list_size = NN_;
    //cout<<"inside "<<N<<" \n";

    //printBst();

    int temp;
    for (int i=0; i<list_size; i++){
        temp = ran.randInt(Bst.size() - 1); //random bond operator
        O_list.push_back(temp);
    }


}//constructor

//Projector::Projector(const int N){ //overloaded constructor 2
//
//
//}//constructor

void Projector::print(){

    cout<<"List size is: "<<list_size<<endl;
    //cout<<O_list.size()<<endl;
    for (int i=0; i<list_size; i++){
        cout<<O_list.at(i)<<" ";
        Bst[O_list.at(i)].print();
    }

}//print


void Projector::Sample_Ops(MTRand& ran){

   int element, oper, old;
   for (int i=0; i<sample_; i++){
       element = ran.randInt(O_list.size()-1);
       old = O_list.at(element);
       do {
          oper = ran.randInt(Bst.size() - 1); //random bond operator
       } while (oper == old);
       O_list.at(element) = oper;
   }


}//Sample_Ops

Projector Projector::operator=(const Projector & Pj){

    for (int i=0; i<list_size; i++)
        O_list.at(i) = Pj.O_list.at(i);

    return *this;

}

void Projector::filewrite(const int & num){

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
	fname[3] = 'c';
	fname[4] = 'o';
	fname[5] = 'n';
	fname[6] = 'f';
	fname[7] = '\0';


	ofstream cfout;
	cfout.open(fname);

    cfout<<list_size<<endl;
    for (int i=0; i<list_size; i++){
        cfout<<O_list.at(i)<<" ";
    }

	cfout.close();

}//filewrite


void Projector::fileread(const int & num){

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
	fname[3] = 'c';
	fname[4] = 'o';
	fname[5] = 'n';
	fname[6] = 'f';
	fname[7] = '\0';


	ifstream cfin;
	cfin.open(fname);

	O_list.clear();

    int temp;
    cfin>>list_size;
    for (int i=0; i<list_size; i++){
		cfin>>temp;
        O_list.push_back(temp);
    }
	if (O_list.size() !=  list_size) cout<<"READ IN ERROR \n";

	cfin.close();

}//fileread

#endif
