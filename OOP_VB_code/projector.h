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

        Projector();           //constructor 1
        Projector(const int);  //constructor 2
        void print();
        void Sample_Ops(); //swap a number of operators

        Projector operator=(const Projector & Pj);


};

Projector::Projector(){ //overloaded constructor 1

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


void Projector::Sample_Ops(){

   int element, oper;
   for (int i=0; i<sample_; i++){
       element = ran.randInt(O_list.size()-1);
       oper = ran.randInt(Bst.size() - 1); //random bond operator
       O_list.at(element) = oper;
   }


}//Sample_Ops

Projector Projector::operator=(const Projector & Pj){

    for (int i=0; i<list_size; i++)
        O_list.at(i) = Pj.O_list.at(i);

    return *this;


}

#endif
