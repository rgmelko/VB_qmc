#ifndef BASIS_H
#define BASIS_H

#include "head_proj.h"
#include "simparam.h"
#include "projector.h"


class Basis: public PARAMS
{

    private: 

        //vector<index2>  VBasis;   //VB basis
        vector<int>  VBasis;   //VB basis

    public:

        double Weight;

        Basis();
        void print();
        void Propogate(const Projector&, Basis&);

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
        //x = a%nX_; y = a/nX_;
        //if ((x+y)%2 == 0) temp.set(a,b);
        //else temp.set(b,a);
        //VBasis.push_back(temp);
        VBasis.push_back(b);  //0 connected to 1
        VBasis.push_back(a);  //1 connected to 0
    }
};

void Basis::print(){

    cout<<"VB basis: \n";
    for (int i=0;  i<VBasis.size(); i++)
        //VBasis[i].print();
        cout<<VBasis[i]<<endl;


};//print

void Basis::Propogate(const Projector& P, Basis& beta){

    for (int i=0; i<(*this).VBasis.size(); i++)
        beta.VBasis.at(i) = (*this).VBasis.at(i);

    Weight = 1;
    
    int a,b;
    int bond1, bond2;
    for (int i=0; i<P.list_size; i++){
        a=Bst.at(P.O_list[i]).A;
        b=Bst.at(P.O_list[i]).B;

        bond1 = beta.VBasis[a]; 
        bond2 = beta.VBasis[b];

        if (a == bond2) {
            if (b != bond1) cout<<"VB connection error \n";
        }
        else{
            beta.VBasis[a] = b;
            beta.VBasis[b] = a;
            beta.VBasis[bond1] = bond2;
            beta.VBasis[bond2] = bond1;
            Weight *= 0.5;
        }

    }//i

}//Propogate


#endif 

