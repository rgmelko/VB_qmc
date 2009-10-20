#ifndef BASIS_H
#define BASIS_H

#include "head_proj.h"
#include "simparam.h"
#include "projector.h"


class Basis: public PARAMS
{

    private: 

        iMatrix is_neighbor; //alternate lookup array

    public:
        vector<int>  VBasis;   //VB basis

        long int Weight;
        double Energy;

        Basis();
        void print();
        void Propogate(const Projector&, Basis&);
        double Calc_Energy();

        int operator|(const Basis & ); //returns number of loops in overlap



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

    //initialize neighbor list (for energy calc)
    is_neighbor.resize(numSpin,numSpin);
    for (int i=0; i<numSpin; i++)
        for (int j=0; j<numSpin; j++)
            is_neighbor(i,j) = 0;

    for (int i=0; i<Bst.size(); i++){
        is_neighbor(Bst[i].A,Bst[i].B) = 1;
        is_neighbor(Bst[i].B,Bst[i].A) = 1;
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

double Basis::Calc_Energy(){

    double n_ij = 0;
    int a,b;
    for (int i=0; i<VBasis.size(); i++){
        a = i;
        b = VBasis.at(i);

        if (is_neighbor(a,b) == 1)
           n_ij += 1; 
    }//i

    n_ij *= 0.5;   //think I've double counted
    n_ij += 1; //add the +1

    return n_ij;

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
    


#endif 

