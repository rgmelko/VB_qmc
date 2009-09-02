//
// VB_oop.cpp: the main body of the program
// A OOP VB basis projector code 
// Roger Melko, Aug 29 2009, Santa Barbara
//
//
//
#include "head_proj.h"
#include "simparam.h"
#include "basis.h"
#include "projector.h"

int main(){

    MTRand met_rand(2343983); //random number for metropolis

    Projector P1; //initialize projector 1
    //P1.print();
    Projector P2; //initialize projector 2
    //P2.print();

    Basis alpha;
    //alpha.print();
  
    long int W_old, W_new; //old and new weights
    double DeltaW;
    P2 = P1;  //set projectors equal
    
    Basis beta;
    alpha.Propogate(P1,beta);
    W_old =  beta.Weight;

    //EQUIL
    for (int i=0; i<2000; i++){
        P2.Sample_Ops();
        alpha.Propogate(P2,beta);
        W_new =  beta.Weight;
        DeltaW = pow(2,W_old - W_new);

        if (DeltaW >= 1){//keep changes
            W_old = W_new;
            P1 = P2;
            //cout<<-beta.Energy<<endl;
        }
        else if (DeltaW > met_rand.rand()){
            W_old = W_new;
            P1 = P2;
            //cout<<-beta.Energy<<endl;
        }
        else P2 = P1;
    }//EQUIL

    //MCS
    int MCS = 250000;
    double E_new, E_old;
    double energy = 0;
    alpha.Propogate(P1,beta);
    W_old =  beta.Weight;
    E_old = beta.Energy;

    for (int i=0; i<MCS; i++){
        P2.Sample_Ops();
        alpha.Propogate(P2,beta);
        W_new =  beta.Weight;
        E_new = beta.Energy;
        //E_new = beta.Calc_Energy();

        if (W_new > W_old){//keep changes
            W_old = W_new;
            E_old = E_new;
            P1 = P2;
        }
        else if (W_new/W_old > met_rand.rand()){
            W_old = W_new;
            E_old = E_new;
            P1 = P2;
            //cout<<-beta.Energy<<endl;
        }
        else P2 = P1;
        energy -= E_old;
    }//MCS
    cout<<energy/MCS<<endl;


  return 0;
};