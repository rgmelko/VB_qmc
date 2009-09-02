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
    
    Basis beta;

    //MCS
    int MCS = 500000;
    double E1_new, E1_old;
    double E2_new, E2_old;
    double energy1, energy2;
    for (int EQMC = 0; EQMC <2; EQMC++) { //EQL and MCS run loop
        P2 = P1;  //set projectors equal
        energy1 = 0;
        energy2 = 0;
        alpha.Propogate(P1,beta);
        W_old =  beta.Weight;
        E1_old = -beta.Energy;
        E2_old = beta.Calc_Energy();

        for (int i=0; i<MCS; i++){

            P2.Sample_Ops();
            alpha.Propogate(P2,beta);
            W_new =  beta.Weight;
            DeltaW = pow(2,W_old - W_new);
            //cout<<W_old<<" "<<W_new<<" "<<DeltaW<<" "<<met_rand.rand()<<endl;
            E1_new = - beta.Energy;
            E2_new = beta.Calc_Energy();

            if (DeltaW >= 1){//keep changes
                W_old = W_new;
                E1_old = E1_new;
                E2_old = E2_new;
                P1 = P2;
            }
            else if (DeltaW > met_rand.rand()){
                W_old = W_new;
                E1_old = E1_new;
                E2_old = E2_new;
                P1 = P2;
                //cout<<-beta.Energy<<endl;
            }
            else P2 = P1;

            energy1 += E1_old;
            energy2 += E2_old;
        }//MCS

    }//EQMC
    energy1 /= MCS;
    energy2 /= MCS;
    energy2 += beta.numLattB/4.0;

    cout<<energy1<<endl;
    cout<<energy2<<endl;


  return 0;
};