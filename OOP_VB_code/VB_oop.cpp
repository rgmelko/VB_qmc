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

    PARAMS param;
    MTRand mrand(param.SEED_); //random number for metropolis

    Projector P1(mrand); //initialize projector 1
    //P1.print();
    Projector P2(mrand); //initialize projector 2
    //P2.print();

    Basis alpha;
    //alpha.print();
  
    long int W_old, W_new; //old and new weights
    double DeltaW;
    
    Basis beta;
    //beta.print();
    //return 0;

    //MCS
    int MCS = 100000;
    double E1_new, E1_old;
    double E2_new, E2_old;
    double energy1, energy2;
    double bondtot_new, bondtot_old;
    for (int EQMC = 0; EQMC <2; EQMC++) { //EQL and MCS run loop
        MCS *= 2;
        P2 = P1;  //set projectors equal
        energy1 = 0;
        energy2 = 0;
		alpha.Propogate(P1,beta);
		W_old =  beta.Weight;

		E1_new = beta.Calc_Energy()-1;
		E2_old = beta.Calc_Energy();

        for (int i=0; i<MCS; i++){

            P2.Sample_Ops(mrand);
            alpha.Propogate(P2,beta);
            W_new =  beta.Weight;
            DeltaW = pow(2,W_old - W_new);
            //cout<<W_old<<" "<<W_new<<" "<<DeltaW<<" "<<met_rand.rand()<<endl;
            E1_new = beta.Calc_Energy()-1;
            E2_new = beta.Calc_Energy();

            if (DeltaW > mrand.rand()){
                W_old = W_new;
                E1_old = E1_new;
                E2_old = E2_new;
                P1 = P2;
                //cout<<-beta.Energy<<endl;
            }
            else P2 = P1;

            energy1 += E1_old;
            energy2 -= E2_old;
        }//MCS

    }//EQMC
    energy1 /= 1.0*MCS;
    energy2 /= 1.0*MCS;
    //energy2 += beta.num/4.0;

    cout<<0.5*energy1+ beta.numLattB/4.0<<endl;
    cout<<-0.5*energy1/beta.numSpin - 0.5<<endl;
    cout<<0.5*energy2/beta.numSpin + 0.5<<endl;


  return 0;
};