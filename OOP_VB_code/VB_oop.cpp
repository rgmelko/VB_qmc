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
    int MCS = 200000;
    double E1_new, E1_old;
    double energy1;
    for (int EQMC = 0; EQMC <2; EQMC++) { //EQL and MCS run loop
        MCS *= 2;
        P2 = P1;  //set projectors equal
        energy1 = 0;
		alpha.Propogate(P1,beta);
        //cout<<(alpha|beta)<<endl; //calculate overlap

		W_old =  beta.Weight;

		E1_new = beta.Calc_Energy()-1; //energy calculation from total nn bond 

        for (int i=0; i<MCS; i++){

            P2.Sample_Ops(mrand);
            alpha.Propogate(P2,beta);
            W_new =  beta.Weight;
            DeltaW = pow(2,W_old - W_new);

            E1_new = beta.Calc_Energy()-1; //energy

            if (DeltaW > mrand.rand()){ //Accept the move
                W_old = W_new;
                E1_old = E1_new;
                P1 = P2;
                //cout<<-beta.Energy<<endl;
            }
            else P2 = P1;  //reject the move

            energy1 += E1_old;

        }//MCS

    }//EQMC
    energy1 /= 1.0*MCS;
    //energy2 += beta.num/4.0;

    cout<<-0.5*energy1/beta.numSpin - 0.5<<endl;


  return 0;
};