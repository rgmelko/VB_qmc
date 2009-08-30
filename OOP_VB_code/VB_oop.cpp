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

    //MTran met_rand(2343983); //random number for metropolis

    Projector P1; //initialize projector 1
    //P1.print();
    Projector P2; //initialize projector 2
    //P2.print();

    Basis alpha;
    //alpha.print();
  
    double W_old, W_new; //old and new weights
    P2 = P1;  //set projectors equal
    
    Basis beta;
    alpha.Propogate(P1,beta);
    W_old =  beta.Weight;

//    for (int i=0; i<10000; i++){
//        P2.Sample_Ops();
//        alpha.Propogate(P2,beta);
//        W_new =  beta.Weight;
//
//        if (W_new > W_old){//keep changes
//            W_old = W_new;
//            P1 = P2;
//        }
//        else if (W_new/W_old > met_rand.rand();){
//            W_old = W_new;
//            P1 = P2;
//        }
//        else P2 = P1;
//    }//MCS



  return 0;
};