//
// VB_oop.cpp: the main body of the program
// A OOP VB basis projector code 
// Roger Melko, Aug 29 2009, Santa Barbara
//
// Double projector algorithm: October 20, ENS Lyon
//
#include "head_proj.h"
#include "simparam.h"
#include "basis.h"
#include "projector.h"
#include "measure.h"

int main(){

    PARAMS param;
    MTRand mrand(param.SEED_); //random number for metropolis

    Projector P1(mrand); //initialize projector 1
    Projector P2(mrand); //initialize projector 2
    Projector Pold1(mrand); //initialize old projector  (for metropolis)
    Projector Pold2(mrand); //initialize old projector  (for metropolis)

    Basis alpha;  //This is the starting baisis |V>
    Basis beta_1; //This is P1|V> = w1|V1>
    Basis beta_2; //This is P2|V> = w2|V2>

    long int W1_old, W1_new; //old and new weights
    long int W2_old, W2_new; //old and new weights
    long int N_loop_old, N_loop_new;
    double DeltaW;
    
    Measure Observ; //create measurement object

    int MCsteps;
    for (int EQMC = 0; EQMC <2; EQMC++) { //EQL and MCS run loop
        Observ.zero(); //set observable values to zero

        //set to old values
        P1 = Pold1;
        P2 = Pold2;

		alpha.Propogate(P1,beta_1);  //P1|alpha> = W1|beta_1>
		alpha.Propogate(P2,beta_2);  //P1|alpha> = W1|beta_1>
        
        W1_old = beta_1.Weight;
        W2_old = beta_2.Weight;
        N_loop_old = beta_1|beta_2;     // calculate number of loops in <V1 | V2>

        //initialize measurements: two steps
        //Observ.measure_energy(beta_1, beta_2); //make initial measurements (assign "new" values)
        Observ.measure_energy2(beta_1, beta_2); //make initial measurements (assign "new" values)
        //Observ.measure_CL2L2(beta_1, beta_2); 

        if (EQMC == 0) MCsteps = param.EQL_;
        else MCsteps = param.MCS_;

        for (int i=0; i<MCsteps; i++){

            //-----sample projector 1 first---------------------------
            P1.Sample_Ops(mrand);       //sample new operators
            alpha.Propogate(P1,beta_1); //propogate basis
            W1_new =  beta_1.Weight;    //calculate new weight
            N_loop_new = beta_1|beta_2; //calcualte new overlap
            DeltaW = pow(2,W1_old - W1_new + N_loop_new - N_loop_old);
            if (DeltaW > mrand.rand()){ //Accept the move
                W1_old = W1_new;
                Pold1 = P1;
                N_loop_old = N_loop_new;
                //measurements            
                Observ.measure_energy2(beta_1, beta_2); //measure energy
				//Observ.measure_energy(beta_1, beta_2); //measure energy
                //Observ.measure_CL2L2(beta_1, beta_2);  //measure spin-spin correlation function
            }
            else {  //reject the move          
                P1 = Pold1;                    
                N_loop_new = N_loop_old;       
            }                                  
            //--------------------------end sample proj 1 --------------

            Observ.record(); //assign running total

            //-----sample projector 2 first---------------------------
            P2.Sample_Ops(mrand);       //sample new operators
            alpha.Propogate(P2,beta_2); //propogate basis
            W2_new =  beta_2.Weight;    //calculate new weight
            N_loop_new = beta_1|beta_2;         //calcualte new overlap
            DeltaW = pow(2,W2_old - W2_new + N_loop_new - N_loop_old);
            if (DeltaW > mrand.rand()){ //Accept the move
                W2_old = W2_new;
                Pold2 = P2;
                N_loop_old = N_loop_new;
                //measurements
			    Observ.measure_energy2(beta_1, beta_2); //measure energy
			    //Observ.measure_energy(beta_1, beta_2); //measure energy
                //Observ.measure_CL2L2(beta_1, beta_2);  //measure spin-spin correlation function
            }
            else {  //reject the move          
                P2 = Pold2;                    
                N_loop_new = N_loop_old;       
            }                                  
            //--------------------------end sample proj 2 --------------

            Observ.record(); //assign running total

        }//MCS

    }//EQMC


    Observ.output(param); //output observables


  return 0;
};