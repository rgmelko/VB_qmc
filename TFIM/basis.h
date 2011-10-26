#ifndef BASIS_H
#define BASIS_H

#include "head_proj.h"
#include "simparam.h"


class Basis: public PARAMS
{

    public:
      vector <index2> OperatorList; //The operator list of 2m
      // (-2,i): an off-diagonal site operator h(sigma^+_i + sigma^-_i)
      // (-1,i): a diagonal site operator h
      // (i,j):  a diagonal bond operator J(sigma^z_i sigma^z_j + 1)

      vector<int> S_left; //the left and right trial spin state
      vector<int> S_right;

      Basis(MTRand &); //constructor
      void DiagonalUpdate(MTRand &);
      void printBasis();

};


Basis::Basis(MTRand& ran){//constructor

    for (int i=0; i< numSpin; i++){
        S_left.push_back(0);   //start with all spin parallel
        S_right.push_back(0);
    }

    int bond;
    index2 temp;
    for (int i=0; i<2*m_; i++){

        if (ran.randInt(1) == 0) { //flip a coin
            temp.set(-1,ran.randInt(numSpin-1)); //diagonal site operator
        }
        else {
            bond = ran.randInt(numLattB-1);      //diagonal bond operator
            temp = Bst[bond]; //checked overloaded =
            cout<<"temp "<<bond<<" ";
            Bst[bond].print();
            cout<<endl;

        }

        OperatorList.push_back(temp);

    }//i

}//----------------constructor



//----------------DiagonalUpdate function
void Basis::DiagonalUpdate(MTRand& ran){

    vector<int> S_prop; //this is the temporary propagated spin state
    S_prop = S_left; //assign to the left spin state

    //Probability for a single-site diagonal operator
    double hProb = h_x/(1.0+h_x);


    int bond;
    int flag;
    for(int i=0; i<OperatorList.size(); i++){

        if (OperatorList[i].A == -2) //this is a off-diagonal site operator
            S_prop[OperatorList[i].B] = S_prop[OperatorList[i].B]^1 ; //binary spin flip

        else { //sample a new diagonal operator

            flag = 0;
            do{ //repeat until a valid selection is made
                if (hProb > ran.rand() ){ //probability to choose a single-site h operator
                    OperatorList[i].set(-1,ran.randInt(numSpin-1)); //diagonal site operator
                    flag = 1; //successful!
                }
                else{
                    cout<<"B "<<i<<" "<<endl;
                    bond = ran.randInt(numLattB-1);  //new bond for diagonal bond operator
                    if (S_prop[Bst[bond].A] == S_prop[Bst[bond].B]) {
                        //spins are the same on the new bond
                        OperatorList[i] = Bst[bond]; //check overloaded =
                        flag = 1; //successful!
                    }
                } //else the diagonal operator stays unchanged!  repeat
            }while(flag == 0);
        }//insert diagonal

    }//i  the 2*m propagation

    //DEBUG: check if the state was propagated correctly
    for (int i=0; i<S_prop.size(); i++)
        if (S_prop[i] != S_right[i]) cout<<"Basis state prop error: DIAG UPDATE \n";


}//----------------DiagonalUpdate




//----------------print function
void Basis::printBasis(){
    cout<<"Basis "<<endl;
    for (int i=0; i<S_left.size(); i++){
        cout<<S_left[i]<<" ";
    }
    cout<<endl;
    for (int i=0; i<S_right.size(); i++){
        cout<<S_left[i]<<" ";
    }
    cout<<endl;
    for (int i=0; i<OperatorList.size(); i++){
        OperatorList[i].print();
    }
    cout<<endl;
}//print



#endif
