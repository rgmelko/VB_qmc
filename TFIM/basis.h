#ifndef BASIS_H
#define BASIS_H

#include "head_proj.h"
#include "simparam.h"


class Basis: public PARAMS
{
    private:
       vector<int> LinkList;
       vector<int> LinkLegType;

    public:
      vector <index2> OperatorList; //The operator list of 2m
      // (-2,i): an off-diagonal site operator h(sigma^+_i + sigma^-_i)
      // (-1,i): a diagonal site operator h
      // (i,j):  a diagonal bond operator J(sigma^z_i sigma^z_j + 1)

      vector<int> S_left; //the left and right trial spin state
      vector<int> S_right;

      Basis(MTRand &); //constructor
      void DiagonalUpdate(MTRand &);
      void LinkedList();
      void printBasis();
      void printLinkedList();

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



//----------------LinkedList function
void Basis::LinkedList(){

    vector<int> First;
    for (int i=0; i<numSpin; i++){ //the first vertex leg for each spin
        First.push_back(i);
        //below, build the first N vertices from the left-basis
        LinkList.push_back(-99); //unknown what these link to!
        LinkLegType.push_back(S_left[i]); //0 or 1
    }

    vector<int> S_prop; //this is the temporary propagated spin state
    S_prop = S_left; //assign to the left spin state

    int site, site1, site2;
    //The linked list is now size N.  Add the 2m operators each of 4 or 2 legs
    for(int i=0; i<OperatorList.size(); i++){

//        if (OperatorList[i].A == -2){ //1-site off-diagonal operator is encountered
//            site = OperatorList[i].B;
//            //"lower" or leftmost leg
//            LinkList.push_back(First[site]); //site index
//            LinkLegType.push_back(S_prop[site]); //the spin of the leg
//            S_prop[site] = S_prop[site]^1;   //this is off-d: flip it
//            LinkList[First[site]] = LinkList.size()-1; //this leg links backwards...
//            First[site] = LinkList.size(); //update
//            //"upper" or rightmost leg
//            LinkList.push_back(-99); //null site index
//            LinkLegType.push_back(S_prop[site]); //the spin of the leg (flipped)
//        }
        if (OperatorList[i].A == -1){ //1-site diagonal operator is encountered
            site = OperatorList[i].B;
            //"lower" or leftmost leg
            LinkList.push_back(First[site]); //site index
            LinkLegType.push_back(S_prop[site]); //the spin of the leg
            LinkList[First[site]] = LinkList.size()-1; //this leg links backwards...
            First[site] = LinkList.size(); //update
            //"upper" or rightmost leg
            LinkList.push_back(-99); //null site index
            LinkLegType.push_back(S_prop[site]); //the spin of the leg 
        }
        else {//2-site diagonal operator is encountered (4 legs)
            //lower left
            site1 = OperatorList[i].A;
            LinkList.push_back(First[site1]); //site index
            LinkLegType.push_back(S_prop[site1]); //the spin of the leg
            LinkList[First[site1]] = LinkList.size()-1; //this leg links backwards...
            First[site1] = LinkList.size()+1;
            //lower right
            site2 = OperatorList[i].B;
            LinkList.push_back(First[site2]); //site index
            LinkLegType.push_back(S_prop[site2]); //the spin of the leg
            LinkList[First[site2]] = LinkList.size()-1; //this leg links backwards...
            First[site2] = LinkList.size()+1;
            //upper left
            LinkList.push_back(-99); //null site index
            LinkLegType.push_back(S_prop[site1]); //the spin of the leg 
            //upper right
            LinkList.push_back(-99); //null site index
            LinkLegType.push_back(S_prop[site2]); //the spin of the leg 
        }

    }//i

    //now add the legs of the final ("top"or right-hand) spin state
    for (int i=0; i<numSpin; i++){ 
        LinkList.push_back(First[i]);
        //if (First[i] == i) //there have been no operators!
        //    LinkList[i] = i;
        LinkList[First[i]] = LinkList.size()-1;
        LinkLegType.push_back(S_prop[i]); //0 or 1
    }


}//----------------LinkedList

//----------------print LinkedList
void Basis::printLinkedList(){

    for (int i=0; i<LinkList.size(); i++){ 
        cout<<i<<" ";
        cout<<LinkList[i]<<" ";
        cout<<LinkLegType[i]<<"\n";
    }


}



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
