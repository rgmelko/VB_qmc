#ifndef BASIS_H
#define BASIS_H

#include "head_proj.h"
#include "simparam.h"


class Basis: public PARAMS
{
    private:
       vector<int> LinkList;
       vector<int> LegType;
       vector<index3> Associates;

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
      void ClusterUpdate(MTRand &);
      void printBasis();
      void printLinkedList();

};


Basis::Basis(MTRand& ran){//constructor

    for (int i=0; i< numSpin; i++){
        S_left.push_back(1);   //start with all spin parallel
        S_right.push_back(1);
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
            //cout<<"temp "<<bond<<" ";
            //Bst[bond].print();
            //cout<<endl;
        }

        OperatorList.push_back(temp);

    }//i
     
    //S_right.at(S_right.size()-1) = S_right.at(S_right.size()-1)^1; 
    //temp.set(-2,numSpin-1);
    //OperatorList.at(m_) =  temp;

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
                    //cout<<"B "<<i<<" "<<endl;
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

    index3 empty(-1,-1,-1);
    index3 temp3;

    vector<int> First;
    for (int i=0; i<numSpin; i++){ //the first vertex leg for each spin
        First.push_back(i);
        //below, build the first N vertices from the left-basis
        LinkList.push_back(-99); //unknown what these link to!
        LegType.push_back(S_left[i]); //0 or 1
        Associates.push_back(empty); //these have no associates
    }

    vector<int> S_prop; //this is the temporary propagated spin state
    S_prop = S_left; //assign to the left spin state

    int count = numSpin;  
    int site, site1, site2;
    //The linked list is now size N.  Add the 2m operators each of 4 or 2 legs
    for(int i=0; i<OperatorList.size(); i++){

        //assign non-trivial associates
        if (OperatorList[i].A != -2 && OperatorList[i].A != -1){
            temp3.set(count+1,count+2,count+3); Associates.push_back(temp3);
            temp3.set(count,count+2,count+3); Associates.push_back(temp3);
            temp3.set(count,count+1,count+3); Associates.push_back(temp3);
            temp3.set(count,count+1,count+2); Associates.push_back(temp3);
            count += 4;
        }//done assigning associates
        else{
            Associates.push_back(empty);
            Associates.push_back(empty);
            count += 2;
        }

        if (OperatorList[i].A == -2){ //1-site off-diagonal operator is encountered
            site = OperatorList[i].B;
            //"lower" or leftmost leg
            LinkList.push_back(First[site]); //site index
            LegType.push_back(S_prop[site]); //the spin of the leg
            S_prop[site] = S_prop[site]^1;   //this is off-d: flip it
            LinkList[First[site]] = LinkList.size()-1; //this leg links backwards...
            First[site] = LinkList.size(); //update
            //"upper" or rightmost leg
            LinkList.push_back(-99); //null site index
            LegType.push_back(S_prop[site]); //the spin of the leg (flipped)
        }
        else if (OperatorList[i].A == -1){ //1-site diagonal operator is encountered
            site = OperatorList[i].B;
            //"lower" or leftmost leg
            LinkList.push_back(First[site]); //site index
            LegType.push_back(S_prop[site]); //the spin of the leg
            LinkList[First[site]] = LinkList.size()-1; //this leg links backwards...
            First[site] = LinkList.size(); //update
            //"upper" or rightmost leg
            LinkList.push_back(-99); //null site index
            LegType.push_back(S_prop[site]); //the spin of the leg 
        }
        else {//2-site diagonal operator is encountered (4 legs)
            //lower left
            site1 = OperatorList[i].A;
            LinkList.push_back(First[site1]); //site index
            LegType.push_back(S_prop[site1]); //the spin of the leg
            LinkList[First[site1]] = LinkList.size()-1; //this leg links backwards...
            First[site1] = LinkList.size()+1;
            //lower right
            site2 = OperatorList[i].B;
            LinkList.push_back(First[site2]); //site index
            LegType.push_back(S_prop[site2]); //the spin of the leg
            LinkList[First[site2]] = LinkList.size()-1; //this leg links backwards...
            First[site2] = LinkList.size()+1;
            //upper left
            LinkList.push_back(-99); //null site index
            LegType.push_back(S_prop[site1]); //the spin of the leg 
            //upper right
            LinkList.push_back(-99); //null site index
            LegType.push_back(S_prop[site2]); //the spin of the leg 
        }

    }//i

    //now add the legs of the final ("top"or right-hand) spin state
    for (int i=0; i<numSpin; i++){ 
        LinkList.push_back(First[i]);
        LinkList[First[i]] = LinkList.size()-1;
        LegType.push_back(S_prop[i]); //0 or 1
        Associates.push_back(empty);
    }

    cout<<"Ass size :"<<Associates.size()<<endl;
    cout<<"LL size :"<<LinkList.size()<<endl;
    cout<<"LT size :"<<LegType.size()<<endl;

    //DEBUG: check if the state was propagated correctly
    for (int i=0; i<S_prop.size(); i++)
        if (S_prop[i] != S_right[i]) cout<<"Basis state prop error: LINKED LIST\n";

}//----------------LinkedList



//----------------ClusterUpdate
void Basis::ClusterUpdate(MTRand& ran){

    vector<int> inCluster(LinkList.size(), 0); //nothing in clusters yet
    stack<int> cluster;

    int leg, assoc;
    bool flip;

    int ccount = 0;
    for (int i=0; i<LinkList.size(); i++){ //loop to find all clusters

        //add a new leg
        if (inCluster[i] == 0 && Associates[i].A == -1){ //spins and site ops only
            ccount ++; //cluster counter

            cluster.push(i);
            inCluster[cluster.top()] = ccount;

            if (ran.rand() < 0.5) flip = true; else flip = false; //flip a coin for SW

            if (flip == true) LegType[cluster.top()] = LegType[cluster.top()]^1;

            while(!cluster.empty()){ //build the cluster associated with this leg

                //first follow the link
                leg = LinkList[cluster.top()];
                cluster.pop();

                if (inCluster[leg] == 0){
                    inCluster[leg] = ccount; //add the linked leg
                    if (flip == true) LegType[leg] = LegType[leg]^1;
                    //now check all associates
                    assoc = Associates[leg].A;
                    if (assoc != -1) { 
                        cluster.push(assoc); inCluster[assoc] = ccount; 
                        if (flip == true) LegType[assoc] = LegType[assoc]^1;
                        assoc = Associates[leg].B;
                        cluster.push(assoc); inCluster[assoc] = ccount; 
                        if (flip == true) LegType[assoc] = LegType[assoc]^1;
                        assoc = Associates[leg].C;
                        cluster.push(assoc); inCluster[assoc] = ccount; 
                        if (flip == true) LegType[assoc] = LegType[assoc]^1;
                    }
                }

            }//while

        }//if building a new cluster

    }//i

    for (int i=0; i<inCluster.size(); i++)
        cout<<inCluster[i]<<" ";
    cout<<endl;

    //map back basis states and operator list
    for(int i=0; i<numSpin; i++){
        S_left[i] = LegType[i];
    }//i
    int count = numSpin;  
    for(int i=0; i<OperatorList.size(); i++){

        //assign non-trivial associates
        if (OperatorList[i].A != -2 && OperatorList[i].A != -1){
            count += 4;
        }//done assigning associates
        else{
            if (LegType[count] == LegType[count+1])
                OperatorList[i].A = -1;
            else
                OperatorList[i].A = -2;
            count += 2;
        }
    }//i
    for(int i=0; i<numSpin; i++){
        S_right[i] = LegType[LegType.size()-numSpin + i];
    }//i


}//----------------ClusterUpdate


//----------------print LinkedList
void Basis::printLinkedList(){

    for (int i=0; i<LinkList.size(); i++){ 
        cout<<i<<" ";
        cout<<LinkList[i]<<" ";
        cout<<LegType[i]<<"\n";
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
        cout<<S_right[i]<<" ";
    }
    cout<<endl;
    for (int i=0; i<OperatorList.size(); i++){
        OperatorList[i].print();
    }
//    cout<<endl;
//    for (int i=0; i<Associates.size(); i++){
//        Associates[i].print();
//    }
//    cout<<endl;

}//print



#endif
