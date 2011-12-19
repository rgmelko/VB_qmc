#ifndef MEASURE_H
#define MEASURE_H

#include "head_proj.h"
#include "basis.h"
#include <fstream>
#include <cmath>

class Measure: public PARAMS
{
    private: 
      //observables here
      double Energy;
	  double Mag1;
	  double Mag2;
	  double Mag3;
      vector<double> Renyi; //naiive direct estimator (obsolete)
      vector<double> Renyi2;//improved LR cluster estimator

      //This is the number of clusters of the unswapped simulation
      int ClustNumber;

    public:

      Measure();
      void zero();
      void measure_E(const Basis &);
      void measure_M(const Basis &, const int &);
      void measure_M_mod(const vector<int>&, const vector<int>&);
      void Renyi_direct(const int& , const int& , const int& );
      void Renyi_LRclust(const vector<int>& ,const vector<int>&, const vector<int>&);
      vector<int> LRoverlap(const vector<int>&, const vector<int>&, int &);
      void output();
  
};

Measure::Measure() {//constructor
    Energy = 0.0; 
    Mag1 = 0.0; 
    Mag2 = 0.0;
    Mag3 = 0.0;
    Renyi.assign(numSpin/2-1,0.0);
    Renyi2.assign(numSpin/2-1,0.0);
};

void Measure::zero(){

    Energy = 0.0;
    Mag1 = 0.0;
    Mag2 = 0.0;
    Mag3 = 0.0;
    Renyi.assign(numSpin/2-1,0.0);
    Renyi2.assign(numSpin/2-1,0.0);

}//zero

void Measure::Renyi_direct(const int& index, const int& numer, const int& denom){

    int frac_s = numer-denom;

    Renyi[index] += pow(2.0,frac_s);

}//Renyi_direct


//This function calculates the swap operator directly from the Left and Right
// "clusters" that have been calculated in the linked list
//*** NOTE always calculate measure_M_mod *FIRST*
void Measure::Renyi_LRclust(const vector<int>& inA,
                            const vector<int>& Left, const vector<int>& Right){

    vector<int> Mtemp;
    int max_index; //this is the maximum cluster index in the overlap

    int numRealSpin = numSpin/2;

    int frac_s, numer;
    int denom = ClustNumber; //global variable calculated in measure_M_mod

    vector<int> RightSwap(Right); //copy constructor?
    int temp;
    int renyi_index; //1D solution only
    //---Swap the Right projector here
    for (int i=0; i<inA.size(); i++){
        if (inA[i] != 0){
            renyi_index = i; 
            temp = RightSwap[i];
            RightSwap[i] = RightSwap[i+numRealSpin];
            RightSwap[i+numRealSpin] = temp;
        }//inA
    }//i

    Mtemp = LRoverlap(Left,RightSwap,max_index);

    vector<int> MidClustsNum(max_index+1,0);
    for (int k=0; k<Mtemp.size(); k++)
        MidClustsNum[Mtemp[k]] = 1;

    int counter = 0; //number of clusters
    for (int k=0; k<MidClustsNum.size(); k++)
        counter += MidClustsNum[k];

    numer = counter;
    frac_s = numer-denom;
    Renyi2.at(renyi_index) += pow(2.0,frac_s);

    
}//Renyi_LRclust


void Measure::measure_E(const Basis & basis){

    int n_0=0;

    for(int i=0; i<basis.OperatorList.size(); i++){

        if (basis.OperatorList[i].A == -1)
            n_0++;
    }//i

    Energy += 1.0*n_0;
    //cout<<Energy<<endl;

}//measure_E


void Measure::measure_M(const Basis & basis, const int & L2){

	int m_0 = 0;

    vector<int> S_prop; 
    S_prop = basis.S_left; 

    for(int i=0; i<basis.OperatorList.size()/2; i++){

        if (basis.OperatorList[i].A == -2) //this is a off-diagonal site operator
            S_prop[basis.OperatorList[i].B] = S_prop[basis.OperatorList[i].B]^1 ; //spin flip

	}

	for (int i=0; i<basis.numSpin; i++)
		m_0 += (2*S_prop[i]-1);


    Mag1 += 1.0*m_0*m_0; //m^2

    Mag2 += 1.0*L2;

}//measure_M

void Measure::measure_M_mod(const vector<int>& Left, const vector<int>& Right){

    vector<int> Mtemp;
    int max_index; //this is the maximum cluster index in the overlap
    Mtemp = LRoverlap(Left,Right,max_index);

    vector<int> MidClustsNum(max_index+1,0);
    vector<int> MidClustsSize(max_index+1,0);
    for (int k=0; k<Mtemp.size(); k++){
        MidClustsNum[Mtemp[k]] = 1;
        MidClustsSize[Mtemp[k]] += 1;
    }

    int counter = 0; //number of clusters
    int sizesquared = 0;
    for (int k=0; k<MidClustsNum.size(); k++){
        counter += MidClustsNum[k];
        sizesquared += MidClustsSize[k]*MidClustsSize[k];
    }

    //cout<<"new clust #: "<<counter<<endl;
    //cout<<counter<<endl;

    ClustNumber = counter; //private global variable

    Mag3 += 1.0*sizesquared;

}//measure_M_mod


//a function that takes the L and R "center" cluster vectors and calculates
//the overlap vector 
vector<int> Measure::LRoverlap(const vector<int>& Left, const vector<int>& Right, int & max){

    int Nspin = Left.size();

    vector<int> Mtemp;

    stack<int> Rstack;
    stack<int> Lstack;
    int current;
    //cout<<Nspin<<endl;
    Mtemp.assign(Nspin,0);
    bool keepgoing = true;
    //for (int ii=1; ii<=Nspin; ii++){
    int ii = 0;
    while(keepgoing == true){
        ii++;
        keepgoing = false;

        current = ii;

        for (int j=0; j<Nspin; j++){
            if (Left[j] == current && Mtemp[j] == 0){
                Mtemp[j] = ii;
                Rstack.push(Right[j]);
            }
            if (Mtemp[j] == 0) keepgoing = true;
        }

        do{
            while(!Rstack.empty()){ 
                current = Rstack.top();
                Rstack.pop();
                for (int j=0; j<Nspin; j++)
                    if (Right[j] == current && Mtemp[j] == 0){
                        Mtemp[j] = ii;
                        Lstack.push(Left[j]);
                    }
            }//while Rstack

            while(!Lstack.empty()){ 
                current = Lstack.top();
                Lstack.pop();
                for (int j=0; j<Nspin; j++)
                    if (Left[j] == current && Mtemp[j] == 0){
                        Mtemp[j] = ii;
                        Rstack.push(Right[j]);
                    }
            }//while Rstack

        }while(!Lstack.empty() || !Rstack.empty());

        //cout<<ii<<": ";
        //for (int k=0; k<Mtemp.size(); k++)
        //    cout<<Mtemp[k]<<" ";
        //cout<<endl;

    }//ii

    //for (int k=0; k<Mtemp.size(); k++)
    //    cout<<Mtemp[k]<<" ";
    //cout<<endl;

    max = ii; //this is the maximum cluster index contained in the overlap
    return Mtemp;

}//LRoverlap


void Measure::output(){

    double one_over_n;
	ofstream cfout;
	cfout.open("00.data",ios::app);
    cfout<<setprecision(8);

    one_over_n = Energy/(1.0*MCS_);
    cfout<<numSpin*h_x*2.0*m_ / one_over_n<<" ";
    cfout<<-(numSpin*h_x*2.0*m_ / one_over_n - numSpin*h_x - numLattB)/numSpin<<" ";
    //cfout<<-(-1.0-h_x + h_x*2.0*m_ / one_over_n)<<" "; //wrong for OBC?
	cfout<<Mag1/(1.0*MCS_*1.0*numSpin*numSpin)<<" ";
	cfout<<Mag2/(1.0*MCS_*1.0*numSpin*numSpin)<<" ";
	cfout<<Mag3/(1.0*MCS_*1.0*numSpin*numSpin);
    cfout<<endl;

	cfout.close();

	cfout.open("01.data",ios::app);

    for (int i=0; i<Renyi.size(); i++){
        cfout<<i<<" "<<-log(Renyi[i]/(1.0*MCS_))<<" ";
        cfout<<i<<" "<<-log(Renyi2[i]/(1.0*MCS_))<<endl;
    }
    //cout<<endl;

	cfout.close();

}//output


#endif
