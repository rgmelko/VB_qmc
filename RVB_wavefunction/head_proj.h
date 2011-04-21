#ifndef HEAD_PROJ_H
#define HEAD_PROJ_H

#include <iostream>
using namespace std;

#include <iomanip>
#include <vector>
#include "MersenneTwister.h"
#include "matrix.h"

class index2{

    public:

        int A;
        int B;

        index2() {A=0; B=0;};
        index2(const int a, const int b) {A=a; B=b;};
        void set(const int a, const int b) {A=a; B=b;};
        void print(){ cout<<"("<<A<<","<<B<<")"<<endl; };
};


class index4{

    public:

        int A;
        int B;
        int C;
        int D;

        index4() {A=0; B=0; C=0; D=0;};
        index4(const int a, const int b, const int c, const int d) 
              {A=a; B=b; C=c; D=d;};
        void set (const int a, const int b, const int c, const int d) 
              {A=a; B=b; C=c; D=d;};
        void print(){ cout<<"("<<A<<","<<B<<","<<C<<","<<D<<")"<<endl; };
};


#endif
