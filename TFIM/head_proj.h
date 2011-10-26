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
        index2 operator=(const index2 &z) {A=z.A; B=z.B; return *this;};
        void set(const int a, const int b) {A=a; B=b;};
        void print(){ cout<<"("<<A<<","<<B<<")"<<endl; };
};


#endif
