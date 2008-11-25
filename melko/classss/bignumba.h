#ifndef bignumba
#define bignumba

class bignum {
 public:
  bignum();
  bignum(long unsigned int lhs, long unsigned int rhs);
  long unsigned int left, right;
  
  bool operator> (bignum);
  bool operator>= (bignum);
  bool operator< (bignum);
  bool operator<= (bignum);
  bool operator== (bignum);
  bignum operator+ (bignum);
  bignum operator+ (int);
  bignum operator++(int);
  bignum operator% (bignum);
  double doublify();
  //  string print();



};

bignum::bignum(){
  left=0, right=0;
}

bignum::bignum(long unsigned int lhs, long unsigned int rhs){
  left=lhs, right=rhs;
}

bool bignum::operator> (bignum bigcompare) {
  if(left > bigcompare.left){return 1;}
  if(left < bigcompare.left){return 0;}
  else{
    if(right > bigcompare.right){return 1;}
    else{return 0;}
  }
}
bool bignum::operator>= (bignum bigcompare) {
  if(left > bigcompare.left){return 1;}
  if(left < bigcompare.left){return 0;}
  else{
    if(right >= bigcompare.right){return 1;}
    else{return 0;}
  }
}

bool bignum::operator< (bignum bigcompare) {
  if(left < bigcompare.left){return 1;}
  if(left > bigcompare.left){return 0;}
  else{
    if(right < bigcompare.right){return 1;}
    else{return 0;}
  }
}
bool bignum::operator<= (bignum bigcompare) {
  if(left < bigcompare.left){return 1;}
  if(left > bigcompare.left){return 0;}
  else{
    if(right <= bigcompare.right){return 1;}
    else{return 0;}
  }
}

bool bignum::operator== (bignum bigcompare) {
  if( (left == bigcompare.left)&(right == bigcompare.right) ) {
    return 1;
  }
  else{return 0;}
}


bignum bignum::operator+ (bignum bigadd) {
  bignum temp;
  long unsigned int runover;
  runover = (right + bigadd.right);
  if(runover>1000000000){runover = runover/1000000000;}
  else{runover = 0;}
  temp.right  = right + bigadd.right - runover*1000000000;
  temp.left = left + bigadd.left + runover;
  return(temp);
}

bignum bignum::operator+ (int smalladd) {
  bignum temp;
  long unsigned int runover;
  runover = (right + smalladd);
  if(runover>1000000000){runover = runover/1000000000;}
  else{runover = 0;}
  temp.right  = right + smalladd - runover*1000000000;
  temp.left = left + runover;
  return(temp);
}

bignum bignum::operator++ (int) {
  bignum temp;
  long unsigned int runover;
  runover = (right + 1);
  if(runover>1000000000){runover = runover/1000000000;}
  else{runover = 0;}
  temp.right  = right + 1 - runover*1000000000;
  temp.left = left + runover;
  return(temp);
}

bignum bignum::operator% (bignum bigmod) {
  bignum temp;
  

}

double bignum::doublify () {
  double temp;
  temp = 1000000000.0*left + right*1.0;
  return(temp);
}

//char bignum::print () {
//  string temp;
//}  

#endif
