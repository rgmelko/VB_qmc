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

#endif
