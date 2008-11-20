#ifndef bignumba
#define bignumba

class bignum {
 public:
  bignum();
  bignum(long unsigned int lhs, long unsigned int rhs);
  long unsigned int left, right;
  
  bignum operator+ (bignum);


};

bignum::bignum(){
  left=0, right=0;
}

bignum::bignum(long unsigned int lhs, long unsigned int rhs){
  left=lhs, right=rhs;
}
  
bignum bignum::operator+ (bignum bigadd) {
  bignum temp;
  long unsigned int runover = (right + bigadd.right)%1000000000;
  temp.right  = right + bigadd.right - runover*1000000000;
  temp.left = left + bigadd.left + runover;
  return(temp);
}

#endif
