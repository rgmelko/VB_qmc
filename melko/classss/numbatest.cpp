
#include<iostream>
#include "bignumba.h"

using namespace std;

main()
{
  bignum ones(1,111111111);
  bignum twos(2,222222222);
 
  bignum threes;

  threes = ones + twos;

  cout << threes.left << threes.right << endl;

  return 0;
}
