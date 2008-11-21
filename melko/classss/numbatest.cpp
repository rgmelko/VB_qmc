
#include<iostream>
#include "bignumba.h"

using namespace std;

main()
{
  bignum first(8888,888888888);
  bignum second(4444,444444444);
 
  bignum answer;

  answer = first + second;
 
  cout << first.left << first.right << " + " <<
    second.left << second.right << " = " <<
    answer.left << answer.right << endl;

  return 0;
}
