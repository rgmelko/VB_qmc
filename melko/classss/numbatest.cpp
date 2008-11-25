
#include<iostream>
#include "bignumba.h"

using namespace std;

main()
{
  bignum first(8,888888888);
  bignum second(4444,444444444);
  long int hat = 2000000000;
  double sectron = 4444444444444.0;
  bignum answer;
  double hatt;


  answer = second;
  cout.precision(8);
  answer = answer + 1;

  hatt = first.doublify();
  
 
  cout << answer.left << answer.right << endl;
  cout << sectron << endl;
  cout << first.left << first.right<< endl;
  cout << hatt << endl;
  return 0;
}
