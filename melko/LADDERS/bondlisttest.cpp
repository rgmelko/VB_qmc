#include<iostream>
#include<iomanip>
using namespace std;

int main()
{
  int legs, length;

  cout << "# legs?" << endl;
  cin >> legs;
  cout << "length?" << endl;
  cin >> length;

  const int number_of_nnbonds = 2*legs*length - legs - length;
  int nnbonds[number_of_nnbonds][2];
  int bondnum=0;

  
  //the first bonds are of the form (0,1),(1,2),(2,3) etc
  for(bondnum; bondnum < legs*length-1; bondnum++)
    {
      nnbonds[bondnum][0] = bondnum;
      nnbonds[bondnum][1] = nnbonds[bondnum][0]+1;
    }

  //the rest are more complicated
  /* 
    

   */
  while(bondnum < number_of_nnbonds)
    {
      int sitenum = 0;
      for(int legcounter=legs*2-1; legcounter>1; legcounter-=2)
	{
	  for(int i00 = sitenum; i00<(length-1)*legs; i00+=legs)
	    {
	      nnbonds[bondnum][0]= i00;
	      nnbonds[bondnum][1]= i00+legcounter;
	      bondnum++;
	    }
	  sitenum++;
	}
    }
  

  cout << endl << "The nn bonds are:\n\n" ;

  for(int i01 = 0; i01<number_of_nnbonds; i01++)
    {
      cout << setw(3) << i01 << ")"<< setw(3) << nnbonds[i01][0] << " ,"
	   << setw(2) << nnbonds[i01][1] << endl;
    }

  cout << endl;
  
  return 0;
}
