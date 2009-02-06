#ifndef ladder
#define ladder

class LADDER 
{

 public:
  int legs, length, number_of_bondops;
  const int number_of_nnbonds = 2*legs*length - legs - length;
  const int number_of_bonds = legs*length/2;
  int nnbonds[number_of_nnbonds][2];
  int bonds[number_of_bonds][2];
  bool nn[number_of_bonds];
  
  LADDER(const int legs,const int length,const in number_of_bondops);
  
  void nnbondlist();
  void apply_ops(int bondops[number_of_bondops]);
  
}
#endif


LADDER::LADDER(const int a, const int b, const int c)
{
  legs = a;
  length = b;
  number_of_bondops = c;
}

LADDER::nnbondlist()
{
  int bondnum = 0;
  
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
}
