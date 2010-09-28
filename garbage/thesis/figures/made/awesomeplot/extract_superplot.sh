#!/bin/sh
awk '{if($1==2) print "2   ",  $2  }'  3leg*
awk '{if($1==3) print "3   ",  $2  }'  3leg*
awk '{if($1==4) print "4   ",  $2  }'  4leg* 
awk '{if($1==5) print "5   ", $2  }'   5leg* 
awk '{if($1==6) print "6   ", $2  }'   6leg* 
awk '{if($1==7) print "7   ", $2  }'   7leg* 
awk '{if($1==8) print "8   ", $2  }'   8leg* 

awk '{if($1==3) print "2   ",  $2  }'  3leg* 
awk '{if($1==5) print "3   ",  $2  }'  3leg* 
awk '{if($1==7) print "4   ",  $2  }'  4leg* 
awk '{if($1==9) print "5   ", $2 }'   5leg* 
awk '{if($1==11) print "6   ", $2  }'  6leg* 
awk '{if($1==13) print "7   ", $2  }'  7leg* 
awk '{if($1==15) print "8   ", $2  }'  8leg* 

awk '{if($1==4) print "2   ",  $2   }'  3leg* 
awk '{if($1==6) print "3   ",  $2   }'  3leg* 
awk '{if($1==8) print "4   ",  $2   }'  4leg* 
awk '{if($1==10) print "5   ", $2  }'  5leg* 
awk '{if($1==12) print "6   ", $2  }'  6leg* 
awk '{if($1==14) print "7   ", $2  }'  7leg* 
awk '{if($1==16) print "8   ", $2  }'  8leg* 
