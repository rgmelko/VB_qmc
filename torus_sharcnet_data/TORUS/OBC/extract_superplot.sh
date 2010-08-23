#!/bin/sh

awk '{if(NR==6) print "3   ",  $2/3 , $3/3 }'  3leg/ratiod 
awk '{if(NR==8) print "4   ",  $2/4 , $3/4 }'  4leg/ratiod 
awk '{if(NR==10) print "5   ", $2/5 , $3/5 }'  5leg/ratiod 
awk '{if(NR==12) print "6   ", $2/6 , $3/6 }'  6leg/ratiod 
awk '{if(NR==14) print "7   ", $2/7 , $3/7 }'  7leg/ratiod 
awk '{if(NR==16) print "8   ", $2/8 , $3/8 }'  8leg/ratiod 
awk '{if(NR==18) print "9   ", $2/9 , $3/9 }'  9leg/ratiod 
awk '{if(NR==20) print "10  ", $2/10, $3/10}'  10leg/ratiod 
awk '{if(NR==22) print "11  ", $2/11, $3/11}'  11leg/ratiod 
