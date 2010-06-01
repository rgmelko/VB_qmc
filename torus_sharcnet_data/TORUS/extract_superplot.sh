#!/bin/sh

awk '{if(NR==6) print "3   ", $2, $3}'  3leg/ratiod 
awk '{if(NR==8) print "4   ", $2, $3}'  4leg/ratiod 
awk '{if(NR==10) print "5   ", $2, $3}' 5leg/ratiod 
awk '{if(NR==12) print "6   ", $2, $3}' 6leg/ratiod 
awk '{if(NR==14) print "7   ", $2, $3}' 7leg/ratiod 
awk '{if(NR==16) print "8   ", $2, $3}' 8leg/ratiod 
awk '{if(NR==18) print "9   ", $2, $3}' 9leg/ratiod 
awk '{if(NR==20) print "10  ", $2, $3}' 10leg/ratiod 
