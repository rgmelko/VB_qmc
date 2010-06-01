#!/bin/sh

awk '{if(NR==6) print "3   ", $2, $3}' plot3;
awk '{if(NR==8) print "4   ", $2, $3}' plot4;
awk '{if(NR==10) print "5   ", $2, $3}' plot5;
awk '{if(NR==12) print "6   ", $2, $3}' plot6;
awk '{if(NR==14) print "7   ", $2, $3}' plot7;
awk '{if(NR==16) print "8   ", $2, $3}' plot8;
awk '{if(NR==18) print "9   ", $2, $3}' plot9;
awk '{if(NR==20) print "10  ", $2, $3}' plot10;
