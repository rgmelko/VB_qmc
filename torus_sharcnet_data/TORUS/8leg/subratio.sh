#!/bin/sh

mkdir 01; cp torus* 01/torus01;  awk '{ if(NR!=6) print $0; if(NR==6) print " 1"}' param.txt > 01/param.txt; cd 01; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 02; cp torus* 02/torus02;  awk '{ if(NR!=6) print $0; if(NR==6) print " 2"}' param.txt > 02/param.txt; cd 02; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 03; cp torus* 03/torus03;  awk '{ if(NR!=6) print $0; if(NR==6) print " 3"}' param.txt > 03/param.txt; cd 03; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 04; cp torus* 04/torus04;  awk '{ if(NR!=6) print $0; if(NR==6) print " 4"}' param.txt > 04/param.txt; cd 04; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 05; cp torus* 05/torus05;  awk '{ if(NR!=6) print $0; if(NR==6) print " 5"}' param.txt > 05/param.txt; cd 05; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 06; cp torus* 06/torus06;  awk '{ if(NR!=6) print $0; if(NR==6) print " 6"}' param.txt > 06/param.txt; cd 06; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 07; cp torus* 07/torus07;  awk '{ if(NR!=6) print $0; if(NR==6) print " 7"}' param.txt > 07/param.txt; cd 07; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 08; cp torus* 08/torus08;  awk '{ if(NR!=6) print $0; if(NR==6) print " 8"}' param.txt > 08/param.txt; cd 08; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 09; cp torus* 09/torus09;  awk '{ if(NR!=6) print $0; if(NR==6) print " 9"}' param.txt > 09/param.txt; cd 09; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 10; cp torus* 10/torus10;  awk '{ if(NR!=6) print $0; if(NR==6) print "10"}' param.txt > 10/param.txt; cd 10; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 11; cp torus* 11/torus11;  awk '{ if(NR!=6) print $0; if(NR==6) print "11"}' param.txt > 11/param.txt; cd 11; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 12; cp torus* 12/torus12;  awk '{ if(NR!=6) print $0; if(NR==6) print "12"}' param.txt > 12/param.txt; cd 12; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 13; cp torus* 13/torus13;  awk '{ if(NR!=6) print $0; if(NR==6) print "13"}' param.txt > 13/param.txt; cd 13; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 14; cp torus* 14/torus14;  awk '{ if(NR!=6) print $0; if(NR==6) print "14"}' param.txt > 14/param.txt; cd 14; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 15; cp torus* 15/torus15;  awk '{ if(NR!=6) print $0; if(NR==6) print "15"}' param.txt > 15/param.txt; cd 15; sqsub -o 0.log -r 7d ./torus*; cd ..; 
