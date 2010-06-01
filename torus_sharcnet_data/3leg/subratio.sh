#!/bin/sh

mkdir 01; cp torus* 01/torus01;  awk '{ if(NR!=6) print $0; if(NR==6) print " 1"}' param.txt > 01/param.txt; cd 01; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 02; cp torus* 02/torus02;  awk '{ if(NR!=6) print $0; if(NR==6) print " 2"}' param.txt > 02/param.txt; cd 02; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 03; cp torus* 03/torus03;  awk '{ if(NR!=6) print $0; if(NR==6) print " 3"}' param.txt > 03/param.txt; cd 03; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 04; cp torus* 04/torus04;  awk '{ if(NR!=6) print $0; if(NR==6) print " 4"}' param.txt > 04/param.txt; cd 04; sqsub -o 0.log -r 7d ./torus*; cd ..; 
mkdir 05; cp torus* 05/torus05;  awk '{ if(NR!=6) print $0; if(NR==6) print " 5"}' param.txt > 05/param.txt; cd 05; sqsub -o 0.log -r 7d ./torus*; cd ..; 
