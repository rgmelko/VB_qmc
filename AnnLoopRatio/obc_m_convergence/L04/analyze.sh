#!/bin/sh

x=5000;
for s in 1 2 3 4 5
    do
    average entrpy04_m$s\swap00.dat 4 $x > temp;
    awk '{print $1, -log($2), $3/$2}' temp >  m0$s\.dat;
    echo averaged $s;
    done
   
 
for s in 10 20 30 40
    do
    average entrpy04_m$s\swap00.dat 4 $x > temp;
    awk '{print $1, -log($2), $3/$2}' temp >  m$s\.dat;
    echo averaged $s;
    done   

cat m*.dat > allm.dat;

awk '{if($1==1) print NR, $2, $3}' allm.dat > allm1.dat;
awk '{if($1==2) print NR, $2, $3}' allm.dat > allm2.dat;
awk '{if($1==3) print NR, $2, $3}' allm.dat > allm3.dat;
awk '{if($1==4) print NR, $2, $3}' allm.dat > allm4.dat;


for s in 1 2 3 4
    do
    awk '{if (NR<6) {print NR, $2, $3} else {print ((NR-5)*10), $2, $3}}' allm$s\.dat > temp;
    mv temp allm$s\.dat;
    done

awk '{print $1, sqrt((($2-0.607568331649827)/0.607568331649827)*(($2-0.607568331649827)/0.607568331649827)), $3/0.607568331649827}' allm2.dat > 04plot.dat;
