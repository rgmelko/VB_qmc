#!/bin/sh

x=5000;
for s in 1 2 3 4 5
    do
    average 08_m0$s\_swap00.dat 8 $x > temp;
    awk '{print $1, -log($2), $3/$2}' temp >  m0$s\.dat;
    echo averaged $s;
    done
   
 
for s in 10 20 30 40
    do
    average 08_m$s\_swap00.dat 8 $x > temp;
    awk '{print $1, -log($2), $3/$2}' temp >  m$s\.dat;
    echo averaged $s;
    done
    

cat m*.dat > allm.dat;

awk '{if($1==1) print NR, $2, $3}' allm.dat > allm1.dat;
awk '{if($1==2) print NR, $2, $3}' allm.dat > allm2.dat;
awk '{if($1==3) print NR, $2, $3}' allm.dat > allm3.dat;
awk '{if($1==4) print NR, $2, $3}' allm.dat > allm4.dat;
awk '{if($1==5) print NR, $2, $3}' allm.dat > allm5.dat;
awk '{if($1==6) print NR, $2, $3}' allm.dat > allm6.dat;
awk '{if($1==7) print NR, $2, $3}' allm.dat > allm7.dat;
awk '{if($1==8) print NR, $2, $3}' allm.dat > allm8.dat;


for s in 1 2 3 4 5 6 7 8
    do
    awk '{if (NR<6) {print NR, $2, $3} else {print ((NR-5)*10), $2, $3}}' allm$s\.dat > temp;
    mv temp allm$s\.dat;
    done

awk '{print $1, sqrt((($2-1.52047512447381)/1.52047512447381)*(($2-1.52047512447381)/1.52047512447381)), $3/1.52047512447381}' allm4.dat > 08plot.dat;
