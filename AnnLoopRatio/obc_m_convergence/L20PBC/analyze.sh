#!/bin/sh

x=10;
y=2;
exact=1.20144;

for s in 1 2 3 4 5
    do
    average entrpy20swap0m0$s\.dat 20 $x > temp;
    awk '{if($2!=0) print $1, -log($2), $3/$2}' temp >  m0$s\.dat;
    echo averaged $s;
    done
   
for s in 10 20 30
    do
    average entrpy20swap0m$s\.dat 20 $x > temp;
    awk '{if ($2!=0) print $1, -log($2), $3/$2}' temp >  m$s\.dat;
    echo averaged $s;
    done   

average entrpy20swap0m40.dat 20 $x > temp; 
echo 1 0.693147180559945 0 > m40.dat;
awk '{if ($2!=0) print $1, -log($2)+log(2.0), $3/($2)}' temp >>  m40.dat;
echo averaged 40;

average entrpy20swap0m10b.dat 20 $x > temp; 
echo 1 0.693147180559945 0 > bm10.dat;
awk '{if ($2!=0) print $1, -log($2)+log(2.0), $3/($2)}' temp >>  bm10.dat;
echo averaged 10b;

cat m*.dat > allm.dat;

awk '{if($1==1) print NR, $2, $3}' allm.dat > allm1.dat;
awk '{if($1==2) print NR, $2, $3}' allm.dat > allm2.dat;
awk '{if($1==3) print NR, $2, $3}' allm.dat > allm3.dat;
awk '{if($1==4) print NR, $2, $3}' allm.dat > allm4.dat;
awk '{if($1==5) print NR, $2, $3}' allm.dat > allm5.dat;
awk '{if($1==6) print NR, $2, $3}' allm.dat > allm6.dat;
awk '{if($1==7) print NR, $2, $3}' allm.dat > allm7.dat;
awk '{if($1==8) print NR, $2, $3}' allm.dat > allm8.dat;
awk '{if($1==9) print NR, $2, $3}' allm.dat > allm9.dat;
awk '{if($1==10) print NR, $2, $3}' allm.dat > allm10.dat;


for s in 1 2 3 4 5 6 7 8 9 10
    do
    awk '{if (NR<6) {print NR, $2, $3} else {print ((NR-5)*10), $2, $3}}' allm$s\.dat > temp;
    mv temp allm$s\.dat;
    done

awk -v v1=$exact '{print $1, sqrt((($2-v1)/v1)*(($2-v1)/v1)), $3/v1}' allm$y\.dat > 20plot.dat;
