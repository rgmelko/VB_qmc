#!/bin/sh

scp akallin@whale.sharcnet.ca:/work/akallin/TORUS/3leg/entropy.dat ./entropy03.dat;
scp akallin@whale.sharcnet.ca:/work/akallin/TORUS/4leg/entropy.dat ./entropy04.dat;
scp akallin@whale.sharcnet.ca:/work/akallin/TORUS/5leg/entropy.dat ./entropy05.dat;
scp akallin@whale.sharcnet.ca:/work/akallin/TORUS/6leg/entropy.dat ./entropy06.dat;
scp akallin@whale.sharcnet.ca:/work/akallin/TORUS/7leg/entropy.dat ./entropy07.dat;
scp akallin@whale.sharcnet.ca:/work/akallin/TORUS/8leg/entropy.dat ./entropy08.dat;
scp akallin@whale.sharcnet.ca:/work/akallin/TORUS/9leg/entropy.dat ./entropy09.dat;
scp akallin@whale.sharcnet.ca:/work/akallin/TORUS/10leg/entropy.dat ./entropy10.dat;

average entropy03.dat 11 0 > plot3;
average entropy04.dat 15 0 > plot4;
average entropy05.dat 19 0 > plot5;
average entropy06.dat 23 0 > plot6;
average entropy07.dat 27 0 > plot7;
average entropy08.dat 31 0 > plot8;
average entropy09.dat 35 0 > plot9;
average entropy10.dat 39 0 > plot10;

