#!/bin/sh

rsync -turav --exclude '*.log' akallin@whale.sharcnet.ca:/work/akallin/TORUS/ .;

average    3leg/entropy.dat 11 0 > 3leg/avgent;
average 3leg/01/entropy.dat 11 0 > 3leg/01/avgent;
average 3leg/02/entropy.dat 11 0 > 3leg/02/avgent;
average 3leg/03/entropy.dat 11 0 > 3leg/03/avgent;
average 3leg/04/entropy.dat 11 0 > 3leg/04/avgent;
average 3leg/05/entropy.dat 11 0 > 3leg/05/avgent;
cd 3leg; ./extract.sh > ratiod; cd ..;

average    4leg/entropy.dat 15 0 > 4leg/avgent;
average 4leg/01/entropy.dat 15 0 > 4leg/01/avgent;
average 4leg/02/entropy.dat 15 0 > 4leg/02/avgent;
average 4leg/03/entropy.dat 15 0 > 4leg/03/avgent;
average 4leg/04/entropy.dat 15 0 > 4leg/04/avgent;
average 4leg/05/entropy.dat 15 0 > 4leg/05/avgent;
average 4leg/06/entropy.dat 15 0 > 4leg/06/avgent;
average 4leg/07/entropy.dat 15 0 > 4leg/07/avgent;
cd 4leg; ./extract.sh > ratiod; cd ..;

average    5leg/entropy.dat 19 0 > 5leg/avgent;
average 5leg/01/entropy.dat 19 0 > 5leg/01/avgent;
average 5leg/02/entropy.dat 19 0 > 5leg/02/avgent;
average 5leg/03/entropy.dat 19 0 > 5leg/03/avgent;
average 5leg/04/entropy.dat 19 0 > 5leg/04/avgent;
average 5leg/05/entropy.dat 19 0 > 5leg/05/avgent;
average 5leg/06/entropy.dat 19 0 > 5leg/06/avgent;
average 5leg/07/entropy.dat 19 0 > 5leg/07/avgent;
average 5leg/08/entropy.dat 19 0 > 5leg/08/avgent;
average 5leg/09/entropy.dat 19 0 > 5leg/09/avgent;
cd 5leg; ./extract.sh > ratiod; cd ..;



average    6leg/entropy.dat 23 0 > 6leg/avgent;
average 6leg/01/entropy.dat 23 0 > 6leg/01/avgent;
average 6leg/02/entropy.dat 23 0 > 6leg/02/avgent;
average 6leg/03/entropy.dat 23 0 > 6leg/03/avgent;
average 6leg/04/entropy.dat 23 0 > 6leg/04/avgent;
average 6leg/05/entropy.dat 23 0 > 6leg/05/avgent;
average 6leg/06/entropy.dat 23 0 > 6leg/06/avgent;
average 6leg/07/entropy.dat 23 0 > 6leg/07/avgent;
average 6leg/08/entropy.dat 23 0 > 6leg/08/avgent;
average 6leg/09/entropy.dat 23 0 > 6leg/09/avgent;
average 6leg/10/entropy.dat 23 0 > 6leg/10/avgent;
average 6leg/11/entropy.dat 23 0 > 6leg/11/avgent;
cd 6leg; ./extract.sh > ratiod; cd ..;


average    7leg/entropy.dat 27 0 > 7leg/avgent;
average 7leg/01/entropy.dat 27 0 > 7leg/01/avgent;
average 7leg/02/entropy.dat 27 0 > 7leg/02/avgent;
average 7leg/03/entropy.dat 27 0 > 7leg/03/avgent;
average 7leg/04/entropy.dat 27 0 > 7leg/04/avgent;
average 7leg/05/entropy.dat 27 0 > 7leg/05/avgent;
average 7leg/06/entropy.dat 27 0 > 7leg/06/avgent;
average 7leg/07/entropy.dat 27 0 > 7leg/07/avgent;
average 7leg/08/entropy.dat 27 0 > 7leg/08/avgent;
average 7leg/09/entropy.dat 27 0 > 7leg/09/avgent;
average 7leg/10/entropy.dat 27 0 > 7leg/10/avgent;
average 7leg/11/entropy.dat 27 0 > 7leg/11/avgent;
average 7leg/12/entropy.dat 27 0 > 7leg/12/avgent;
average 7leg/13/entropy.dat 27 0 > 7leg/13/avgent;
cd 7leg; ./extract.sh > ratiod; cd ..;


average    8leg/entropy.dat 31 0 > 8leg/avgent;
average 8leg/01/entropy.dat 31 0 > 8leg/01/avgent;
average 8leg/02/entropy.dat 31 0 > 8leg/02/avgent;
average 8leg/03/entropy.dat 31 0 > 8leg/03/avgent;
average 8leg/04/entropy.dat 31 0 > 8leg/04/avgent;
average 8leg/05/entropy.dat 31 0 > 8leg/05/avgent;
average 8leg/06/entropy.dat 31 0 > 8leg/06/avgent;
average 8leg/07/entropy.dat 31 0 > 8leg/07/avgent;
average 8leg/08/entropy.dat 31 0 > 8leg/08/avgent;
average 8leg/09/entropy.dat 31 0 > 8leg/09/avgent;
average 8leg/10/entropy.dat 31 0 > 8leg/10/avgent;
average 8leg/11/entropy.dat 31 0 > 8leg/11/avgent;
average 8leg/12/entropy.dat 31 0 > 8leg/12/avgent;
average 8leg/13/entropy.dat 31 0 > 8leg/13/avgent;
average 8leg/14/entropy.dat 31 0 > 8leg/14/avgent;
average 8leg/15/entropy.dat 31 0 > 8leg/15/avgent;
cd 8leg; ./extract.sh > ratiod; cd ..;


average    9leg/entropy.dat 35 0 > 9leg/avgent;
average 9leg/01/entropy.dat 35 0 > 9leg/01/avgent;
average 9leg/02/entropy.dat 35 0 > 9leg/02/avgent;
average 9leg/03/entropy.dat 35 0 > 9leg/03/avgent;
average 9leg/04/entropy.dat 35 0 > 9leg/04/avgent;
average 9leg/05/entropy.dat 35 0 > 9leg/05/avgent;
average 9leg/06/entropy.dat 35 0 > 9leg/06/avgent;
average 9leg/07/entropy.dat 35 0 > 9leg/07/avgent;
average 9leg/08/entropy.dat 35 0 > 9leg/08/avgent;
average 9leg/09/entropy.dat 35 0 > 9leg/09/avgent;
average 9leg/10/entropy.dat 35 0 > 9leg/10/avgent;
average 9leg/11/entropy.dat 35 0 > 9leg/11/avgent;
average 9leg/12/entropy.dat 35 0 > 9leg/12/avgent;
average 9leg/13/entropy.dat 35 0 > 9leg/13/avgent;
average 9leg/14/entropy.dat 35 0 > 9leg/14/avgent;
average 9leg/15/entropy.dat 35 0 > 9leg/15/avgent;
average 9leg/16/entropy.dat 35 0 > 9leg/16/avgent;
average 9leg/17/entropy.dat 35 0 > 9leg/17/avgent;
cd 9leg; ./extract.sh > ratiod; cd ..;


average    10leg/entropy.dat 39 0 > 10leg/avgent;
average 10leg/01/entropy.dat 39 0 > 10leg/01/avgent;
average 10leg/02/entropy.dat 39 0 > 10leg/02/avgent;
average 10leg/03/entropy.dat 39 0 > 10leg/03/avgent;
average 10leg/04/entropy.dat 39 0 > 10leg/04/avgent;
average 10leg/05/entropy.dat 39 0 > 10leg/05/avgent;
average 10leg/06/entropy.dat 39 0 > 10leg/06/avgent;
average 10leg/07/entropy.dat 39 0 > 10leg/07/avgent;
average 10leg/08/entropy.dat 39 0 > 10leg/08/avgent;
average 10leg/09/entropy.dat 39 0 > 10leg/09/avgent;
average 10leg/10/entropy.dat 39 0 > 10leg/10/avgent;
average 10leg/11/entropy.dat 39 0 > 10leg/11/avgent;
average 10leg/12/entropy.dat 39 0 > 10leg/12/avgent;
average 10leg/13/entropy.dat 39 0 > 10leg/13/avgent;
average 10leg/14/entropy.dat 39 0 > 10leg/14/avgent;
average 10leg/15/entropy.dat 39 0 > 10leg/15/avgent;
average 10leg/16/entropy.dat 39 0 > 10leg/16/avgent;
average 10leg/17/entropy.dat 39 0 > 10leg/17/avgent;
average 10leg/18/entropy.dat 39 0 > 10leg/18/avgent;
average 10leg/19/entropy.dat 39 0 > 10leg/19/avgent;
cd 10leg; ./extract.sh > ratiod; cd ..;

./extract_superplot.sh > superplot;


cd OBC;
./justmakedata.sh;
cd ..
