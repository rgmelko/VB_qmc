#!/bin/sh

# first run sisipuk:~/VB_qmc/clogl_sharcnet/getnewratiodata.sh then
# run sisipuk:~/VB_qmc/clogl_sharcnet/superduperfiteverything.sh 
# and then run this script.  Ann Kallin - Oct 28, 2011

scp akallin@sisipuk.uwaterloo.ca:~/VB_qmc/clog*/rec[0-6][0-8].dat .;
scp akallin@sisipuk.uwaterloo.ca:~/VB_qmc/clog*/fit[0-6][0-8].dat .;
scp akallin@sisipuk.uwaterloo.ca:~/VB_qmc/clog*/slopes_2.dat ./slopes.dat;
