#!/bin/sh


       sqsub -o 0.log -r 7d ./obc_torus*;        
cd 01; sqsub -o 0.log -r 7d ./obc_torus*; cd ..; 
cd 02; sqsub -o 0.log -r 7d ./obc_torus*; cd ..; 
cd 03; sqsub -o 0.log -r 7d ./obc_torus*; cd ..; 
cd 04; sqsub -o 0.log -r 7d ./obc_torus*; cd ..; 
cd 05; sqsub -o 0.log -r 7d ./obc_torus*; cd ..; 
