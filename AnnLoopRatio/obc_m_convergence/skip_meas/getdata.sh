#!/bin/sh

rm en*swap*dat;

rsync -turav akallin@brown.sharcnet.ca:/work/akallin/loopratio/rectanglem10/obc_dmrg_test/code/skip_measurements/2test/00/en* .;
mv ene* 2_energy.dat;
mv ent* 2_entrpy.dat;

rsync -turav akallin@brown.sharcnet.ca:/work/akallin/loopratio/rectanglem10/obc_dmrg_test/code/skip_measurements/3test/en* .;
mv ene* 3_energy.dat;
mv ent* 3_entrpy.dat;

rsync -turav akallin@brown.sharcnet.ca:/work/akallin/loopratio/rectanglem10/obc_dmrg_test/code/skip_measurements/4test/en* .;
mv ene* 4_energy.dat;
mv ent* 4_entrpy.dat;
