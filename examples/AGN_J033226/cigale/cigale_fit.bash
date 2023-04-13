# !/bin/bash

cd filters/

pcigale-filters add ACS_F775W.txt
pcigale-filters add ACS_F850LP.txt
pcigale-filters add ISAAC_Ks.txt
pcigale-filters add xraybandpass1.dat
pcigale-filters add xraybandpass2.dat
pcigale-filters add xraybandpass3.dat
pcigale-filters add xraybandpass4.dat
pcigale-filters add xraybandpass5.dat
pcigale-filters add xraybandpass6.dat
pcigale-filters add xraybandpass7.dat
pcigale-filters add xraybandpass8.dat
pcigale-filters add xraybandpass9.dat
pcigale-filters add xraybandpass10.dat
pcigale-filters add xraybandpass11.dat
pcigale-filters add xraybandpass12.dat
pcigale-filters add xraybandpass13.dat
pcigale-filters add xraybandpass14.dat
pcigale-filters add xraybandpass15.dat

cd ..

pcigale run

pcigale-plots sed