#! /bin/bash

nden=0.0001260
echo $nden

nden==0.0000892803189156


paste dr3_DsDi_mag_cut dr3_DsRi_mag_cut dr3_RsDi_mag_cut dr3_RsRi_mag_cut >tmp 
awk '{wp=$2/$5-$8/$11; {print $1,wp,wp/'$nden'}}'<tmp >wp_dr3_mag_cut





