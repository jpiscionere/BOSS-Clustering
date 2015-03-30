#!/bin/bash

path=/hd0/Research/Clustering/Boss/dr11/dr11v2
imaging_randoms=$path/Ri.dr11v2.out
mask=$path/mask-cmass-dr11v2-N-Anderson.ply
area_tot=2.050670
imaging_data=$path/Di.dr11v2_selection.out

nbins=20
r_min=0.01


for bin in bin1 bin2 bin3 bin4 
do


data=$path/${bin}_Ds_dr11v2_selection.out
randoms=$path/${bin}_Rs_dr11v2.out




time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp $data $imaging_data $r_min 10 $nbins 1 $area_tot >$path/${bin}_DsDi_nojack.out
time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp $data $imaging_randoms  $r_min 10 $nbins 2 $area_tot >$path/${bin}_DsRi_nojack.out
time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp $randoms $imaging_data  $r_min 10 $nbins 1 $area_tot >$path/${bin}_RsDi_nojack.out
time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp $randoms $imaging_randoms $r_min 10 $nbins 2 $area_tot >$path/${bin}_RsRi_nojack.out

done
