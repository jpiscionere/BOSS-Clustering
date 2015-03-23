#! /bin/bash

data=/hd0/Research/Clustering/Boss/DR3/dr3_Ds_zrange_trim2
randoms=/hd0/Research/Clustering/Boss/DR3/dr3_Rs_zrange_trim
imaging_data=/hd0/Research/Clustering/Boss/DR3/dr3_Di_primtarget_no
imaging_randoms=/hd0/Research/Clustering/Boss/DR3/dr3_Ri_short2
mask=window.sample14.trimmed.ply 

nbins=17
area_tot=0.9076429
r_min=0.005

#time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_openmp $data $imaging_data $mask outfile_DsDi $r_min 10 $nbins 50 1 >dr3_DsDi.out 
#time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_openmp $data $imaging_randoms $mask outfile_DsRi  $r_min 10 $nbins 50 2 >dr3_DsRi.out 
#time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_openmp $randoms $imaging_data $mask outfile_RsDi  $r_min 10 $nbins 50 1 >dr3_RsDi.out 
#time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_openmp $randoms $imaging_randoms $mask outfile_RsRi  $r_min 10 $nbins 50 2 >dr3_RsRi.out 

#time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp $data $imaging_data $r_min 10 $nbins 1 $area_tot >dr3_DsDi_nojack_15.out
#time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp_mag_cut $data $imaging_randoms  $r_min 10 $nbins 2 $area_tot >dr3_DsRi_nojack_15.out 
#time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp $randoms $imaging_data  $r_min 10 $nbins 1 $area_tot >dr3_RsDi_nojack_15.out
#time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp_mag_cut $randoms $imaging_randoms $r_min 10 $nbins 2 $area_tot >dr3_RsRi_nojack_15.out

time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp_mag_cut $data dr3_Di_mag_extinction $r_min 10 $nbins 1 $area_tot >dr3_DsDi_mag_cut
time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp_mag_cut $randoms dr3_Di_mag_extinction $r_min 10 $nbins 1 $area_tot >dr3_RsDi_mag_cut
time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp_mag_cut $data $imaging_randoms  $r_min 10 $nbins 2 $area_tot >dr3_DsRi_mag_cut
time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp_mag_cut $randoms $imaging_randoms $r_min 10 $nbins 2 $area_tot >dr3_RsRi_mag_cut




