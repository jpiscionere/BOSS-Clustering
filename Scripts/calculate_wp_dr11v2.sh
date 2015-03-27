#! /bin/bash





nbin=20
n_jack=150

path=/hd0/Research/Clustering/Boss/dr11/dr11v2


nden[1]=0.0002494737166611
nden[2]=0.0003150319712696
nden[3]=0.0002238406773599
nden[4]=0.0000912964784524
nden[5]=0.0001847687542961


i=1

for bin in bin1 bin2 bin3 bin4 bin_all
do

nden1=${nden[$i]}
echo $nden1

~/Clustering/Boss/Source/Vpac_Codes/calculate_wp $path/${bin}_DsDi_with_norm.out.new $path/${bin}_DsRi_with_norm.out.new $path/${bin}_RsDi_with_norm.out.new $path/${bin}_RsRi_with_norm.out.new $nbin $n_jack $nden1 1 $path/${bin}_covar_inv_jack150 >$path/${bin}_wp_dr11v2_new


i=`expr $i + 1`

done


#~/Clustering/Boss/Source/Vpac_Codes/calculate_wp bin1_DsDi_with_norm.out.new_jack150 bin1_DsRi_with_norm.out.new_jack150 bin1_RsDi_with_norm.out.new_jack150 bin1_RsRi_with_norm.out.new_jack150 20 $n_jack $nden1 1  bin1_covar_inv_with_norm.dr11.out.new_jack150 >bin1_wp_dr11v2_new
#~/Clustering/Boss/Source/Vpac_Codes/calculate_wp bin2_DsDi_with_norm.out.new_jack150 bin2_DsRi_with_norm.out.new_jack150 bin2_RsDi_with_norm.out.new_jack150 bin2_RsRi_with_norm.out.new_jack150 20 $n_jack $nden2  1  bin2_covar_inv_with_norm.dr11.out.new_jack150 >bin2_wp_dr11v2_new
#~/Clustering/Boss/Source/Vpac_Codes/calculate_wp bin3_DsDi_with_norm.out.new_jack150 bin3_DsRi_with_norm.out.new_jack150 bin3_RsDi_with_norm.out.new_jack150 bin3_RsRi_with_norm.out.new_jack150 20 $n_jack $nden3 1 bin3_covar_inv_with_norm.dr11.out.new_jack150 >bin3_wp_dr11v2_new
#~/Clustering/Boss/Source/Vpac_Codes/calculate_wp bin4_DsDi_with_norm.out.new_jack150 bin4_DsRi_with_norm.out.new_jack150 bin4_RsDi_with_norm.out.new_jack150 bin4_RsRi_with_norm.out.new_jack150 20 $n_jack $nden4 1 bin4_covar_inv_with_norm.dr11.out.new_jack150 >bin4_wp_dr11v2_new
#~/Clustering/Boss/Source/Vpac_Codes/calculate_wp bin_all_DsDi_with_norm.out.new_jack150 bin_all_DsRi_with_norm.out.new_jack150 bin_all_RsDi_with_norm.out.new_jack150 bin_all_RsRi_with_norm.out.new_jack150 20 $n_jack 0.000185 1 bin_all_covar_inv_with_norm.dr11.out.new_jack150 >bin_all_wp_dr11v2_full_map_jack150_littleh.out






