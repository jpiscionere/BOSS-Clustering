#! /bin/bash





nbin=20
n_jack=150

path=/hd0/Research/Clustering/Boss/dr11/dr11v2


nden[1]=0.0002494737166611
nden[2]=0.0003150319712696
nden[3]=0.0002238406773599
nden[4]=0.0000912964784524
nden[5]=0.0001847687542961
nden[6]=0.0002796705422133
nden[7]=0.0001307780091214

i=1

for bin in bin1 bin2 bin3 bin4 bin_all bin12 bin34
do

nden1=${nden[$i]}
echo $nden1

~/Clustering/Boss/Source/Vpac_Codes/calculate_wp $path/${bin}_DsDi_with_norm.out.new $path/${bin}_DsRi_with_norm.out.new $path/${bin}_RsDi_with_norm.out.new $path/${bin}_RsRi_with_norm.out.new $nbin $n_jack $nden1 1 $path/${bin}_covar_inv_jack150 >$path/${bin}_wp_dr11v2_new


i=`expr $i + 1`

done













