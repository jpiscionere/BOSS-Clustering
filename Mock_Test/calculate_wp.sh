#! /bin/bash

for i in $(seq 2021 1 2021)
do
ND1=`wc -l boss_data_test_$i.out | cut -f1 -d' '`
nden=`echo "$ND1/(3.0/(4*3.14159265359)*11210000000/2)" |bc -l `
volume=`echo "4.0*3.14159265359/3.0*(417.1*417.1*417.1 )"|bc -l`
nden2=`echo "$ND1/$volume"|bc -l`
im_norm=`echo "2500000/(4*3.14159265359)*(0.7281599)" | bc -l`
echo "$nden $nden2"
echo "$im_norm"

paste ${i}_DsDi_nojack.out.short ${i}_DsRi_nojack.out.retry ${i}_RsDi_nojack.out.short RsRi_nojack.out.retry >tmp
#paste ${i}_DsDi_nojack.out.short ${i}_DsRi_nojack.out.short.norm_2 ${i}_RsDi_nojack.out.short RsRi_nojack.out.norm_2 >tmp
#awk '{wp=$2/($5/'$im_norm')-$8/($11/'$im_norm'); {print $1,wp,wp/'$nden'}}'<tmp >wp_nojack_short_$i.norm_1.retry.out
awk '{wp=$2/$5-$8/$11; {print $1,wp,wp/'$nden'}}'<tmp >wp_nojack_short_$i.norm_2.retry.out



done

