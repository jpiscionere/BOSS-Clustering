#! /bin/bash
polygon=/hd0/Research/BOSS-dr9/BOSSPix/mask-CMASS-DR9-completeness.pol

for bin in bin1 bin2 bin3 bin4 
do
time ./measure_Boss_wp Data/spectro_sample_${bin}_new.txt Data/imaging_sample_NS_dr9.new.txt $polygon outfile_DsDi_$bin 0.01 10 20 100 1 >${bin}_DsDi_with_norm.new.out &
time ./measure_Boss_wp Data/spectro_sample_${bin}_new.txt Data/imaging_randoms.txt $polygon outfile_DsRi_$bin 0.01 10 20 100 2 >${bin}_DsRi_with_norm.new.out &
time ./measure_Boss_wp Data/randoms_spectro_$bin.inv.out Data/imaging_sample_NS_dr9.new.txt $polygon outfile_RsDi_$bin 0.01 10 20 100 1 >${bin}_RsDi_with_norm.new.out & 
time ./measure_Boss_wp Data/randoms_spectro_$bin.inv.out Data/imaging_randoms.txt $polygon outfile_RsRi_$bin 0.01 10 20 100 2 >${bin}_RsRi_with_norm.new.out &
done


#for bin in bin1 bin2 bin3 bin4
#do
#DD
#time ./measure_Boss_wp_ls Data/spectro_sample_$bin.txt Data/spectro_sample_$bin.txt /hd0/Research/BOSS-dr9/BOSSPix/mask-CMASS-DR9-completeness.pol outfile_DsDs_$bin 0.01 1 20 20 1 >${bin}_DsDs_ls.out &
#DR
#time ./measure_Boss_wp_ls Data/spectro_sample_$bin.txt Data/randoms_spectro_$bin.out /hd0/Research/BOSS-dr9/BOSSPix/mask-CMASS-DR9-completeness.pol outfile_DsRs_$bin 0.01 1 20 20 2 >${bin}_DsRs_ls.out &
#RR
#time  ./measure_Boss_wp_ls Data/randoms_spectro_$bin.out Data/randoms_spectro_$bin.out /hd0/Research/BOSS-dr9/BOSSPix/mask-CMASS-DR9-completeness.pol outfile_RsRs_$bin 0.01 1 20 20 2 >${bin}_RsRs_ls.out &


#done


 
