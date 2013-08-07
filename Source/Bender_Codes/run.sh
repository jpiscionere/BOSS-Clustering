#! /bin/bash

for bin in bin1 bin2 bin3 bin4 
do
time ~/for_deep/measure_Boss_wp ${bin}_Ds.out Di.out mask-cmass-dr11v0-N-Anderson.ply outfile_DsDi_$bin 0.01 10 20 100 1 >${bin}_DsDi_with_norm.out.new_jack &
time ~/for_deep/measure_Boss_wp ${bin}_Ds.out Ri.out mask-cmass-dr11v0-N-Anderson.ply outfile_DsRi_$bin 0.01 10 20 100 2 >${bin}_DsRi_with_norm.out.new_jack &
time ~/for_deep/measure_Boss_wp ${bin}_Rs.out Di.out mask-cmass-dr11v0-N-Anderson.ply outfile_RsDi_$bin 0.01 10 20 100 1 >${bin}_RsDi_with_norm.out.new_jack & 
time ~/for_deep/measure_Boss_wp ${bin}_Rs.out Ri.out mask-cmass-dr11v0-N-Anderson.ply outfile_RsRi_$bin 0.01 10 20 100 2 >${bin}_RsRi_with_norm.out.new_jack &
done


#for bin in bin1 bin2 bin3 bin4
#do
#time ~/for_deep/measure_Boss_wp ${bin}_Ds.out Di.out  /data2/jap/SDSSPix/mask-CMASS-DR9-completeness.pol outfile_DsDi_$bin 0.01 10 20 20 1 >${bin}_DsDi_with_norm.out &
#time ~/for_deep/measure_Boss_wp ${bin}_Ds.out Ri.out  /data2/jap/SDSSPix/mask-CMASS-DR9-completeness.pol outfile_DsRi_$bin 0.01 10 20 20 2 >${bin}_DsRi_with_norm.out &
#time ~/for_deep/measure_Boss_wp ${bin}_Rs.out Di.out  /data2/jap/SDSSPix/mask-CMASS-DR9-completeness.pol outfile_RsDi_$bin 0.01 10 20 20 1 >${bin}_RsDi_with_norm.out & 
#time ~/for_deep/measure_Boss_wp ${bIn}_Rs.out Ri.out  /data2/jap/SDSSPix/mask-CMASS-DR9-completeness.pol outfile_RsRi_$bin 0.01 10 20 20 2 >${bin}_RsRi_with_norm.out &
#done










 
