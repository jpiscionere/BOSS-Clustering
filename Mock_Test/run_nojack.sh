#! /bin/bash
#mask=/hd0/Research/Clustering/Boss/dr11/dr11v2/allsky.pol
mask=/hd0/Research/Clustering/Boss/dr11/dr11v2/mask-cmass-dr11v2-N-Anderson.ply
bin=data_sphere_test
#time measure_Boss_wp random_mock_test random_mock_test_imaging $mask outfile_RsRi_$bin 0.01 10 20 100 2 >${bin}_RsRi_with_norm.out &
Rs=/hd0/Research/Clustering/Boss/Mock_Test/Rs_0
Ri=/hd0/Research/Clustering/Boss/Mock_Test/Ri_0
/home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojackopenmp $Rs $Ri 0.01 10 20 2 > RsRi_nojack.out.littleh & 



for i in $(seq 2021 1 2031)
do

galaxy_file=boss_mock_$i.sphere.galaxies

#/hd0/Research/Mangle/mangle2.2/bin/polyid $mask  $galaxy_file outfile
#awk '{if(NR>1) {if(NF==3){print $0} else {print $0,0}}}'< outfile >no_header
#paste $galaxy_file no_header > pasted
#awk '{if($6 >0 ){print $1,$2,$3,1,$6}}'<pasted > boss_data_test_$i.out
#awk '{if($6 >0 ){print $1,$2,$6}}'<pasted > boss_imaging_test_$i.out

awk '{print $0,1.0,1.0}' <$galaxy_file >boss_data_test_$i.out
awk '{print $1,$2,1.0}' <$galaxy_file > boss_imaging_test_$i.out

Ds=boss_data_test_$i.out
Di=boss_imaging_test_$i.out

#R=`echo "$i - 2021" | bc -l `

#Rs=/hd0/Research/Clustering/Boss/Mock_Test/Rs_$R
#Ri=/hd0/Research/Clustering/Boss/Mock_Test/Ri_$R



#done
/home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp $Ds $Di 0.01 10 20 1 > ${i}_DsDi_nojack.out.littleh &
/home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp $Ds $Ri 0.01 10 20 2 > ${i}_DsRi_nojack.out.littleh &
/home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack_openmp $Rs $Di 0.01 10 20 1 > ${i}_RsDi_nojack.out.littleh &
#/home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack $Rs $Ri 0.01 10 20 2 > RsRi_nojack.out.new_norm & 
#/home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_nojack $Rs $Ri 0.01 10 20 1 > RsRi_nojack.out & 

#time /home/piscioja/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp $Ds $Di $mask outfile_DsDi_$i 0.01 10 20 100 1 >${i}_DsDi_with_norm.out &
#time measure_Boss_wp $Ds $Ri $mask outfile_DsRi_$i 0.01 10 20 100 2 >${i}_DsRi_with_norm.out &
#time measure_Boss_wp $Rs $Di $mask outfile_RsDi_$i 0.01 10 20 100 1 >${i}_RsDi_with_norm.out & 
#time measure_Boss_wp boss_random_data_test.out /hd0/Research/Clustering/Boss/dr11/dr11v2/Ri.dr11v2.out $mask outfile_RsRi_$bin 0.01 10 20 100 2 >${bin}_RsRi_with_norm.out &

done













