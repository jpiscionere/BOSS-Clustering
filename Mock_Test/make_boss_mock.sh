#!/bin/bash

logMmin=13.31
logM0=0.57
siglogM=0.58
logM1=14.30
alpha=1.56
Mstar=2E12
fgal=1
gamma=1
delta_vir=264.
redshift=0.5
halo_file=/ssd1/Research/halo_files/Carmen/High_Z
map_dir=/home/piscioja/SDSSPix/Maps
map_full=$map_dir/mask-cmass-dr11v2-N-Anderson.pix
polygon_file=/hd0/Research/Clustering/Boss/dr11/dr11v2/mask-cmass-dr11v2-N-Anderson.ply
zmin=0.02
zmax=0.225
zspace=0
trim=$1
i=$2
fibcol=$3

#for i in $(seq 2021 1 2031)
#do


halobias_fof_nfw_quiet	 3 4 1 $logMmin $siglogM $logM0 $logM1 $alpha 1 $gamma $fgal -21 $delta_vir $Mstar $redshift  trashfile.out $i < $halo_file/Carmen_${i}_z0p520_fof_b0p2.dpp.halos >output_$trim.$i 

	if [[ $trim == 2 ]]
	then
    		echo "....................Using Polygon File To Trim"
	
		tag="sdss_mangle" 
		
		 ~/halobias/makemock_double 0 $zspace 0 0 0 0 $zmin $zmax $map_full 0.6 0.52 < output_$trim.$i > boss_mock_$i.$tag.galaxies		

		/hd0/Research/Mangle/mangle2.2/bin/polyid $polygon_file boss_mock_$i.$tag.galaxies outfile.$i

		awk '{if(NR>1) {if(NF==3){print $0} else {print $0,0}}}'< outfile.$i >no_header.$i 

		paste boss_mock_$i.$tag.galaxies no_header.$i > pasted.$i 

		awk '{if($6 >0 ){print $1,$2,$3}}'<pasted.$i > boss_mock_$i.$tag.galaxies
		awk '{if($6 >0 ){print $1,$2,$3,1,$6}}'<pasted.$i > boss_mock_$i.$tag.sector.galaxies
		awk '{if($6 >0 ){print $1,$2,$6}}'<pasted.$i > boss_mock_$i.$tag.imaging.galaxies
			
	
		ND1=`wc -l boss_mock_$i.$tag.galaxies | cut -f1 -d' '`
		echo "$ND1 Number of Mangle Filtered Galaxies Used"

		rm outfile_$trim.$i no_header.$i pasted.$i outfile.$i
	
	elif [[ $trim == 1 ]];
	then
    		echo "......................Trimmed using Pixmap"
		tag="sdss_pix_map"
		 ~/halobias/makemock_double $trim $zspace 0 0 0 0 $zmin $zmax $map_full 0.6 0.52 < output_$trim.$i > boss_mock_$i.$tag.galaxies	
	else
		echo ".......................Not Trimming, Output Sphere"
		tag="sphere"
		 ~/halobias/makemock_double $trim $zspace 0 0 0 0 $zmin $zmax $map_dir/lss_geometry.fullsphere_single_index.pix 0.6 0.52 < output_$trim.$i > boss_mock_$i.$tag.galaxies	
	fi

	if [[ $fibcol == 1 ]]
	then
		tagfibcol=$tag.fibcol
		echo "......................Adding Fiber Collisions"
		IDfib 0 1 162.0 boss_mock_$i.$tag.galaxies > tmp.$i
#		awk '{if($4==-1){print $1,$2,$3} else {print $1,$2,$4}}' <tmp.$i >boss_mock_$i.$tag.galaxies
		awk '{if($4==-1){print $1,$2,$3}}' <tmp.$i >boss_mock_$i.$tagfibcol.galaxies
	fi	

DDrppi 0.01 10 20 boss_mock_$i.$tag.galaxies a boss_mock_$i.$tag.galaxies a 500 >DD_$tag.$i  

rm -r output_$trim.$i 

#done
#./clf 1 13.31 1 1 -3004  clf_pnm hb_pnmfile  /net/bender/data0/LasDamas/Carmen/2021/Carmen_2021_z0p520_fof_b0p2.*.bgc < /net/bender/data0/LasDamas/Carmen/2021/Carmen_2021_z0p520_fof_b0p2.dpp.halos > fff
#makerandommock/make_random_mock 1000 25000000 0 2 0 0 0 0 0.45 0.6 ~/SDSSPix/Maps/mask-cmass-dr11v2-N-Anderson.pix 0.6 -300 0 0.45 >random_box

#./data_wp/DDrppi 0.01 10 10 /hd0/Research/Clustering/Boss/dr11/dr11v2/rands_for_mock_test /hd0/Research/Clustering/Boss/dr11/dr11v2/rands_for_mock_test 40 >RR_$tag 

Rand_file=/hd0/Research/Clustering/Boss/Mock_Test/Randoms_$tag
NR1=`wc -l $Rand_file | cut -f1 -d' '`
echo "$NR1 Number of Randoms Galaxies Used"

#DDrppi 0.01 10 20 $Rand_file $Rand_file 40 >RR_$tag 

echo "$NR1 Number of Randoms Used"
#for i in $(seq 2021 1 2031)

#do

DDrppi 0.01 10 20 boss_mock_$i.$tag.galaxies a $Rand_file a 500  >D_$tag.${i}.R_$tag 


#done



#for i in $(seq 2021 1 2031)
#do

ND1=`wc -l boss_mock_$i.$tag.galaxies | cut -f1 -d' '`

echo "wprp 20 $ND1 $ND1 $NR1 $NR1 DD_$tag.$i D_$tag.${i}.R_$tag D_$tag.${i}.R_$tag RR_$tag 40 >wprp_$i.$tag.$fibcol.out"

wprp 20 $ND1 $ND1 $NR1 $NR1 DD_$tag.$i D_$tag.${i}.R_$tag D_$tag.${i}.R_$tag RR_$tag 500 >test_$i.$tag.$fibcol.out.test 

#rm boss_mock_$i.$tag.galaxies
#done


