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
halo_file=/net/bender/data0/LasDamas/Carmen
map_dir=/home/piscioja/SDSSPix/Maps
map_full=$map_dir/mask-cmass-dr11v2-N-Anderson.pix
polygon_file=/hd0/Research/Clustering/Boss/dr11/dr11v2/mask-cmass-dr11v2-N-Anderson.ply
zmin=0.02
zmax=0.225
zspace=2
trim=$1
#i=$2
fibcol=$3
i=2021

echo "../../FibCol/clf 1 13.31 1 1 -$i  clf_pnm$i hb_pnmfile$i $halo_file/$i/Carmen_${i}_z0p520_fof_b0p2.*.bgc <$halo_file/$i/Carmen_${i}_z0p520_fof_b0p2.dpp.halos  > fff_$i"	
../../FibCol/clf 1 13.31 1 1 -$i  clf_pnm$i hb_pnmfile$i $halo_file/$i/Carmen_${i}_z0p520_fof_b0p2.*.bgc <$halo_file/$i/Carmen_${i}_z0p520_fof_b0p2.dpp.halos  > fff_$i
#makemock 0 0 0 0 0 0 0.02 0.225 /home/piscioja/SDSSPix/Maps/mask-cmass-dr11v2-N-Anderson.pix 0.6 0.52 < fff_$i > clf_galaxies_$i
#IDfib 0 1 55.0 ../Boss/Mock_Test/boss_mock_$i.sphere.galaxies clf_galaxies_$i > collided_clf_$i
#awk '{if($4==-1){print $1,$2,$3} else {print $1,$2,$4}}' <collided_clf_$i >boss_mock_with_clf_$i
#DDrppi 0.01 10 20 boss_mock_with_clf_$i boss_mock_with_clf_$i >DD_clf_$i
#DDrppi 0.01 10 20 boss_mock_with_clf_$i /hd0/Research/Clustering/Boss/Mock_Test/Randoms_sphere  >D_collided_R_sphere_$i
ND1=`wc -l boss_mock_$i.$tag.galaxies | cut -f1 -d' '`
#wprp 20 $ND1 $ND1 25000000 25000000 DD_clf_$i D_collided_R_sphere_$i D_collided_R_sphere_$i ../Boss/Mock_Test/RR_sphere 40 >wprp_clf_$i
