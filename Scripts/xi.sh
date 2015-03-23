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
redshift=0.186
#halo_file=/net/bender/data0/LasDamas/Carmen
halo_file=/ssd1/Research/halo_files/Carmen/Fof
map_dir=/home/piscioja/SDSSPix/Maps
map_full=$map_dir/mask-cmass-dr11v2-N-Anderson.pix
polygon_file=/hd0/Research/Clustering/Boss/dr11/dr11v2/mask-cmass-dr11v2-N-Anderson.ply
zmin=0.02
zmax=0.225
zspace=0
trim=sphere



for i in $(seq 2020 1 2029)
do

halobias_fof_nfw 3 4 1 $logMmin $siglogM $logM0 $logM1 $alpha 1 $gamma $fgal -21 $delta_vir $Mstar $redshift $i trashfile.out  < $halo_file/Carmen_${i}_z0p132_fof_b0p2.dpp.halos >output_$trim.$i
covar3 0.01 42.3 13 1000 0 1000 1 output_$trim.$i f 0 1 auto > xi_lowz_${i}.out
wp2 11.0 < xi_lowz_${i}.out >wp2_$i.out 
rm -r output_$trim.$i
done
