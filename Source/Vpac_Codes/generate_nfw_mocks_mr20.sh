#! /bin/bash

output=stomp_sdssmock_main20_esmeralda
maxtheta=0.1
binning=10
#fof=/home/piscioja/halobias/halobias_fof_nfw
fof=halobias_fof_nfw_quiet
so=/home/piscioja/halobias/halobias_so_nfw2
map_dir=/home/piscioja/SDSSPix/Maps
map_sdss=$map_dir/window.dr72brvoid0.stripe_trim.pix
map_lss=$map_dir/lss_geometry.dr72.stripe_trim.pix
map_full=$map_dir/lss_geometry.fullsphere.pix
map_overlap1=$map_dir/window.dr72brvoid1.overlap.balkanized.noweight.pix
map_overlap2=$map_dir/window.dr72brvoid1.overlap.noweight.pix
stompdirectory=/home/piscioja/astrostomp/astro-stomp-read-only/examples
logMmin=11.8067544233
siglogM=0.0963472124455
logM0=12.6369512043
logM1=13.1533827128
alpha=1.06123071572
delta_vir=377 #linking length of 0.2
Mstar=2.29E12
redshift=0.082


for i in $(seq 3020 1 3029) 
do


#for gamma in 0.206254091967 
#do

#for fgal in 1.79335961349 
#do



time $fof 3 4 1 $logMmin $siglogM $logM0 $logM1 $alpha 1 $gamma $fgal -20 $delta_vir $Mstar $redshift  trashfile.out $i < /ssd1/Research/halo_files/Esmeralda/Esmeralda_${i}_z0p082_fof_b0p2.dpp.halos > halobias_fof_nfw_${i}_fff

time /home/piscioja/halobias/makemock 1 1 0 0 0 0 0.02 0.106 $map_lss 0.6 0 < halobias_fof_nfw_${i}_fff > halobias_fof_nfw_${i}_galaxies

time $stompdirectory/stomp_galaxy_autocorrelation --map_file=$map_lss --galaxy_file=halobias_fof_nfw_${i}_galaxies  -output_tag=${output}.$i.nfw.fof.lss_geometry2  --theta_max=$maxtheta --n_bins_per_decade=$binning -single_index 

time $stompdirectory/stomp_galaxy_autocorrelation --map_file=$map_overlap1 --galaxy_file=halobias_fof_nfw_${i}_galaxies  -output_tag=${output}.$i.nfw.fof.overlap2  --theta_max=$maxtheta --n_bins_per_decade=$binning -single_index

rm halobias*


done
#done
#done



