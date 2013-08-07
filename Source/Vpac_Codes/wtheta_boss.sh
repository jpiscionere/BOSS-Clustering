#!/bin/bash
for bin in bin1 bin2 bin3 bin4
do
DDtheta -3 -1 20 /home/piscioja/Clustering/Boss/Data/spectro_sample_$bin.txt /home/piscioja/Clustering/Boss/Data/spectro_sample_$bin.txt >${bin}_DD_theta.out &
DDtheta -3 -1 20 /home/piscioja/Clustering/Boss/Data/spectro_sample_$bin.txt /home/piscioja/Clustering/Boss/Data/randoms_spectro_$bin.out >${bin}_DR_theta.out &
DDtheta -3 -1 20 /home/piscioja/Clustering/Boss/Data/randoms_spectro_$bin.out /home/piscioja/Clustering/Boss/Data/randoms_spectro_$bin.out >${bin}_RR_theta.out &
done 
