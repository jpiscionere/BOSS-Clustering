Ds=/hd0/Research/Clustering/Boss/DR3/dr3_Ds_zrange_trim
Rs=/hd0/Research/Clustering/Boss/DR3/dr3_Rs_zrange_trim
Di=/hd0/Research/Clustering/Boss/DR3/dr3_Di_mag_extinction
Ri=/hd0/Research/Clustering/Boss/DR3/dr3_Ri_short
mask=/hd0/Research/Clustering/Boss/DR3/window.sample14.trimmed.ply

Ds=dr3_Ds_zrange_trim
Di=dr3_Di_mag_extinction
mask=window.sample14.trimmed.ply 

nbins=17
area_tot=0.9076429
r_min=0.005
n_jack=150


#~/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_openmp_mag_cut $Ds $Di $mask outfile_DsDi $r_min 10 $nbins $n_jack 1 >dr3_DsDi_jackknife
~/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_openmp_mag_cut $Ds $Ri $mask outfile_DsDi $r_min 10 $nbins $n_jack 2 >dr3_DsRi_jackknife
#~/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_openmp_mag_cut $Rs $Di $mask outfile_DsDi $r_min 10 $nbins $n_jack 1 >dr3_RsDi_jackknife
#~/Clustering/Boss/Source/Vpac_Codes/measure_Boss_wp_openmp_mag_cut $Rs $Ri $mask outfile_DsDi $r_min 10 $nbins $n_jack 2 >dr3_RsRi_jackknife
