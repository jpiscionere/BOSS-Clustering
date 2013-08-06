#! /bin/bash

~/Clustering/Boss/calculate_wp bin1_DsDi_with_norm.out.new_jack bin1_DsRi_with_norm.out.new_jack bin1_RsDi_with_norm.out.new_jack bin1_RsRi_with_norm.out.new_jack 20 100 0.0002897403 1  bin1_covar_inv_with_norm.dr11.out.new_jack >bin1_wp_with_norm.new_jackknife
~/Clustering/Boss/calculate_wp bin2_DsDi_with_norm.out.new_jack bin2_DsRi_with_norm.out.new_jack bin2_RsDi_with_norm.out.new_jack bin2_RsRi_with_norm.out.new_jack 20 100 0.0003750743 1  bin2_covar_inv_with_norm.dr11.out.new_jack >bin2_wp_with_norm.new_jackknife
~/Clustering/Boss/calculate_wp bin3_DsDi_with_norm.out.new_jack bin3_DsRi_with_norm.out.new_jack bin3_RsDi_with_norm.out.new_jack bin3_RsRi_with_norm.out.new_jack 20 100 0.0002662834 1 bin3_covar_inv_with_norm.dr11.out.new_jack >bin3_wp_with_norm.new_jackknife
~/Clustering/Boss/calculate_wp bin4_DsDi_with_norm.out.new_jack bin4_DsRi_with_norm.out.new_jack bin4_RsDi_with_norm.out.new_jack bin4_RsRi_with_norm.out.new_jack 20 100 0.0001091332 1 bin4_covar_inv_with_norm.dr11.out.new_jack >bin4_wp_with_norm.new_jackknife
