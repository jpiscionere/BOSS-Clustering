#!/bin/bash
gcc -O3 -Wall `gsl-config --cflags` calculate_wp.c `gsl-config --libs` -lm  -o calculate_wp
icc -std=c99 measure_Boss_wp_openmp_mag_cut.c jackknife_it.c gridlink1D.c utils.c ../utils/progressbar.c -I../utils/ -o measure_Boss_wp_openmp_mag_cut  -lgsl -lgslcblas -Wall -Wextra -g -xhost -ipo -opt-prefetch -openmp
icc -std=c99 measure_Boss_wp_nojack_openmp.c gridlink1D.c utils.c ../utils/progressbar.c -I../utils/ -o measure_Boss_wp_nojack_openmp  -lgsl -lgslcblas -Wall -Wextra -g -xhost -ipo -opt-prefetch -openmp
icc -std=c99 measure_Boss_wp_nojack.c gridlink1D.c utils.c -o measure_Boss_wp_nojack  -lgsl -lgslcblas -Wall -Wextra -g -xhost -ipo -opt-prefetch
icc -std=c99 jackknife_it.c measure_Boss_wp.c gridlink1D.c utils.c -o measure_Boss_wp  -lgsl -lgslcblas -Wall -Wextra -g -xhost -ipo -opt-prefetch
gcc -O3 -Wall `gsl-config --cflags` filter_polygon3.c `gsl-config --libs` -lm  -o filter_polygon3
icc -std=c99 jackknife_it_standalone.c gridlink1D.c utils.c -o jackknife_it_standalone  -lgsl -lgslcblas -Wall -Wextra -g -xhost -ipo -opt-prefetch
