#!/bin/bash

icc -std=c99 jackknife_it.c measure_Boss_wp.c gridlink1D.c utils.c -o measure_Boss_wp  -lgsl -lgslcblas -Wall -Wextra -g -xhost -ipo -opt-prefetch
