#!/bin/bash

for i in $(seq 2021 1 2031)
do



        ./make_boss_mock.sh 0 $i 0 &

done

