#! /bin/bash



#nden[1]=0.0002569
#nden[2]=0.0003243
#nden[3]=0.0002327
#nden[4]=0.0000942
#nden[5]=0.0001848



nden[1]=0.0002494737166611
nden[2]=0.0003150319712696
nden[3]=0.0002238406773599
nden[4]=0.0000912964784524
nden[5]=0.0001847687542961

i=1

for bin in bin1 bin2 bin3 bin4 bin_all 
do

nden1=${nden[$i]}
echo $nden1

paste ${bin}_DsDi_nojack.out ${bin}_DsRi_nojack.out ${bin}_RsDi_nojack.out ${bin}_RsRi_nojack.out >tmp 
awk '{wp=$2/$5-$8/$11; {print $1,wp,wp/'$nden1'}}'<tmp >wp_${bin}_nojack.out

i=`expr $i + 1`

done



