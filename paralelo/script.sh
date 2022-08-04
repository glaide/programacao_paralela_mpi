#!/bin/bash
for i in {2,1}
do

    ./lcs_par fileE.in fileF.in i

done
n=2
while [ $n -lt 3 ]
do 
    while [ $j -lt 20 ]
    do
    ./lcs_par fileE.in fileF.in $n
    j=$((j+1))

    done
  n=$((n+1))
done