#!/bin/bash


for i in {1,2,4,6}
do
    n=0
    while [ $n -lt 20 ]
    do
       ./lcs_par fileE.in fileF.in i
       n=$((n+1))
    done
       echo "Número de threads: " $i
        
done