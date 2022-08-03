#!/bin/bash


# for i in {4,2,1}
# do
#     n=0
#        echo "Número de threads: " $i
#     while [ $n -lt 20 ]
#     do
#        ./lcs_par fileE.in fileF.in i
#        n=$((n+1))
#     done
        
# done
n=0
while [ $n -lt 20 ]
do
   ./lcs_par fileE.in fileF.in 
   n=$((n+1))
done