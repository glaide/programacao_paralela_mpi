#!/bin/bash
make
for i in {1,2,4,6}
do
    ./lcs_par fileE.in fileF.in i
done