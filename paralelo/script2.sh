#!/bin/bash

n=0
    while [ $n -lt 20 ]
    do
       ./lcs_par fileE.in fileF.in 6
       n=$((n+1))
    done
