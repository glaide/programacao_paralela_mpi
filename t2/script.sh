#!/bin/bash

echo "Número de processos 2"
	for i in {1,2,3,4,5}
	do
	    mpirun -np 2 ./lcs_mpi 5t.txt 6t.txt
	done


echo "Número de processos 5"
	for i in {1,2,3,4,5}
	do
	    mpirun -np 4 ./lcs_mpi 5t.txt 6t.txt
	done