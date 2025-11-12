#!/bin/bash
#sh run.sh
#for j in 1 2 3 4 5 6 7 8 9 10
#do
#  for i in 1 2 4 6 8 10 12 14 16 18 20 22 24
for i in 1 2 3 4 5 6 7 8 
do
#	export OMP_NUM_THREADS=$i
	#time ./cylinder_solver.out 2>> erro204Time.txt 1>> saida204Time.txt 
	time ./cylinder_solver.out 2>> erro153Time.txt 1>> saida153Time.txt 
#  done
done


