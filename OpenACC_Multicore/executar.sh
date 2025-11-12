#!/bin/bash
#sh run.sh
#for j in 1 2 3 4 5 6 7 8 9 10
#do
#  for i in 1 2 4 6 8 10 12 14 16 18 20 22 24
for i in 1 2 4 6 8 10 12
do
	export ACC_NUM_CORES=$i
	time ./cylinder_solver.out 2>> erro51Time.txt 1>> saida51Time.txt 
#  done
done


