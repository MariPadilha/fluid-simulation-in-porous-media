#!/bin/bash
#sh run.sh
#for j in 1 2 3 4 5 6 7 8 9 10
#do
  for i in 12 10 8 6 4 2 1
  do
  echo "Executando com $i threads"
	export OMP_NUM_THREADS=$i
	time ./cylinder_solver.out 2>> erro16.10.txt 1>> saida16.10.txt
  done
#done