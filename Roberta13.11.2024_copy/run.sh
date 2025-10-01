rm -f *.mod  # remove arquivos .mod, sem erro se não existir
rm -f *.out  # remove arquivos .out
sh cleaning.sh  # executa script para limpeza extra


#gfortran -O3 -fopenmp \
#/opt/nvidia/hpc_sdk/Linux_x86_64/24.9/compilers/bin/pgf90 -O3 -pg \
#/opt/nvidia/hpc_sdk/Linux_x86_64/24.9/compilers/bin/pgf90 -O3 -mp -pg \

#/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/compilers/bin/pgf90 
gcc -O2 -I. \
  comum.c \
  max.c \
  boundary.c \
  initial.c \
  nonsymetric_mesh.c \
  comp_mean.c \
  main.c  \
  equations.c \
  convergence.c \
  transient.c \
  output.c \
  probe.c \
  flametip.c \
  -o cylinder_solver.out -lm



#export OMP_NUM_THREADS=8
#nohup time ./cylinder_solver.out | tee archivelog.log

#sh posv.sh

