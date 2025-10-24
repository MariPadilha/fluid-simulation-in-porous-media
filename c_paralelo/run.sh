rm -f *.mod  # remove arquivos .mod, sem erro se não existir
rm -f *.out  # remove arquivos .out
sh cleaning.sh  # executa script para limpeza extra


#gfortran -O3 -fopenmp \
#/opt/nvidia/hpc_sdk/Linux_x86_64/24.9/compilers/bin/pgf90 -O3 -pg \
#/opt/nvidia/hpc_sdk/Linux_x86_64/24.9/compilers/bin/pgf90 -O3 -mp -pg \
gcc -c probe.c -o probe.o
gcc -c convergence.c -o convergence.o 
gcc -c flametip.c -o flametip.o

#/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/compilers/bin/pgf90 
nvcc -arch=sm_86 -O2 -std=c++11 -I. -diag-suppress 1650\
  xm_ym.cu \
  x_y.cu \
  vol.cu \
  solve_C.cu \
  solve_P.cu \
  solve_U.cu \
  solve_V.cu \
  solve_Z.cu \
  upwind_Vj.cu \
  upwind_Vi.cu \
  upwind_Uj.cu \
  upwind_Ui.cu \
  transient.c \
  resz.cu \
  resv.cu \
  resu.cu \
  resc.cu \
  nonsymetric_mesh.cu \
  max_reduce.cu \
  main.cu \
  initial.cu \
  grids.cu \
  dx_dy.cu \
  comum.cu \
  comp_mean.cu \
  bcZ.cu \
  bcUV.cu \
  bcP.cu \
  bcC.cu \
  area_das_faces.cu \
  output.cu \
  probe.o \
  convergence.o \
  flametip.o \
  -o cylinder_solver.out -lm



#export OMP_NUM_THREADS=8
#nohup time ./cylinder_solver.out | tee archivelog.log

#sh posv.sh

