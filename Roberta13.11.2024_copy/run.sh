rm -f *.mod  # remove arquivos .mod, sem erro se não existir
rm -f *.out  # remove arquivos .out
sh cleaning.sh  # executa script para limpeza extra


#gfortran -O3 -fopenmp \
#/opt/nvidia/hpc_sdk/Linux_x86_64/24.9/compilers/bin/pgf90 -O3 -pg \
#/opt/nvidia/hpc_sdk/Linux_x86_64/24.9/compilers/bin/pgf90 -O3 -mp -pg \

#/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/compilers/bin/pgf90 
nvcc -arch=sm_86 -O2 -std=c++11 -I. \
  xm_ym.cu \
  x_y.cu \
  upwind_Vj.cu \
  upwind_Vi.cu \
  upwind_Uj.cu \
  upwind_Ui.cu \
  transient.c \
  resz.cu \
  resv.cu \
  resu.cu \
  resc.cu \
  probe.c \
  output.c  \
  nonsymetric_mesh.cu \
  max_reduce.cu \
  main.cu \
  initial.cu \
  grids.cu \
  flametip.c \
  dx_dy.cu \
  convergence.c \
  comum.cu \
  comp_mean.c \
  bcZ.cu \
  bcUV.cu \
  bcP.cu \
  bcC.cu \
  area_das_faces.cu \
  -o cylinder_solver.out -lm



#export OMP_NUM_THREADS=8
#nohup time ./cylinder_solver.out | tee archivelog.log

#sh posv.sh

