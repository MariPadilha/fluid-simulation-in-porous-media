rm -f *.mod  # remove arquivos .mod, sem erro se não existir
rm -f *.out  # remove arquivos .out
sh cleaning.sh  # executa script para limpeza extra

NVCC=nvcc
ARCH=sm_86
OPT="-O2"
STD="-std=c++11"

gcc -c probe.c -o probe.o
gcc -c convergence.c -o convergence.o 
gcc -c flametip.c -o flametip.o

$NVCC -arch=$ARCH $OPT $STD -lineinfo \
  xm_ym.cu \
  x_y.cu \
  vol.cu \
  solve_C.cu \
  solve_P.cu \
  solve_U.cu \
  solve_V.cu \
  solve_Z.cu \
  upwind_U_pair.cu \
  upwind_V_pair.cu \
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

