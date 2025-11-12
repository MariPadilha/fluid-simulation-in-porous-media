rm *.mod # deleta arquivos módulos
rm *.out # deleta executável
sh cleaning.sh  # apaga arquivos que não serão utilizados

#gfortran -O3 -fopenmp \
#/opt/nvidia/hpc_sdk/Linux_x86_64/24.9/compilers/bin/pgf90 -O3 -pg \
#/opt/nvidia/hpc_sdk/Linux_x86_64/24.9/compilers/bin/pgf90 -O3 -mp -pg \
#/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/compilers/bin/pgf90 -O3 -stdpar=gpu -acc -Minfo=accel \
#/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/compilers/bin/pgf90 -O3 -stdpar=multicore -acc -Minfo=accel \
/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/compilers/bin/pgf90 -O3 -stdpar=multicore -acc=multicore -fast -Minfo=all \
	comum.f90 \
properties_CH4.f90 \
       initial.f90 \
nonsymetric_mesh.f90 \
     comp_mean.f90 \
         main.f90  \
     equations.f90 \
   convergence.f90 \
     transient.f90 \
        output.f90 \
         probe.f90 \
      flametip.f90 \
-o cylinder_solver.out

#boundary.f90 \

#export OMP_NUM_THREADS=8
#nohup time ./cylinder_solver.out | tee archivelog.log

#sh posv.sh

