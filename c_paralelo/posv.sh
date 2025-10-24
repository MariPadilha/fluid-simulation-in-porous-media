gfortran \
         comum.f90 \
properties_CH4.f90 \
   pos_process.f90 \
-o pos_process.out

./pos_process.out

#python2.7 -W ignore analytical.py
python3 -W ignore graphicsv.py
python3 -W ignore malha_teste.py
# ffmpeg -framerate 6 -start_number 102 -i output/streamlines%d.png output/s_animation.mp4
#xdg-open output/animation2.mp4

