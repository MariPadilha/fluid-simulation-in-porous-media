gfortran \
         comum.f90 \
properties_CH4.f90 \
   pos_process.f90 \
-o pos_process.out

./pos_process.out

#python2.7 -W ignore analytical.py
python3 -W ignore graphicsC.py
ffmpeg -framerate 6 -start_number 102 -i output/concentration%d.png output/animation_C.mp4
#xdg-open output/animation2.mp4
