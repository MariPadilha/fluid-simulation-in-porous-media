gfortran \
         comum.f90 \
properties_CH4.f90 \
   pos_process.f90 \
-o pos_process.out

./pos_process.out

#python2.7 -W ignore analytical.py
python3 -W ignore graphicsT.py
#ffmpeg -framerate 12 -start_number 102 -i output/temperature%d.png output/animation_T.mp4
#xdg-open output/animation2.mp4
