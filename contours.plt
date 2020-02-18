#! /usr/bin/gnuplot
   set termoption noenhanced
   set term postscript enhanced color
   set output 'output.ps'
   set pm3d
   set pm3d map
   set pm3d interpolate 5,5
   set palette rgbformulae 22,13,-31
   set contour
   set cntrparam levels 0
   set style line 2 lc rgb "#000000"
   set style line 3 lc rgb "#000000"
   set style increment user
   set size square
   set xr [ -1 : 1 ]
   set yr [ -1 : 1 ]
   splot 'IntpltdFx.dat' u 1:2:3 notit
   splot 'IntpltdFr.dat' u 1:2:3 notit
   splot 'IntpltdFt.dat' u 1:2:3 notit

   unset pm3d
   set contour
   set nosurface
   set cntrparam levels 0
   set style line 2 lc rgb "#ff0000"
   set style line 4 lc rgb "#0000ff"
   set style increment user
   set size square
   set xr [ -1 : 1 ]
   set yr [ -1 : 1 ]
   splot 'IntpltdFr.dat' u 1:2:3 w l notit, \
         'IntpltdFt.dat' u 1:2:3 w l notit
