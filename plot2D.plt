set terminal postscript eps enhanced color
set output filename.'.eps'

stats './dat/'.filename.'.dat' nooutput

set xlabel "x"
set xlabel "y"
set ylabel "T"
#set yrange [0:1.1]

set pm3d at bs
set palette
set hidden3d
set view 90, 0, 1, 1

set zrange [0.0:1.1]

splot './dat/'.filename.'.dat' i (20) using 1:2:3 with pm3d
