set terminal postscript eps enhanced color

set output './results/'.filename.'.eps'

set size ratio 1/1

set xrange [0:0.1]
set xtics 0,0.05,0.1
set yrange [0:0.1]
set ytics 0,0.05,0.1

set ticslevel 0.5

set view map

set pm3d at bs
set style fill transparent solid 0.8 noborder
set hidden3d

set palette rgb 33,13,10

set parametric
set isosamples 51, 51

set title "Lid driven cavity p-v fields after ".filename." iterations"

set xlabel '{/:Italic x}' font ",12"
set ylabel '{/:Italic y}' font ",12"
set cblabel '{/:Italic p}' font ",12"
set cbrange [-0.05:0.06]
plot './dat/'.filename.'.dat' i (1) u 1:2:6 w image notitle, \
	'./dat/'.filename.'.dat' i (1) u 1:2:($3/50):($4/50) every 2:2 w vec \
	lc -1 head filled size screen 0.01,30,45 notitle
