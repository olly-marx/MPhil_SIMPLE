set terminal gif animate size 800,600 delay 20
set output filename.'.gif'

stats './dat/'.filename.'.dat' nooutput

set xlabel "x"
set ylabel "y"
set zlabel "T"

set pm3d at bs
set palette
set hidden3d

set zrange [0.0:3.5]
#set xrange [145:155]
do for [i=1:int(STATS_blocks)] {
	set view (75), (10), 1, 1
	splot './dat/'.filename.'.dat' index (i-1) using 1:2:3 with lines
}
